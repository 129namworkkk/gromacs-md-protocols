#!/usr/bin/env bash
# =============================================================================
# prepare_system.sh — PDB → simulation-ready system, with automatic repair
#
#   1. split ligand (holo) / remove ligand (apo)
#   2. fix_structure.py: model missing residues, rebuild missing atoms,
#      replace non-standard residues, remove heterogens   (PDBFixer)
#   3. pdb2gmx (non-interactive)  — if it still fails on missing atoms,
#      the structure is repaired again and pdb2gmx is retried
#   4. holo only: merge the ligand (LigParGen lig.gro + lig.itp) into
#      complex.gro / topol.top, add ligand position restraints
#   5. box → solvate → neutralise + 0.15 M NaCl → index.ndx
#
# Output (in --workdir): solv_ions.gro, topol.top, *.itp, index.ndx,
# fix_report.json, prepare.log — ready for runpod/run_md.sh or the guides' EM step.
#
# Examples
#   bash scripts/prepare_system.sh --system pp   --pdb 3bik.pdb          --workdir runs/3bik
#   bash scripts/prepare_system.sh --system apo  --pdb complex.pdb       --workdir runs/apo
#   bash scripts/prepare_system.sh --system holo --pdb complex.pdb --ligand UNL \
#        --lig-gro lig.gro --lig-itp lig.itp --workdir runs/holo
#   ... add --fasta target.fasta when the PDB has no SEQRES records (docking / AlphaFold / Chimera output)
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=../runpod/common.sh
source "$REPO_DIR/runpod/common.sh"

usage() {
cat <<'EOF'
Usage: prepare_system.sh --system {pp|apo|holo} --pdb FILE --workdir DIR [options]

Structure repair
  --fasta FILE       full sequence per chain (needed to model gaps when the PDB has no SEQRES)
  --max-gap N        model internal gaps up to N residues (default 15); longer gaps are split (TER)
  --add-terminal     also model missing terminal tails (default: no)
  --no-fix           skip PDBFixer (pdb2gmx errors on missing atoms still trigger a repair)
  --strict           stop if chain breaks remain after repair (instead of splitting with TER)

Force field (see docs/FORCE_FIELDS.md)
  --ff NAME          pdb2gmx force field (default: gromos45a3 for pp, oplsaa for apo/holo)
  --water NAME       water model (default: spc for pp, tip3p for apo/holo)

Ligand (holo / apo)
  --ligand RES       ligand residue name in the PDB (default: UNL)
  --lig-gro FILE     ligand coordinates from LigParGen (.gro)   } holo: if missing, the script
  --lig-itp FILE     ligand topology from LigParGen (.itp)      } stops after the protein topology

Box / ions
  --box TYPE         cubic | triclinic | dodecahedron (default: cubic for pp, triclinic otherwise)
  --distance NM      solute–box distance (default 1.0)
  --conc M           NaCl concentration (default 0.15)
EOF
}

SYSTEM=""; PDB=""; WORKDIR=""; FASTA=""; MAXGAP=15; ADDTERM=0; NOFIX=0; STRICT=0
FF=""; WATER=""; LIG="UNL"; LIGGRO=""; LIGITP=""; BOX=""; DIST=1.0; CONC=0.15
while [[ $# -gt 0 ]]; do
    case "$1" in
        --system) SYSTEM="$2"; shift 2 ;;       --pdb) PDB="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;     --fasta) FASTA="$2"; shift 2 ;;
        --max-gap) MAXGAP="$2"; shift 2 ;;      --add-terminal) ADDTERM=1; shift ;;
        --no-fix) NOFIX=1; shift ;;             --strict) STRICT=1; shift ;;
        --ff) FF="$2"; shift 2 ;;               --water) WATER="$2"; shift 2 ;;
        --ligand) LIG="$2"; shift 2 ;;          --lig-gro) LIGGRO="$2"; shift 2 ;;
        --lig-itp) LIGITP="$2"; shift 2 ;;      --box) BOX="$2"; shift 2 ;;
        --distance) DIST="$2"; shift 2 ;;       --conc) CONC="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) usage; die "unknown option: $1" ;;
    esac
done
[[ -n "$SYSTEM" && -n "$PDB" && -n "$WORKDIR" ]] || { usage; die "--system, --pdb and --workdir are required"; }
[[ -f "$PDB" ]] || die "PDB not found: $PDB"
case "$SYSTEM" in
    pp)       : "${FF:=gromos45a3}"; : "${WATER:=spc}";   : "${BOX:=cubic}";     MDPDIR="$REPO_DIR/protein-protein/mdp" ;;
    apo|holo) : "${FF:=oplsaa}";     : "${WATER:=tip3p}"; : "${BOX:=triclinic}"; MDPDIR="$REPO_DIR/protein-ligand/$SYSTEM/mdp" ;;
    *) die "--system must be pp, apo or holo" ;;
esac
IONS_MDP="$MDPDIR/ions.mdp"
PY="${PYTHON:-python3}"
FIXER="$SCRIPT_DIR/fix_structure.py"

mkdir -p "$WORKDIR"
PDB="$(cd "$(dirname "$PDB")" && pwd)/$(basename "$PDB")"
[[ -n "$FASTA" ]] && FASTA="$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")"
[[ -n "$LIGGRO" ]] && LIGGRO="$(cd "$(dirname "$LIGGRO")" && pwd)/$(basename "$LIGGRO")"
[[ -n "$LIGITP" ]] && LIGITP="$(cd "$(dirname "$LIGITP")" && pwd)/$(basename "$LIGITP")"
cd "$WORKDIR"
exec > >(tee -a prepare.log) 2>&1
load_gmx
export GMX_MAXBACKUP=-1
log "system=$SYSTEM  pdb=$PDB  ff=$FF  water=$WATER  box=$BOX  d=$DIST nm  NaCl=$CONC M"

# ---------------------------------------------------------------- 1. ligand split
cp "$PDB" input.pdb
if grep -qE "^(ATOM  |HETATM).{11}$LIG" input.pdb; then
    grep -E "^(ATOM  |HETATM).{11}$LIG" input.pdb > LIG.pdb
    grep -vE "^(ATOM  |HETATM).{11}$LIG" input.pdb > protein_raw.pdb
    log "ligand $LIG: $(grep -c . LIG.pdb) atoms → LIG.pdb (removed from the protein)"
else
    cp input.pdb protein_raw.pdb
    [[ "$SYSTEM" == "holo" ]] && log "WARNING: no residue named $LIG in the PDB (use --ligand RES)"
fi

# ---------------------------------------------------------------- 2. repair
run_fixer() {   # extra args...
    local args=(protein_raw.pdb -o protein_fixed.pdb --report fix_report.json --max-gap "$MAXGAP")
    [[ -n "$FASTA" ]] && args+=(--fasta "$FASTA")
    [[ $ADDTERM -eq 1 ]] && args+=(--add-terminal)
    local rc=0
    "$PY" "$FIXER" "${args[@]}" "$@" || rc=$?
    if [[ $rc -eq 1 ]]; then
        die "structure repair failed (PDBFixer missing? conda install -c conda-forge pdbfixer openmm; or set PYTHON=/path/to/python)"
    fi
    if [[ $rc -eq 3 && $STRICT -eq 1 ]]; then
        die "chain breaks remain (see above). Provide --fasta, raise --max-gap, or build the loop with MODELLER/AlphaFold."
    fi
    return 0
}
STRUCT=protein_raw.pdb
if [[ $NOFIX -eq 0 ]]; then
    log "repairing structure (PDBFixer)"
    run_fixer
    STRUCT=protein_fixed.pdb
else
    "$PY" "$FIXER" protein_raw.pdb --check-only || true
    grep -vE '^HETATM' protein_raw.pdb > protein_nohet.pdb; STRUCT=protein_nohet.pdb
fi

# ---------------------------------------------------------------- 3. pdb2gmx (+ automatic retry)
pdb2gmx_run() {
    gmx pdb2gmx -f "$1" -o processed.gro -p topol.top -ff "$FF" -water "$WATER" -ignh > logs_pdb2gmx.txt 2>&1
}
if ! pdb2gmx_run "$STRUCT"; then
    err=$(grep -A6 -E "Fatal error" logs_pdb2gmx.txt | tr '\n' ' ' || true)
    log "pdb2gmx failed: ${err:0:400}"
    if grep -qE "not found in the input file|is missing in the input|not found in residue topology database|Atom .* not found" logs_pdb2gmx.txt \
       && [[ "$STRUCT" != "protein_fixed.pdb" ]]; then
        log "→ repairing the structure (PDBFixer) and retrying pdb2gmx"
        run_fixer
        STRUCT=protein_fixed.pdb
        pdb2gmx_run "$STRUCT" || {
            tail -n 25 logs_pdb2gmx.txt
            die "pdb2gmx still fails after repair. Common causes: a ligand/ion still in the file (use --ligand), a residue the force field does not know, or a wrong --ff/--water pair (docs/FORCE_FIELDS.md)."
        }
    else
        tail -n 25 logs_pdb2gmx.txt
        [[ "$STRUCT" == "protein_fixed.pdb" ]] && log "the structure was already repaired; see docs/TROUBLESHOOTING.md and docs/FORCE_FIELDS.md"
        tail -n 25 logs_pdb2gmx.txt
        die "pdb2gmx failed (see logs_pdb2gmx.txt)"
    fi
fi
log "pdb2gmx OK → processed.gro, topol.top ($(grep -c 'moleculetype' topol.top *.itp 2>/dev/null | awk -F: '{s+=$2} END{print s}') molecule types)"
COORD=processed.gro

# ---------------------------------------------------------------- 4. ligand merge (holo)
if [[ "$SYSTEM" == "holo" ]]; then
    if [[ -z "$LIGGRO" || -z "$LIGITP" ]]; then
        log "Protein topology ready. Next: make lig.gro + lig.itp for LIG.pdb with LigParGen (add H in Chimera first),"
        log "then re-run this command with --lig-gro lig.gro --lig-itp lig.itp (PDBFixer/pdb2gmx are simply redone)."
        exit 0
    fi
    cp "$LIGITP" lig.itp; cp "$LIGGRO" lig.gro
    LIGMOL=$(awk '/^\[ *moleculetype *\]/{f=1;next} f && $0!~/^[ \t]*;/ && NF>=1 {print $1; exit}' lig.itp)
    [[ -n "$LIGMOL" ]] || die "cannot read the [ moleculetype ] name from lig.itp"
    NLIG=$(sed -n '2p' lig.gro | tr -d ' ')
    LIGRES=$(sed -n '3p' lig.gro | cut -c6-10 | tr -d ' ')
    NPRO=$(sed -n '2p' processed.gro | tr -d ' ')
    # complex.gro = protein atoms + ligand atoms + protein box
    {
        echo "Protein-ligand complex"
        echo "$((NPRO + NLIG))"
        sed -n "3,$((NPRO + 2))p" processed.gro
        sed -n "3,$((NLIG + 2))p" lig.gro
        tail -n 1 processed.gro
    } > complex.gro
    # ligand position restraints on heavy atoms
    printf '0 & ! a H*\nq\n' | gmx make_ndx -f lig.gro -o index_lig.ndx > logs_posre_lig.txt 2>&1
    NG=$(grep -c '^\[' index_lig.ndx)
    echo $((NG - 1)) | gmx genrestr -f lig.gro -n index_lig.ndx -o posre_lig.itp -fc 1000 1000 1000 >> logs_posre_lig.txt 2>&1 \
        || die "genrestr failed (logs_posre_lig.txt)"
    # topol.top: include ligand right after the force field, add it to [ molecules ]
    "$PY" - "$LIGMOL" <<'PYEOF'
import re, sys
lig = sys.argv[1]
top = open("topol.top").read().splitlines()
out, done = [], False
for line in top:
    out.append(line)
    if not done and re.match(r'\s*#include\s+".*forcefield\.itp"', line):
        out += ["", "; Include ligand topology", '#include "lig.itp"', "",
                "; Ligand position restraints", "#ifdef POSRES", '#include "posre_lig.itp"', "#endif"]
        done = True
if not done:
    sys.exit("forcefield.itp include not found in topol.top")
text = "\n".join(out).rstrip() + "\n"
if not re.search(r'^\s*' + re.escape(lig) + r'\s+1\s*$', text, re.M):
    text += f"{lig:<20}1\n"
open("topol.top", "w").write(text)
PYEOF
    # coordinates sanity check: LigParGen must keep the docked pose
    if [[ -s LIG.pdb ]]; then
        "$PY" - <<'PYEOF' || true
import math
def cen_pdb(f):
    xs=[(float(l[30:38]),float(l[38:46]),float(l[46:54])) for l in open(f) if l[:6] in ("ATOM  ","HETATM")]
    return [sum(c)/len(xs)/10 for c in zip(*xs)]
def cen_gro(f):
    L=open(f).read().splitlines(); n=int(L[1])
    xs=[(float(l[20:28]),float(l[28:36]),float(l[36:44])) for l in L[2:2+n]]
    return [sum(c)/len(xs) for c in zip(*xs)]
d=math.dist(cen_pdb("LIG.pdb"), cen_gro("lig.gro"))
print(f"[check] ligand centre shift docked pose → lig.gro: {d:.2f} nm" + ("  ⚠ ligand moved! check that LigParGen kept the coordinates" if d>0.2 else "  (OK)"))
PYEOF
    fi
    COORD=complex.gro
    log "ligand $LIGMOL merged: complex.gro ($((NPRO + NLIG)) atoms), lig.itp + posre_lig.itp included in topol.top"
fi

# ---------------------------------------------------------------- 5. box, water, ions, index
gmx editconf -f "$COORD" -o newbox.gro -bt "$BOX" -d "$DIST" -c > logs_editconf.txt 2>&1
case "$WATER" in tip4p*) WBOX=tip4p.gro ;; tip5p*) WBOX=tip5p.gro ;; *) WBOX=spc216.gro ;; esac
gmx solvate -cp newbox.gro -cs "$WBOX" -p topol.top -o solv.gro > logs_solvate.txt 2>&1
gmx grompp -f "$IONS_MDP" -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 > logs_ions_grompp.txt 2>&1 \
    || { tail -n 30 logs_ions_grompp.txt; die "grompp (ions) failed — see logs_ions_grompp.txt"; }
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc "$CONC" \
    > logs_genion.txt 2>&1 || { tail -n 20 logs_genion.txt; die "genion failed"; }

if [[ "$SYSTEM" == "holo" ]]; then
    printf '"Protein" | r %s\nname %s Protein_UNL\nq\n' "$LIGRES" \
        "$(echo q | gmx make_ndx -f solv_ions.gro -o /dev/null 2>/dev/null | grep -cE '^ *[0-9]+ ')" \
        | gmx make_ndx -f solv_ions.gro -o index.ndx > logs_index.txt 2>&1
    grep -q "Protein_UNL" index.ndx || die "could not create the Protein_UNL index group (logs_index.txt)"
else
    echo q | gmx make_ndx -f solv_ions.gro -o index.ndx > logs_index.txt 2>&1
fi

log "molecules:"; sed -n '/^\[ *molecules *\]/,$p' topol.top | sed 's/^/    /'
log "READY: $(pwd)/solv_ions.gro + topol.top + index.ndx"
log "next:  bash $REPO_DIR/runpod/run_md.sh --system $SYSTEM --workdir $(pwd) --quick-test"
