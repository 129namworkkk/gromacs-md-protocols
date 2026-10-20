# Troubleshooting

| Symptom | Usual cause | Fix |
|---------|-------------|-----|
| `pdb2gmx`: *atom CG used in that entry is not found* | residue with missing atoms | **automatic:** `scripts/prepare_system.sh` · manual: [PDBFixer](#1-missing-atoms--residues-pdb2gmx) |
| EM blows up / LINCS warnings at one spot, a very long bond | **missing residues** (gap) — pdb2gmx bonds the two ends of the gap directly | **automatic:** `scripts/prepare_system.sh` (models the gap or splits the chain with TER) |
| `Residue 'XXX' not found in residue topology database` | a non-standard residue (MSE, HOH, ligand, ion) is still in the file | `scripts/fix_structure.py` replaces / removes it; split the ligand off (`--ligand`) |
| `grompp`: group `Water_and_ions` / `Protein_UNL` not found | the index has no group matching `tc-grps` | [Recreate the groups](#2-tc-grps-group-errors) |
| `trjconv`: atom-count mismatch | output group *System* chosen while the `.xtc` contains only Protein | choose an output group that matches the `.xtc` (`gmx check -f file.xtc`) |
| Huge RMSD spikes while SASA does not change | PBC artefact | run `-pbc nojump`, then `-pbc mol -center -ur compact` (see [ANALYSIS.md](../protein-ligand/ANALYSIS.md)) |
| `gmx_MMPBSA`: *ChamberParm does not support all potential terms…* | the GROMACS topology cannot be converted to Amber | [Build an Amber topology for the ligand](#3-gmx_mmpbsa-chamberparm-error) |
| `gmx_MMPBSA`: `AMBERHOME` error | AmberTools not activated | `conda activate gmxMMPBSA` and run again |

---

## 1. Missing atoms / residues (pdb2gmx)

### Automatic (recommended)

```bash
# check only (no extra dependencies)
python scripts/fix_structure.py pro.pdb --check-only

# repair + pdb2gmx + box + water + ions + index in one command
bash scripts/prepare_system.sh --system apo --pdb pro.pdb --workdir runs/apo
#   no SEQRES records (docking, AlphaFold, Chimera output)? add --fasta target.fasta
```

`prepare_system.sh` always repairs the structure before `pdb2gmx`; if `pdb2gmx` still reports missing atoms, the script repairs the structure again and retries `pdb2gmx` (also with `--no-fix`). A detailed report is written to `fix_report.json` and `prepare.log`.

| Situation | Result |
|-----------|--------|
| Internal gap ≤ `--max-gap` (default 15) and the sequence is known (SEQRES or `--fasta`) | missing residues modelled |
| Longer gap, or no sequence available | chain split with `TER` + warning (`--strict` stops instead of splitting) |
| Missing N-/C-terminal tail | skipped (usually disordered); `--add-terminal` models it |
| Truncated side chains, missing OXT | rebuilt |

> PDBFixer builds long loops (> 15 residues) with simple geometry only. For an important region, build it with MODELLER or graft it from an AlphaFold model before running `prepare_system.sh`.

### Manual

```text
Fatal error:
Residue 72 named LYS of a molecule in the input file was mapped
to an entry in the topology database, but the atom CG used in
that entry is not found in the input file. Perhaps your atom
and/or residue naming needs to be fixed.
```

```bash
conda install -c conda-forge pdbfixer openmm -y
pdbfixer pro.pdb --add-atoms=heavy --add-residues --output pro_fixed.pdb

# run again
gmx pdb2gmx -f pro_fixed.pdb -o processed.gro -p topol.top -ter -ff oplsaa -water spc
```

## 2. tc-grps group errors

Check whether `solv_ions.gro` contains ions. If `grompp` still reports a group error, recreate the groups with `make_ndx`:

```bash
gmx make_ndx -f em.gro -o index.ndx
```

```text
> 1 | 13                     # Protein | UNL   → Protein_UNL   (HOLO)
> r SOL | r NA | r CL        # → rename to Water_and_ions
> name XX Water_and_ions
> q
```

## 3. gmx_MMPBSA: ChamberParm error

Build an Amber (GAFF2) topology for the ligand:

```bash
antechamber -i LIG.pdb -fi pdb -o LIG.mol2 -fo mol2 -c bcc -nc YOUR_CHARGE -at gaff2
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod
```

Create `tleap.in`:

```bash
cat > tleap.in << 'EOF'
source leaprc.protein.ff14SB
source leaprc.gaff2

LIG = loadmol2 LIG.mol2
loadAmberParams LIG.frcmod

check LIG
saveOff LIG LIG.lib
saveAmberParm LIG LIG.prmtop LIG.inpcrd
quit
EOF

ls -l tleap.in LIG.mol2 LIG.frcmod     # check that the files exist
tleap -f tleap.in
```

It worked when all of these exist:

```bash
ls -l LIG.lib LIG.prmtop LIG.inpcrd
```
