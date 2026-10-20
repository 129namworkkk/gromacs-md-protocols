#!/usr/bin/env bash
# =============================================================================
# run_md.sh — resumable EM → NVT → NPT → production MD pipeline
#
# Runs the same steps as the guides in this repository, non-interactively,
# so a long simulation (100 ns or 500 ns) can run unattended on a RunPod GPU
# pod (or any Linux machine with GROMACS). Re-running the same command after
# an interruption skips finished stages and resumes from checkpoints.
#
# Usage:
#   bash runpod/run_md.sh --system holo --length 500 --workdir /workspace/runs/my_complex
#
# Run `bash runpod/run_md.sh --help` for all options.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
cat <<'EOF'
Usage: run_md.sh --system {pp|apo|holo} [options]

Required
  --system S        pp   = protein–protein   (protein-protein/mdp)
                    apo  = protein only      (protein-ligand/apo/mdp)
                    holo = protein + ligand  (protein-ligand/holo/mdp)

Common options
  --length N        production length: 100 (default) or 500 [ns]
  --workdir DIR     folder with the prepared system (default: current folder)
  --start FILE      solvated + ionised structure
                    (default: solv_ions.gro, or the only *_ion.pdb / *_ion.gro)
  --top FILE        topology (default: topol.top, or LM.top for pp)
  --index FILE      index file (default: index.ndx if present;
                    required for holo because of the Protein_UNL group)
  --stages "..."    stages to run, in order (default: "em nvt npt md")
  --grompp-only     for the LAST stage in --stages, only build the .tpr
                    (useful before benchmark.sh)

Advanced
  --mdp-dir DIR     folder with the .mdp files (default: the repo folder for --system)
  --md-mdp FILE     custom production .mdp (overrides --length)
  --ntomp N         OpenMP threads (default: cores available to the container)
  --gpu-id ID       GPU to use (default: let GROMACS choose)
  --maxwarn N       grompp -maxwarn for every stage
                    (default: 1 for pp = GROMOS warning; 1 for apo/holo NPT =
                    Berendsen warning; 0 otherwise. Always read logs/*.grompp.log)
  --quick-test      run every stage for only 500 steps in WORKDIR/quicktest
                    to check that inputs are correct before the real run
  -h, --help        show this help
EOF
}

# ------------------------------ arguments ------------------------------------
SYSTEM=""; LENGTH=100; WORKDIR="$PWD"; START=""; TOP=""; INDEX=""
MDP_DIR=""; MD_MDP=""; NTOMP=""; GPU_ID=""; MAXWARN=""
STAGES="em nvt npt md"; GROMPP_ONLY=0; QUICK=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --system)      SYSTEM="$2"; shift 2 ;;
        --length)      LENGTH="$2"; shift 2 ;;
        --workdir)     WORKDIR="$2"; shift 2 ;;
        --start)       START="$2"; shift 2 ;;
        --top)         TOP="$2"; shift 2 ;;
        --index)       INDEX="$2"; shift 2 ;;
        --stages)      STAGES="$2"; shift 2 ;;
        --grompp-only) GROMPP_ONLY=1; shift ;;
        --mdp-dir)     MDP_DIR="$2"; shift 2 ;;
        --md-mdp)      MD_MDP="$2"; shift 2 ;;
        --ntomp)       NTOMP="$2"; shift 2 ;;
        --gpu-id)      GPU_ID="$2"; shift 2 ;;
        --maxwarn)     MAXWARN="$2"; shift 2 ;;
        --quick-test)  QUICK=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) usage; die "unknown option: $1" ;;
    esac
done

[[ -n "$SYSTEM" ]] || { usage; die "--system is required"; }
[[ -d "$WORKDIR" ]] || die "workdir not found: $WORKDIR"
WORKDIR="$(cd "$WORKDIR" && pwd)"

# ------------------------- per-system defaults --------------------------------
case "$SYSTEM" in
    pp)
        : "${MDP_DIR:=$REPO_DIR/protein-protein/mdp}"
        MDP_EM=minim.mdp; MDP_NVT=nvt.mdp; MDP_NPT=npt.mdp
        MDP_MD100=md.mdp; MDP_MD500=md_500ns.mdp ;;
    apo|holo)
        : "${MDP_DIR:=$REPO_DIR/protein-ligand/$SYSTEM/mdp}"
        MDP_EM=EM.mdp; MDP_NVT=NVT.mdp; MDP_NPT=NPT.mdp
        MDP_MD100=MD.mdp; MDP_MD500=MD_500ns.mdp ;;
    *) die "--system must be pp, apo or holo (got: $SYSTEM)" ;;
esac

if [[ -n "$MD_MDP" ]]; then
    MDP_MD="$MD_MDP"
else
    case "$LENGTH" in
        100) MDP_MD="$MDP_DIR/$MDP_MD100" ;;
        500) MDP_MD="$MDP_DIR/$MDP_MD500" ;;
        *)   die "--length must be 100 or 500 (or use --md-mdp for a custom length)" ;;
    esac
fi
MDP_EM="$MDP_DIR/$MDP_EM"; MDP_NVT="$MDP_DIR/$MDP_NVT"; MDP_NPT="$MDP_DIR/$MDP_NPT"
for f in "$MDP_EM" "$MDP_NVT" "$MDP_NPT" "$MDP_MD"; do
    [[ -f "$f" ]] || die "mdp file not found: $f"
done

cd "$WORKDIR"

# Inputs: structure, topology, index
if [[ -z "$START" ]]; then
    if [[ -f solv_ions.gro ]]; then
        START=solv_ions.gro
    else
        mapfile -t cands < <(ls -1 *_ion.pdb *_ion.gro 2>/dev/null || true)
        if [[ ${#cands[@]} -eq 1 ]]; then START="${cands[0]}"
        else die "cannot guess the starting structure; pass --start FILE"; fi
    fi
fi
if [[ -z "$TOP" ]]; then
    if [[ -f topol.top ]]; then TOP=topol.top
    elif [[ -f LM.top ]]; then TOP=LM.top
    else die "no topology found (topol.top / LM.top); pass --top FILE"; fi
fi
if [[ -z "$INDEX" && -f index.ndx ]]; then INDEX=index.ndx; fi
[[ -f "$START" ]] || die "start structure not found: $START"
[[ -f "$TOP" ]]   || die "topology not found: $TOP"
if [[ "$SYSTEM" == "holo" && -z "$INDEX" ]]; then
    die "holo needs an index file with the Protein_UNL group (see protein-ligand/holo/README.md, step 11)"
fi
[[ -z "$INDEX" || -f "$INDEX" ]] || die "index file not found: $INDEX"

# Quick test: work on a copy so the real run starts clean
if [[ $QUICK -eq 1 ]]; then
    QDIR="$WORKDIR/quicktest"
    rm -rf "$QDIR"; mkdir -p "$QDIR"
    find "$WORKDIR" -maxdepth 1 -type f \( -name '*.gro' -o -name '*.pdb' -o -name '*.top' \
        -o -name '*.itp' -o -name '*.ndx' \) -exec cp {} "$QDIR/" \;
    for d in "$WORKDIR"/*.ff; do [[ -d "$d" ]] && cp -r "$d" "$QDIR/"; done
    cd "$QDIR"; WORKDIR="$QDIR"
    log "QUICK TEST in $QDIR (500 steps per stage)"
fi

load_gmx
: "${NTOMP:=$(ncpus)}"
GPU_MD_FLAGS="$(gpu_flags_md)"
COMMON_MDRUN=(-ntmpi 1 -ntomp "$NTOMP")
[[ -n "$GPU_ID" ]] && COMMON_MDRUN+=(-gpu_id "$GPU_ID")
[[ $QUICK -eq 1 ]] && COMMON_MDRUN+=(-nsteps 500)

mkdir -p logs mdp_used
cp "$MDP_EM" mdp_used/em.mdp;  cp "$MDP_NVT" mdp_used/nvt.mdp
cp "$MDP_NPT" mdp_used/npt.mdp; cp "$MDP_MD" mdp_used/md.mdp

IDX_ARGS=(); [[ -n "$INDEX" ]] && IDX_ARGS=(-n "$INDEX")

maxwarn_for() {   # stage -> -maxwarn value
    if [[ -n "$MAXWARN" ]]; then echo "$MAXWARN"; return; fi
    # pp: GROMOS force fields always trigger one grompp warning (twin-range note)
    # apo/holo: the Berendsen barostat in NPT.mdp triggers one warning
    if [[ "$SYSTEM" == "pp" ]]; then echo 1
    elif [[ "$1" == "npt" ]]; then echo 1
    else echo 0; fi
}

stage_done() { [[ -f "$1.gro" ]] && grep -q "Finished mdrun" "$1.log" 2>/dev/null; }

status() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" >> STATUS; log "$*"; }

log "system=$SYSTEM  length=${LENGTH} ns  workdir=$WORKDIR"
log "start=$START  top=$TOP  index=${INDEX:-<default groups>}  ntomp=$NTOMP  gpu=${GPU_MD_FLAGS:-none}"
log "production mdp: $MDP_MD"

# ------------------------------- stages ---------------------------------------
run_stage() {
    local name=$1 mdp=$2 last=$3; shift 3
    local grompp_args=("$@")
    local glog="logs/${name}.grompp.log" mlog="logs/${name}.mdrun.log"

    if stage_done "$name"; then status "$name: already finished, skipping"; return 0; fi

    local mdrun_extra=()
    if [[ "$name" != "em" ]]; then
        # shellcheck disable=SC2206
        mdrun_extra=($GPU_MD_FLAGS -cpt 15)
    fi

    if [[ "$name" != "em" && -f "$name.cpt" && -f "$name.tpr" ]]; then
        status "$name: checkpoint found, resuming"
        gmx mdrun -deffnm "$name" -cpi "$name.cpt" "${COMMON_MDRUN[@]}" "${mdrun_extra[@]}" \
            2>&1 | tee -a "$mlog" || die "$name: mdrun failed (log: $WORKDIR/$mlog)"
    else
        status "$name: grompp"
        gmx grompp -f "$mdp" "${grompp_args[@]}" -p "$TOP" "${IDX_ARGS[@]}" \
            -o "$name.tpr" -po "mdp_used/${name}.mdout.mdp" -maxwarn "$(maxwarn_for "$name")" \
            > "$glog" 2>&1 || { tail -n 40 "$glog"; die "$name: grompp failed (full log: $WORKDIR/$glog)"; }
        if [[ $GROMPP_ONLY -eq 1 && $last -eq 1 ]]; then
            status "$name: $name.tpr written (--grompp-only), stopping"
            return 0
        fi
        status "$name: mdrun"
        gmx mdrun -deffnm "$name" "${COMMON_MDRUN[@]}" "${mdrun_extra[@]}" 2>&1 | tee -a "$mlog" \
            || die "$name: mdrun failed (log: $WORKDIR/$mlog)"
    fi

    if [[ $QUICK -eq 1 ]]; then
        [[ -f "$name.gro" ]] || die "$name: no $name.gro written"
    else
        stage_done "$name" || die "$name: mdrun stopped before the end. Re-run the same command to resume."
    fi
    local perf
    perf=$(grep -E '^Performance:' "$name.log" 2>/dev/null | awk '{print $2}' || true)
    status "$name: finished${perf:+ (${perf} ns/day)}"
}

read -r -a STAGE_LIST <<< "$STAGES"
n=${#STAGE_LIST[@]}
for i in "${!STAGE_LIST[@]}"; do
    s=${STAGE_LIST[$i]}; last=0; [[ $i -eq $((n-1)) ]] && last=1
    case "$s" in
        em)  run_stage em  "$MDP_EM"  $last -c "$START" ;;
        nvt) [[ -f em.gro ]] || die "nvt needs em.gro"
             run_stage nvt "$MDP_NVT" $last -c em.gro -r em.gro ;;
        npt) [[ -f nvt.gro && -f nvt.cpt ]] || die "npt needs nvt.gro and nvt.cpt"
             run_stage npt "$MDP_NPT" $last -c nvt.gro -r nvt.gro -t nvt.cpt ;;
        md)  [[ -f npt.gro && -f npt.cpt ]] || die "md needs npt.gro and npt.cpt"
             run_stage md  "$MDP_MD"  $last -c npt.gro -t npt.cpt ;;
        *)   die "unknown stage: $s (use em, nvt, npt, md)" ;;
    esac
done

status "pipeline done: $STAGES"
if [[ $QUICK -eq 1 ]]; then
    log "Quick test passed. Start the real run with the same command without --quick-test."
fi
