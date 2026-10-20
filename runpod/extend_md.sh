#!/usr/bin/env bash
# =============================================================================
# extend_md.sh — continue a finished (or running) production run to a longer total
#
# Example: a 100 ns run is done and the reviewer asks for 500 ns:
#   bash runpod/extend_md.sh --workdir /workspace/runs/my_complex --to 500
#
# The trajectory, energies and log keep being appended to md.xtc / md.edr / md.log.
# Output frequencies stay those of the original .mdp (they cannot be changed
# by convert-tpr). Re-running the same command after an interruption resumes.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

WORKDIR="$PWD"; NAME="md"; TO_NS=""; NTOMP=""; GPU_ID=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --workdir) WORKDIR="$2"; shift 2 ;;
        --name)    NAME="$2"; shift 2 ;;
        --to)      TO_NS="$2"; shift 2 ;;
        --ntomp)   NTOMP="$2"; shift 2 ;;
        --gpu-id)  GPU_ID="$2"; shift 2 ;;
        -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done
[[ -n "$TO_NS" ]] || die "--to <total length in ns> is required"
[[ "$TO_NS" =~ ^[0-9]+$ ]] || die "--to must be an integer number of ns"

cd "$WORKDIR"
[[ -f "$NAME.tpr" && -f "$NAME.cpt" ]] || die "$NAME.tpr and $NAME.cpt are required in $WORKDIR"
load_gmx
: "${NTOMP:=$(ncpus)}"
mkdir -p logs

MARK=".extended_to_${TO_NS}ns"
if [[ ! -f "$MARK" ]]; then
    cp "$NAME.tpr" "${NAME}_before_extend.tpr"
    gmx convert-tpr -s "${NAME}_before_extend.tpr" -until "$((TO_NS * 1000))" -o "$NAME.tpr" \
        > "logs/${NAME}.extend.log" 2>&1 || { cat "logs/${NAME}.extend.log"; die "convert-tpr failed"; }
    touch "$MARK"
    echo "$(date '+%Y-%m-%d %H:%M:%S') $NAME: tpr extended to $TO_NS ns" >> STATUS
    log "$NAME.tpr now ends at $TO_NS ns (original kept as ${NAME}_before_extend.tpr)"
fi

EXTRA=(); [[ -n "$GPU_ID" ]] && EXTRA+=(-gpu_id "$GPU_ID")
# shellcheck disable=SC2046
gmx mdrun -deffnm "$NAME" -cpi "$NAME.cpt" -ntmpi 1 -ntomp "$NTOMP" -cpt 15 \
    $(gpu_flags_md) "${EXTRA[@]}" 2>&1 | tee -a "logs/${NAME}.mdrun.log" \
    || die "mdrun failed — re-run the same command to resume"

grep -q "Finished mdrun" "$NAME.log" || die "mdrun stopped early — re-run the same command to resume"
echo "$(date '+%Y-%m-%d %H:%M:%S') $NAME: extended run finished ($TO_NS ns)" >> STATUS
log "Extension to $TO_NS ns finished."
