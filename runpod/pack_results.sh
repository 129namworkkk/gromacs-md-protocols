#!/usr/bin/env bash
# =============================================================================
# pack_results.sh — shrink and bundle results for download
#
# Writes, next to the run:
#   md_noWater.tpr / md_noWater.xtc   protein (+ ligand, ions) without water —
#                                     typically ~10x smaller; use them as
#                                     holo.tpr/holo.xtc (or apo.*) in ANALYSIS.md
#   md_final_noWater.pdb              last frame, for a quick look in PyMOL/Chimera
#   potential_em.xvg, temperature_nvt.xvg, pressure_npt.xvg, density_npt.xvg
# and bundles them with logs, energies and topology into results_<name>.tar.gz
#
# Usage:
#   bash runpod/pack_results.sh --workdir /workspace/runs/my_complex [--with-water]
#   --with-water also adds the full md.tpr + md.xtc (large) to the archive
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

WORKDIR="$PWD"; NAME="md"; WITH_WATER=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --workdir)    WORKDIR="$2"; shift 2 ;;
        --name)       NAME="$2"; shift 2 ;;
        --with-water) WITH_WATER=1; shift ;;
        -h|--help)    sed -n '2,17p' "$0"; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

cd "$WORKDIR"
[[ -f "$NAME.tpr" && -f "$NAME.xtc" ]] || die "$NAME.tpr and $NAME.xtc are required in $WORKDIR"
load_gmx
mkdir -p logs

# default index groups (System, Protein, non-Water, ...) straight from the tpr
echo q | gmx make_ndx -f "$NAME.tpr" -o pack.ndx > logs/pack.log 2>&1

log "Writing water-free tpr / trajectory / final frame"
echo non-Water | gmx convert-tpr -s "$NAME.tpr" -n pack.ndx -o "${NAME}_noWater.tpr" >> logs/pack.log 2>&1
echo non-Water | gmx trjconv -s "$NAME.tpr" -f "$NAME.xtc" -n pack.ndx -o "${NAME}_noWater.xtc" >> logs/pack.log 2>&1
if [[ -f "$NAME.gro" ]]; then
    echo non-Water | gmx trjconv -s "$NAME.tpr" -f "$NAME.gro" -n pack.ndx -pbc mol \
        -o "${NAME}_final_noWater.pdb" >> logs/pack.log 2>&1
fi

log "Extracting equilibration curves"
energy_term() {   # edr term out.xvg
    [[ -f "$1" ]] || return 0
    printf '%s\n0\n' "$2" | gmx energy -f "$1" -o "$3" >> logs/pack.log 2>&1 \
        || log "  (could not extract $2 from $1)"
}
energy_term em.edr  Potential   potential_em.xvg
energy_term nvt.edr Temperature temperature_nvt.xvg
energy_term npt.edr Pressure    pressure_npt.xvg
energy_term npt.edr Density     density_npt.xvg

ARCHIVE="results_$(basename "$WORKDIR")_${NAME}.tar.gz"
files=()
for f in "${NAME}_noWater.tpr" "${NAME}_noWater.xtc" "${NAME}_final_noWater.pdb" \
         "$NAME.edr" "$NAME.log" em.log nvt.log npt.log em.edr nvt.edr npt.edr \
         potential_em.xvg temperature_nvt.xvg pressure_npt.xvg density_npt.xvg \
         ./*.top ./*.itp ./*.ndx STATUS; do
    [[ -e "$f" ]] && files+=("$f")
done
[[ -d logs ]] && files+=(logs)
[[ -d mdp_used ]] && files+=(mdp_used)
if [[ $WITH_WATER -eq 1 ]]; then files+=("$NAME.tpr" "$NAME.xtc"); fi

tar czf "$ARCHIVE" "${files[@]}"
log "Archive: $WORKDIR/$ARCHIVE ($(du -h "$ARCHIVE" | cut -f1))"
log "Full trajectory stays on the volume: $WORKDIR/$NAME.xtc ($(du -h "$NAME.xtc" | cut -f1))"
