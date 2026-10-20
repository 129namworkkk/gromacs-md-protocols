#!/usr/bin/env bash
# =============================================================================
# benchmark.sh — measure ns/day on this pod and estimate run time and cost
#
# Needs a production .tpr. Create it without starting the long run:
#   bash runpod/run_md.sh --system holo --length 500 --workdir DIR --stages "em nvt npt md" --grompp-only
# then:
#   bash runpod/benchmark.sh --workdir DIR --price 0.69     # price = USD per hour of this pod
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

WORKDIR="$PWD"; TPR="md.tpr"; STEPS=20000; PRICE=""; NTOMP=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --workdir) WORKDIR="$2"; shift 2 ;;
        --tpr)     TPR="$2"; shift 2 ;;
        --steps)   STEPS="$2"; shift 2 ;;
        --price)   PRICE="$2"; shift 2 ;;
        --ntomp)   NTOMP="$2"; shift 2 ;;
        -h|--help) sed -n '2,10p' "$0"; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done

cd "$WORKDIR"
[[ -f "$TPR" ]] || die "$TPR not found in $WORKDIR (see the header of this script)"
load_gmx
: "${NTOMP:=$(ncpus)}"
# shellcheck disable=SC2046
gmx mdrun -s "$TPR" -deffnm bench -nsteps "$STEPS" -resethway -noconfout \
    -ntmpi 1 -ntomp "$NTOMP" $(gpu_flags_md) > logs_bench.txt 2>&1 \
    || { tail -n 30 logs_bench.txt; die "benchmark mdrun failed"; }

NSDAY=$(grep -E '^Performance:' bench.log | awk '{print $2}')
rm -f bench.xtc bench.trr bench.edr bench.cpt bench.gro logs_bench.txt
[[ -n "$NSDAY" ]] || die "could not read performance from bench.log"

awk -v p="$NSDAY" -v price="$PRICE" 'BEGIN {
    printf "\nPerformance on this pod: %.1f ns/day\n\n", p
    printf "  %-8s %10s %10s", "length", "hours", "days"
    if (price != "") printf " %12s", "cost (USD)"
    printf "\n"
    n = split("100 200 500 1000", L, " ")
    for (i = 1; i <= n; i++) {
        h = L[i] / p * 24
        printf "  %-8s %10.1f %10.2f", L[i] " ns", h, h / 24
        if (price != "") printf " %12.2f", h * price
        printf "\n"
    }
    printf "\n(estimate for production only; equilibration adds minutes)\n"
}'
