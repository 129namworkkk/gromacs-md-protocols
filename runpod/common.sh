#!/usr/bin/env bash
# Shared helpers for the RunPod pipeline scripts. Sourced, not executed.

# Timestamped log line
log() { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die() { log "ERROR: $*" >&2; exit 1; }

# Number of CPU cores actually available to this container.
# Pods often report the host's core count in `nproc` while the cgroup quota
# is much smaller; oversubscribing OpenMP threads slows GROMACS down.
ncpus() {
    local n q p c
    n=$(nproc)
    if [[ -r /sys/fs/cgroup/cpu.max ]]; then                 # cgroup v2
        read -r q p < /sys/fs/cgroup/cpu.max
        if [[ "$q" != "max" && "$p" -gt 0 ]]; then
            c=$(( (q + p - 1) / p ))
            (( c > 0 && c < n )) && n=$c
        fi
    elif [[ -r /sys/fs/cgroup/cpu/cpu.cfs_quota_us ]]; then  # cgroup v1
        q=$(cat /sys/fs/cgroup/cpu/cpu.cfs_quota_us)
        p=$(cat /sys/fs/cgroup/cpu/cpu.cfs_period_us)
        if (( q > 0 && p > 0 )); then
            c=$(( (q + p - 1) / p ))
            (( c > 0 && c < n )) && n=$c
        fi
    fi
    echo "$n"
}

# Make `gmx` available: use it from PATH, otherwise source a GMXRC
# from the persistent volume (installed by setup_gromacs.sh).
load_gmx() {
    if command -v gmx >/dev/null 2>&1; then return 0; fi
    local rc
    for rc in "${GMXRC:-}" /workspace/gromacs-*-cuda/bin/GMXRC "$HOME"/gromacs-*-cuda/bin/GMXRC; do
        if [[ -n "$rc" && -f "$rc" ]]; then
            # shellcheck disable=SC1090
            source "$rc"
            command -v gmx >/dev/null 2>&1 && return 0
        fi
    done
    die "gmx not found. Run runpod/setup_gromacs.sh first, or set GMXRC=/path/to/GMXRC."
}

# True when an NVIDIA GPU is visible AND this gmx build has GPU support.
has_gpu() {
    command -v nvidia-smi >/dev/null 2>&1 || return 1
    nvidia-smi -L >/dev/null 2>&1 || return 1
    gmx --version 2>/dev/null | grep -E '^GPU support:' | grep -qvi 'disabled'
}

# mdrun flags that offload work to the GPU for dynamics stages.
gpu_flags_md() {
    if has_gpu; then echo "-nb gpu -pme gpu -bonded gpu"; fi
}
