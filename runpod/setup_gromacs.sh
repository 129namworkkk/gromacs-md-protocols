#!/usr/bin/env bash
# =============================================================================
# setup_gromacs.sh — build GROMACS (CUDA) once on the RunPod network volume
#
# Installs into /workspace so it survives pod restarts: later pods attached to
# the same volume only need `source /workspace/gromacs-2023.3-cuda/bin/GMXRC`.
#
# Requirements: a pod started from a CUDA *devel* image (it must contain nvcc),
# e.g. the official "RunPod PyTorch" template. Build time ≈ 15–30 min.
#
# Environment variables (optional):
#   GMX_VERSION   default 2023.3 (same as the guides in this repository)
#   PREFIX        default /workspace/gromacs-$GMX_VERSION-cuda
#   CUDA_SM       GPU architectures to compile for, e.g. "89" or "80;86;89;90".
#                 Default: the GPU of this pod + 80;86;89 so the build also
#                 runs on A100 / A40 / RTX 4090 / L40S pods later.
# =============================================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "$SCRIPT_DIR/common.sh"

GMX_VERSION="${GMX_VERSION:-2023.3}"
PREFIX="${PREFIX:-/workspace/gromacs-${GMX_VERSION}-cuda}"
SRC_ROOT="${SRC_ROOT:-/workspace/src}"
SUDO=""; [[ $(id -u) -ne 0 ]] && SUDO="sudo"

if [[ -x "$PREFIX/bin/gmx" ]]; then
    log "GROMACS already installed in $PREFIX"
    "$PREFIX/bin/gmx" --version | grep -E 'GROMACS version|GPU support|SIMD instructions'
    log "Load it with: source $PREFIX/bin/GMXRC"
    exit 0
fi

# ---- GPU + CUDA compiler -----------------------------------------------------
nvidia-smi || die "no NVIDIA GPU visible — start a GPU pod"
if ! command -v nvcc >/dev/null 2>&1 && [[ -x /usr/local/cuda/bin/nvcc ]]; then
    export PATH="/usr/local/cuda/bin:$PATH"
fi
command -v nvcc >/dev/null 2>&1 || die "nvcc not found — use a CUDA *devel* image (e.g. RunPod PyTorch template)"
log "CUDA compiler: $(nvcc --version | grep -o 'release [0-9.]*')"

NVCC_VER=$(nvcc --version | grep -oE 'release [0-9]+\.[0-9]+' | awk '{print $2}')
ver_ge() { [[ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n1)" == "$2" ]]; }
if [[ -z "${CUDA_SM:-}" ]]; then
    this_sm=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -n1 | tr -d '. ' || true)
    CUDA_SM="80;86"
    ver_ge "$NVCC_VER" 11.8 && CUDA_SM="$CUDA_SM;89"     # RTX 40xx / L40S / L4 need CUDA >= 11.8
    if [[ -n "$this_sm" && ";$CUDA_SM;" != *";$this_sm;"* ]]; then CUDA_SM="$CUDA_SM;$this_sm"; fi
fi
log "Compiling for CUDA architectures: $CUDA_SM"

# ---- build dependencies --------------------------------------------------------
export DEBIAN_FRONTEND=noninteractive
$SUDO apt-get update -y
$SUDO apt-get install -y --no-install-recommends cmake build-essential wget ca-certificates libhwloc-dev

cmake_ver=$(cmake --version | head -n1 | grep -oE '[0-9]+\.[0-9]+')
if [[ "$(printf '%s\n3.18\n' "$cmake_ver" | sort -V | head -n1)" != "3.18" ]]; then
    log "cmake $cmake_ver is too old, installing a newer one with pip"
    python3 -m pip install --quiet "cmake>=3.18"
    hash -r
fi

# ---- download + build ------------------------------------------------------------
mkdir -p "$SRC_ROOT"; cd "$SRC_ROOT"
if [[ ! -d "gromacs-$GMX_VERSION" ]]; then
    wget -q "https://ftp.gromacs.org/gromacs/gromacs-$GMX_VERSION.tar.gz"
    tar xzf "gromacs-$GMX_VERSION.tar.gz"
fi
rm -rf "build-$GMX_VERSION"; mkdir "build-$GMX_VERSION"; cd "build-$GMX_VERSION"

JOBS=$(ncpus)
cmake "../gromacs-$GMX_VERSION" \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_GPU=CUDA \
    -DGMX_CUDA_TARGET_SM="$CUDA_SM" \
    -DGMX_MPI=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$PREFIX"
make -j"$JOBS"
make install

# shellcheck disable=SC1091
source "$PREFIX/bin/GMXRC"
gmx --version | grep -E 'GROMACS version|GPU support|SIMD instructions'
log "Done. In every new terminal / pod run:  source $PREFIX/bin/GMXRC"
