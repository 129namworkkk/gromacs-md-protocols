# Installation (Ubuntu / WSL)

This guide is shared by the **protein–protein** and **protein–ligand** workflows.

- [1. System packages](#1-system-packages)
- [2. Build GROMACS 2023.3 with CUDA (GPU) support](#2-build-gromacs-20233-with-cuda-gpu-support)
- [3. Load GROMACS in every new Ubuntu session](#3-load-gromacs-in-every-new-ubuntu-session)
- [4. Python environment](#4-python-environment)
- [5. gmx_MMPBSA (MM/GBSA)](#5-gmx_mmpbsa-mmgbsa)
- [6. Supporting software](#6-supporting-software)

---

## 1. System packages

```bash
sudo apt update
sudo apt upgrade
sudo apt install -y gcc cmake build-essential git wget
sudo apt install -y libfftw3-dev
sudo apt-get install -y grace          # xmgrace, quick look at .xvg files
sudo apt-get install -y pymol
```

Quickest option (CPU only, apt package), fine for testing:

```bash
sudo apt install gromacs -y
```

---

## 2. Build GROMACS 2023.3 with CUDA (GPU) support

### 2.1 Check the GPU and driver

```bash
nvidia-smi
```

If the GPU name and driver version are shown, the system sees the GPU.

### 2.2 Install the CUDA toolkit and compilers

```bash
sudo apt update
sudo apt install -y cmake build-essential git wget nvidia-cuda-toolkit
nvcc --version

# GCC 12 avoids CUDA build errors
sudo apt install -y gcc-12 g++-12
```

### 2.3 Download the GROMACS 2023.3 source

```bash
cd ~
wget https://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
tar xzf gromacs-2023.3.tar.gz
```

### 2.4 Configure with CMake (CUDA + FFTW, no MPI)

```bash
rm -rf ~/gromacs-2023.3-build
mkdir -p ~/gromacs-2023.3-build
cd ~/gromacs-2023.3-build

cmake ../gromacs-2023.3 \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_GPU=CUDA \
  -DGMX_CUDA_TARGET_SM=89 \
  -DGMX_SIMD=AVX2_256 \
  -DGMX_MPI=OFF \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-12 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 \
  -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-12 \
  -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2023.3-cuda
```

| Option | Meaning |
|--------|---------|
| `GMX_GPU=CUDA` | enable CUDA GPU support |
| `GMX_CUDA_TARGET_SM=89` | GPU architecture of RTX 40-series cards (change for your GPU: 80 = A100, 86 = RTX 30 / A40, 89 = RTX 40 / L40S, 90 = H100) |
| `GMX_SIMD=AVX2_256` | optimised for modern CPUs |
| `GMX_BUILD_OWN_FFTW=ON` | avoids FFTW problems |
| `GMX_MPI=OFF` | stable, simple (thread-MPI) |
| GCC 12 | good compatibility with CUDA |

### 2.5 Build and install

```bash
make -j$(nproc)      # or make -j16
make install
```

### 2.6 Check the build

```bash
source $HOME/gromacs-2023.3-cuda/bin/GMXRC
gmx --version
```

Expected output:

```text
GROMACS version:    2023.3
Precision:          mixed
Memory model:       64 bit
MPI library:        thread_mpi
OpenMP support:     enabled (GMX_OPENMP_MAX_THREADS = 128)
GPU support:        CUDA
NB cluster size:    8
SIMD instructions:  AVX2_256
```

### 2.7 (Optional) Rebuild with MPI (multi-node)

```bash
rm -rf ~/gromacs-2023.3-build
mkdir -p ~/gromacs-2023.3-build
cd ~/gromacs-2023.3-build

cmake ../gromacs-2023.3 \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DGMX_GPU=CUDA \
  -DGMX_CUDA_TARGET_SM=89 \
  -DGMX_SIMD=AVX2_256 \
  -DGMX_MPI=ON \
  -DCMAKE_INSTALL_PREFIX=$HOME/gromacs-2023.3-cuda

make -j$(nproc)
make install
source $HOME/gromacs-2023.3-cuda/bin/GMXRC
gmx --version        # MPI: ON
```

---

## 3. Load GROMACS in every new Ubuntu session

```bash
# use the GPU build
source ~/gromacs-2023.3-cuda/bin/GMXRC
```

---

## 4. Python environment

`gmx-ff` environment (Anaconda prompt):

```bash
conda create -n gmx-ff
conda activate gmx-ff
conda install -c conda-forge python=3.7.12 numpy=1.21.5 numpy-base=1.21.5 blas=1.0 \
  mkl=2021.4.0 mkl-service=2.4.0 mkl_fft=1.3.1 mkl_random=1.2.2 intel-openmp=2021.4.0 \
  pip=24.0 setuptools=69.0.3 wheel=0.42.0 six=1.16.0 openssl=3.6.0 \
  ca-certificates=2025.11.4 sqlite=3.51.0 libsqlite=3.51.0 networkx=2.3
```

For plotting with [`scripts/plot_compare_metrics.py`](../scripts/plot_compare_metrics.py) and structure repair with [`scripts/fix_structure.py`](../scripts/fix_structure.py):

```bash
pip install -r requirements.txt     # numpy, matplotlib
pip install openmm pdbfixer         # or: conda install -c conda-forge pdbfixer openmm
```

---

## 5. gmx_MMPBSA (MM/GBSA)

```bash
conda update conda
conda create -n gmxMMPBSA python=3.11.8 -y
conda activate gmxMMPBSA
conda install -c conda-forge "mpi4py=4.0.1" "ambertools<=23.3" -y
python -m pip install gmx_MMPBSA

# (recommended) MPI for faster runs
sudo apt update
sudo apt install openmpi-bin libopenmpi-dev openssh-client

# check
gmx_MMPBSA -h
```

> An error saying `AMBERHOME` is not set usually means AmberTools is not activated correctly; gmx_MMPBSA requires `AMBERHOME` to exist.

---

## 6. Supporting software

| Tool | Used for | Link |
|------|----------|------|
| UCSF Chimera | adding hydrogens to ligands, splitting / converting PDB files | https://www.cgl.ucsf.edu/chimera/download.html |
| LigParGen | OPLS-AA ligand topologies | https://zarbi.chem.yale.edu/ligpargen/ |
| PDBFixer (OpenMM) | automatic repair of missing residues / atoms (`scripts/fix_structure.py`, `prepare_system.sh`) | `conda install -c conda-forge pdbfixer openmm -y` or `pip install openmm pdbfixer` |
| PyMOL / VMD | viewing structures and trajectories | — |
| xmgrace | quick look at `.xvg` files | `sudo apt-get install grace` |
