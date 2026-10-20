# Cài đặt môi trường (Ubuntu / WSL)

Hướng dẫn này dùng chung cho cả quy trình **Protein–Protein** và **Protein–Ligand**.

- [1. Gói hệ thống](#1-gói-hệ-thống)
- [2. Build GROMACS 2023.3 có hỗ trợ CUDA (GPU)](#2-build-gromacs-20233-có-hỗ-trợ-cuda-gpu)
- [3. Nạp GROMACS mỗi lần mở Ubuntu](#3-nạp-gromacs-mỗi-lần-mở-ubuntu)
- [4. Môi trường Python](#4-môi-trường-python)
- [5. gmx_MMPBSA (MM/GBSA)](#5-gmx_mmpbsa-mmgbsa)
- [6. Phần mềm phụ trợ](#6-phần-mềm-phụ-trợ)

---

## 1. Gói hệ thống

```bash
sudo apt update
sudo apt upgrade
sudo apt install -y gcc cmake build-essential git wget
sudo apt install -y libfftw3-dev
sudo apt-get install -y grace          # xmgrace để xem nhanh .xvg
sudo apt-get install -y pymol
```

Cách nhanh nhất (CPU, bản apt), dùng để chạy thử:

```bash
sudo apt install gromacs -y
```

---

## 2. Build GROMACS 2023.3 có hỗ trợ CUDA (GPU)

### 2.1 Kiểm tra GPU + driver

```bash
nvidia-smi
```

Nếu lệnh hiển thị tên GPU, driver version… thì hệ thống đã nhận GPU.

### 2.2 Cài CUDA toolkit + compiler

```bash
sudo apt update
sudo apt install -y cmake build-essential git wget nvidia-cuda-toolkit
nvcc --version

# GCC 12 để tránh lỗi build CUDA
sudo apt install -y gcc-12 g++-12
```

### 2.3 Tải source GROMACS 2023.3

```bash
cd ~
wget https://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
tar xzf gromacs-2023.3.tar.gz
```

### 2.4 Cấu hình CMake (CUDA + FFTW, không MPI)

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

| Tuỳ chọn | Ý nghĩa |
|----------|---------|
| `GMX_GPU=CUDA` | bật hỗ trợ GPU CUDA |
| `GMX_CUDA_TARGET_SM=89` | kiến trúc GPU RTX 40-series (đổi theo GPU của bạn) |
| `GMX_SIMD=AVX2_256` | tối ưu cho CPU hiện đại |
| `GMX_BUILD_OWN_FFTW=ON` | tránh lỗi FFTW |
| `GMX_MPI=OFF` | ổn định, dễ dùng (thread-MPI) |
| GCC-12 | tương thích tốt với CUDA |

### 2.5 Build & install

```bash
make -j$(nproc)      # hoặc make -j16
make install
```

### 2.6 Kiểm tra

```bash
source $HOME/gromacs-2023.3-cuda/bin/GMXRC
gmx --version
```

Kết quả đúng cần thấy:

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

### 2.7 (Tuỳ chọn) Build lại với MPI (multi-node)

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

## 3. Nạp GROMACS mỗi lần mở Ubuntu

```bash
# để bật GPU chạy cho nhanh
source ~/gromacs-2023.3-cuda/bin/GMXRC
```

---

## 4. Môi trường Python

Môi trường `gmx-ff` (Anaconda prompt):

```bash
conda create -n gmx-ff
conda activate gmx-ff
conda install -c conda-forge python=3.7.12 numpy=1.21.5 numpy-base=1.21.5 blas=1.0 \
  mkl=2021.4.0 mkl-service=2.4.0 mkl_fft=1.3.1 mkl_random=1.2.2 intel-openmp=2021.4.0 \
  pip=24.0 setuptools=69.0.3 wheel=0.42.0 six=1.16.0 openssl=3.6.0 \
  ca-certificates=2025.11.4 sqlite=3.51.0 libsqlite=3.51.0 networkx=2.3
```

Để vẽ đồ thị bằng [`scripts/plot_compare_metrics.py`](../scripts/plot_compare_metrics.py):

```bash
pip install -r requirements.txt     # numpy, matplotlib
```

---

## 5. gmx_MMPBSA (MM/GBSA)

```bash
conda update conda
conda create -n gmxMMPBSA python=3.11.8 -y
conda activate gmxMMPBSA
conda install -c conda-forge "mpi4py=4.0.1" "ambertools<=23.3" -y
python -m pip install gmx_MMPBSA

# (khuyến nghị) chạy MPI cho nhanh
sudo apt update
sudo apt install openmpi-bin libopenmpi-dev openssh-client

# kiểm tra
gmx_MMPBSA -h
```

> Nếu báo lỗi `AMBERHOME` chưa set thì thường do AmberTools chưa được activate đúng; gmx_MMPBSA yêu cầu `AMBERHOME` phải tồn tại.

---

## 6. Phần mềm phụ trợ

| Công cụ | Dùng để | Link |
|---------|---------|------|
| UCSF Chimera | thêm hydrogen cho ligand, tách/chuyển định dạng PDB | https://www.cgl.ucsf.edu/chimera/download.html |
| LigParGen | sinh topology OPLS-AA cho ligand | https://zarbi.chem.yale.edu/ligpargen/ |
| PDBFixer (OpenMM) | tự sửa missing residues / atoms (`scripts/fix_structure.py`, `prepare_system.sh`) | `conda install -c conda-forge pdbfixer openmm -y` hoặc `pip install openmm pdbfixer` |
| PyMOL / VMD | xem cấu trúc, trajectory | — |
| xmgrace | xem nhanh file `.xvg` | `sudo apt-get install grace` |
