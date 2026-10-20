# Troubleshooting

| Triệu chứng | Nguyên nhân thường gặp | Cách xử lý |
|-------------|------------------------|-----------|
| `pdb2gmx`: *atom CG used in that entry is not found* | residue thiếu atom | **tự động:** `scripts/prepare_system.sh` · thủ công: [PDBFixer](#1-thiếu-atom--residue-pdb2gmx) |
| EM bùng nổ / LINCS warning ở một vị trí, bond rất dài | **thiếu residue** (gap) — pdb2gmx nối thẳng hai đầu gap | **tự động:** `scripts/prepare_system.sh` (model gap hoặc tách chain bằng TER) |
| `Residue 'XXX' not found in residue topology database` | residue không chuẩn (MSE, HOH, ligand, ion) còn trong file | `scripts/fix_structure.py` thay/loại; ligand tách riêng (`--ligand`) |
| `grompp`: không tìm thấy nhóm `Water_and_ions` / `Protein_UNL` | index chưa có nhóm tương ứng với `tc-grps` | [Tạo lại nhóm](#2-lỗi-nhóm-tc-grps) |
| `trjconv`: mismatch số atom | chọn Output *System* khi `.xtc` chỉ có Protein | chọn Output trùng với nội dung `.xtc` (`gmx check -f file.xtc`) |
| RMSD có spike rất lớn, SASA không đổi | PBC artifact | chạy `-pbc nojump` rồi `-pbc mol -center -ur compact` (xem [ANALYSIS.md](../protein-ligand/ANALYSIS.md)) |
| `gmx_MMPBSA`: *ChamberParm does not support all potential terms…* | topology GROMACS không chuyển được sang Amber | [Build topology Amber cho ligand](#3-gmx_mmpbsa-chamberparm-error) |
| `gmx_MMPBSA`: lỗi `AMBERHOME` | AmberTools chưa được activate | `conda activate gmxMMPBSA` rồi chạy lại |

---

## 1. Thiếu atom / residue (pdb2gmx)

### Cách tự động (khuyên dùng)

```bash
# chỉ kiểm tra (không cần cài gì thêm)
python scripts/fix_structure.py pro.pdb --check-only

# sửa + pdb2gmx + box + nước + ion + index trong một lệnh
bash scripts/prepare_system.sh --system apo --pdb pro.pdb --workdir runs/apo
#   file không có SEQRES (docking, AlphaFold, Chimera)? thêm --fasta target.fasta
```

`prepare_system.sh` luôn sửa cấu trúc trước `pdb2gmx`; nếu `pdb2gmx` vẫn báo thiếu atom, script sửa lại và chạy lại `pdb2gmx` (kể cả khi dùng `--no-fix`). Báo cáo chi tiết ở `fix_report.json` và `prepare.log`.

| Tình huống | Kết quả |
|------------|---------|
| Gap trong chuỗi ≤ `--max-gap` (mặc định 15) và biết trình tự (SEQRES hoặc `--fasta`) | model residue thiếu |
| Gap dài hơn, hoặc không có trình tự | tách chain bằng `TER` + cảnh báo (`--strict` để dừng thay vì tách) |
| Đuôi N/C-terminal thiếu | bỏ qua (thường là vùng mất trật tự); `--add-terminal` để model |
| Side chain cụt, thiếu OXT | dựng lại |

> Loop dài (> 15 residue) model bằng PDBFixer chỉ là hình học đơn giản. Với vùng quan trọng, dùng MODELLER hoặc ghép từ mô hình AlphaFold rồi mới chạy `prepare_system.sh`.

### Cách thủ công

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

# chạy lại
gmx pdb2gmx -f pro_fixed.pdb -o processed.gro -p topol.top -ter -ff oplsaa -water spc
```

## 2. Lỗi nhóm `tc-grps`

Kiểm tra `solv_ions.gro` có ion hay không. Nếu `grompp` vẫn báo lỗi nhóm, tạo lại bằng `make_ndx`:

```bash
gmx make_ndx -f em.gro -o index.ndx
```

```text
> 1 | 13                     # Protein | UNL   → Protein_UNL   (HOLO)
> r SOL | r NA | r CL        # → đổi tên thành Water_and_ions
> name XX Water_and_ions
> q
```

## 3. gmx_MMPBSA: ChamberParm error

Build topology Amber (GAFF2) cho ligand:

```bash
antechamber -i LIG.pdb -fi pdb -o LIG.mol2 -fo mol2 -c bcc -nc YOUR_CHARGE -at gaff2
parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod
```

Tạo `tleap.in`:

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

ls -l tleap.in LIG.mol2 LIG.frcmod     # kiểm tra file đã có
tleap -f tleap.in
```

Thành công khi có đủ:

```bash
ls -l LIG.lib LIG.prmtop LIG.inpcrd
```
