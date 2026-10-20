# Đọc và diễn giải kết quả MD (Protein–Ligand)

> Tài liệu hướng dẫn đọc kết quả MD (tiếng Việt).

**MD không "chứng minh thuốc hiệu quả"**, nhưng tạo bộ bằng chứng định lượng rằng ligand bám ổn + tương tác hợp lý + có cơ sở ưu tiên thử nghiệm.

---

## 1. Bức tranh tổng thể: MD đang trả lời câu hỏi gì?

Sau docking, MD thường trả lời 3 câu hỏi lớn:

| # | Câu hỏi | Chỉ số |
|---|---------|--------|
| 1 | Complex có bền hình học không? | RMSD, ligand RMSD, RMSF, Rg, SASA, COM distance, clustering… |
| 2 | Tương tác ligand–pocket có ổn định không? | H-bond (và occupancy), contacts, mindist, water bridges, key residues… |
| 3 | Binding có "thuận lợi về năng lượng" không? | MM/GBSA (ΔG_bind), decomposition theo residue (ai đóng góp), so sánh với ligand khác/pose khác |

Muốn "toàn vẹn" cho quy trình MD → bạn nên có đủ 3 nhóm bằng chứng.

---

## 2. Giải thích và cách đọc từng chỉ số MD

### A) RMSD (protein & ligand)

**Protein RMSD (APO vs HOLO)**
- *Ý nghĩa:* cấu trúc tổng thể của protein có ổn định không.
- *Đọc nhanh:*
  - RMSD ổn định sau giai đoạn đầu → hệ đạt cân bằng tương đối.
  - HOLO RMSD thấp/ổn định hơn APO đôi khi gợi ý ligand "stabilize" cấu trúc.
  - ⚠️ *Red flags:* spike cực lớn mà SASA không đổi → thường là PBC artifact.

**Ligand RMSD (fit theo protein)**
- *Ý nghĩa:* ligand có giữ tư thế trong pocket không.
- *Đọc nhanh:*
  - Ligand RMSD thấp và không drift mạnh → ligand không "trôi khỏi pocket".
  - Nếu ligand RMSD tăng dần rồi plateau cao: ligand đổi pose nhưng vẫn bám; cần clustering để xem đổi sang pose ổn định mới.

### B) RMSF (Cα) (APO vs HOLO)
- *Ý nghĩa:* residue nào linh động; ligand làm pocket "cứng" hay "mềm" hơn.
- *Đọc nhanh:*
  - Pocket residues giảm RMSF ở HOLO → ligand có xu hướng "khóa" vùng pocket.
  - Nếu một loop quan trọng tăng RMSF mạnh ở HOLO → có thể ligand gây bất ổn cục bộ.

### C) SASA (Protein / Ligand / Complex)
- **Protein SASA (APO vs HOLO):** protein có mở/đóng không.
- **Ligand SASA (trong HOLO):** ligand có "buried" ổn định không.
- **Complex SASA:** bổ trợ cho "độ compact" của toàn complex.

*Cách diễn giải hợp lý:*
- Ligand SASA thấp và ổn định → ligand ít lộ dung môi → hay được hiểu là nằm sâu/ổn định trong pocket (tốt).
- Protein SASA không nhảy bất thường → không có unfolding đáng kể.

### D) Rg (Protein)
- *Ý nghĩa:* độ "gọn" toàn protein.
- *Đọc nhanh:*
  - Rg ổn định → không có giãn/mở bất thường.
  - HOLO Rg hơi thấp hơn APO đôi khi gợi ý complex "compact" hơn.

### E) H-bond (Protein–Ligand)
- *Ý nghĩa:* liên kết H-bond có duy trì không.
- *Đọc nhanh:*
  - H-bond count không cần quá nhiều; quan trọng là ổn định và đúng residue chủ chốt.
  - Nên bổ sung: H-bond occupancy theo từng residue (xem mục 3.4).

### F) DSSP (Secondary structure)
- *Ý nghĩa:* protein có giữ alpha/beta ổn định không (apo vs holo).
- *Đọc nhanh:*
  - Nếu holo làm SS ổn định hơn (ít chuyển coil↔helix) → hỗ trợ "stabilization".
  - Nếu pocket hoặc vùng chức năng bị biến đổi SS mạnh → cần cẩn trọng.

### G) MM/GBSA
- *Ý nghĩa:* ước lượng năng lượng liên kết ΔG_bind theo thời gian/mean±SD, và phân rã thành phần.
- *Đọc nhanh:* ΔTOTAL âm ổn định → binding thuận lợi theo mô hình.
- *Thành phần nói "ligand bám kiểu gì":*

| Thành phần | Dấu thường gặp | Ý nghĩa |
|-----------|----------------|---------|
| ΔVDWAALS | âm lớn | đóng gói kỵ nước / shape complementarity là driver chính |
| ΔEGB | dương | polar solvation penalty (thường gặp) |
| ΔESURF | âm | nonpolar solvation thuận lợi (burial) |
| ΔEEL | nhỏ | electrostatic trực tiếp không phải driver chính |

Với kết quả ΔTOTAL âm rõ + driver vdW, cách diễn giải "đúng bài" là: binding được dẫn dắt bởi packing kỵ nước và phù hợp hình dạng, dung môi phân cực triệt một phần nhưng tổng vẫn thuận lợi.

---

## 3. Những phân tích "rất nên thêm" để củng cố luận điểm complex tốt

Khi đã có phần "core", để bài chắc hơn và gần phong cách paper, thường bổ sung các mảnh sau:

1. **Contacts & min distance (interface stability)**: `ncontacts(t)` giữa Protein–Ligand (trong cutoff, ví dụ 0.35–0.45 nm); `mindist(t)` protein–ligand.
   *Vì sao mạnh:* nó trực tiếp chứng minh "ligand vẫn giữ tiếp xúc" chứ không chỉ RMSD.
2. **COM distance (ligand có trôi khỏi pocket không)**: khoảng cách COM(ligand)–COM(pocket) hoặc COM(protein).
   *Vì sao mạnh:* ligand RMSD đôi khi tăng vì xoay/đổi pose; COM distance mới nói rõ "có rời pocket không".
3. **Clustering (tư thế binding chủ đạo)**: chạy clustering cho đoạn ổn định (vd 20–100 ns). Report: top 1–3 cluster + % occupancy + representative pose.
   *Vì sao mạnh:* paper rất thích "dominant binding mode".
4. **H-bond occupancy theo residue**: không chỉ "bao nhiêu H-bond", mà là H-bond nào giữ lâu nhất (%) và residue nào là key.
   *Vì sao mạnh:* giúp kể câu chuyện cơ chế tương tác, chọn residue để đột biến/assay.
5. **Per-residue energy decomposition (MM/GBSA decomp)**: bảng residue đóng góp năng lượng (top favorable residues).
   *Vì sao mạnh:* kết nối năng lượng ↔ pocket residues ↔ rational drug design.
6. **PCA / Essential dynamics (apo vs holo)**: so sánh không gian chuyển động chính; holo có "giảm dao động" ở pocket không.
   *Vì sao mạnh:* cho thấy ligand ảnh hưởng động lực học toàn cục.
7. **Replicate MD (cực quan trọng khi muốn kết luận mạnh)**: chạy 2–3 lần MD (seed khác), so RMSD/Rg/SASA/MMGBSA có đồng thuận không.
   *Vì sao mạnh:* tránh kết luận dựa trên "1 trajectory may be lucky".

---

## 4. "Bộ bằng chứng" tối thiểu để nói: đáng thử in vitro / in vivo

| Lớp | Bằng chứng |
|-----|-----------|
| 1. Ổn định cấu trúc | RMSD/Rg/SASA ổn định; DSSP không biến dạng lớn |
| 2. Ổn định tương tác | ligand RMSD + contacts + mindist + H-bond occupancy |
| 3. Năng lượng | MM/GBSA ΔG_bind âm ổn định + (tốt nhất) residue decomposition |
| 4. Binding mode | clustering cho thấy 1–2 binding mode chiếm ưu thế |
| 5. Độ tin cậy | chạy replicate hoặc so với ligand đối chứng / known inhibitor |

Đủ 5 lớp này thì phần "complex tiềm năng" rất tròn bài.

---

## 5. Lưu ý quan trọng cho người mới

- MD + MM/GBSA **không chứng minh hiệu quả sinh học**, chỉ cung cấp bằng chứng cơ học/động lực học cho binding ổn định.
- Kết luận nên viết kiểu: *"support / consistent with / likely"* hơn là *"confirm"*.
- Khi muốn so sánh ligand, hãy giữ **cùng protocol** (thời gian, interval, cutoff, forcefield, xử lý PBC, đoạn ổn định).
