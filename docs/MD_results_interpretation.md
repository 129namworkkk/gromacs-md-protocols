# Reading and interpreting MD results (protein–ligand)

**MD does not "prove that a drug works"**, but it provides quantitative evidence that the ligand binds stably, interacts plausibly, and deserves priority for experimental testing.

---

## 1. The big picture: which questions does MD answer?

After docking, MD usually addresses three main questions:

| # | Question | Metrics |
|---|----------|---------|
| 1 | Is the complex geometrically stable? | RMSD, ligand RMSD, RMSF, Rg, SASA, COM distance, clustering… |
| 2 | Are the ligand–pocket interactions stable? | H-bonds (and occupancy), contacts, mindist, water bridges, key residues… |
| 3 | Is binding energetically favourable? | MM/GBSA (ΔG_bind), per-residue decomposition (which residues contribute), comparison with other ligands / poses |

A complete MD study should provide all three kinds of evidence.

---

## 2. What each metric means and how to read it

### A) RMSD (protein & ligand)

**Protein RMSD (APO vs HOLO)**
- *Meaning:* is the overall protein structure stable?
- *Quick reading:*
  - RMSD levels off after the initial phase → the system has reached relative equilibrium.
  - A lower / more stable HOLO RMSD than APO sometimes suggests that the ligand "stabilises" the structure.
  - ⚠️ *Red flag:* a very large spike while SASA does not change → usually a PBC artefact.

**Ligand RMSD (fitted on the protein)**
- *Meaning:* does the ligand keep its pose in the pocket?
- *Quick reading:*
  - Low ligand RMSD without strong drift → the ligand does not "drift out of the pocket".
  - Ligand RMSD rising and then plateauing at a high value: the ligand changed pose but stays bound; use clustering to see whether it moved to a new stable pose.

### B) RMSF (Cα) (APO vs HOLO)
- *Meaning:* which residues are flexible; does the ligand make the pocket more rigid or more flexible?
- *Quick reading:*
  - Lower RMSF of pocket residues in HOLO → the ligand tends to "lock" the pocket.
  - A key loop with much higher RMSF in HOLO → the ligand may cause local instability.

### C) SASA (protein / ligand / complex)
- **Protein SASA (APO vs HOLO):** does the protein open or close?
- **Ligand SASA (in HOLO):** does the ligand stay buried?
- **Complex SASA:** complements the picture of how compact the whole complex is.

*Reasonable interpretation:*
- Low and stable ligand SASA → the ligand is little exposed to solvent → usually read as sitting deep / stably in the pocket (good).
- Protein SASA without abnormal jumps → no significant unfolding.

### D) Rg (protein)
- *Meaning:* overall compactness of the protein.
- *Quick reading:*
  - Stable Rg → no abnormal expansion or opening.
  - A slightly lower HOLO Rg than APO sometimes suggests a more compact complex.

### E) H-bonds (protein–ligand)
- *Meaning:* are hydrogen bonds maintained?
- *Quick reading:*
  - The number of H-bonds does not need to be large; what matters is stability and involvement of key residues.
  - Add per-residue H-bond occupancy (see section 3.4).

### F) DSSP (secondary structure)
- *Meaning:* does the protein keep its α/β structure (APO vs HOLO)?
- *Quick reading:*
  - If HOLO keeps the secondary structure more stable (fewer coil↔helix transitions) → supports "stabilisation".
  - Strong secondary-structure changes in the pocket or a functional region → interpret with care.

### G) MM/GBSA
- *Meaning:* estimated binding free energy ΔG_bind over time / as mean ± SD, with its decomposition into terms.
- *Quick reading:* a stable negative ΔTOTAL → favourable binding according to the model.
- *The terms tell you how the ligand binds:*

| Term | Typical sign | Meaning |
|------|--------------|---------|
| ΔVDWAALS | large negative | hydrophobic packing / shape complementarity is the main driver |
| ΔEGB | positive | polar solvation penalty (common) |
| ΔESURF | negative | favourable non-polar solvation (burial) |
| ΔEEL | small | direct electrostatics is not the main driver |

With a clearly negative ΔTOTAL dominated by vdW, the textbook interpretation is: binding is driven by hydrophobic packing and shape complementarity; polar solvation partly opposes it, but the total remains favourable.

---

## 3. Highly recommended additional analyses

Once the core analyses are done, these make the study stronger and closer to publication style:

1. **Contacts & minimum distance (interface stability)**: `ncontacts(t)` between protein and ligand (within a cut-off, e.g. 0.35–0.45 nm); `mindist(t)` protein–ligand.
   *Why it is strong:* it directly shows that "the ligand stays in contact", not only via RMSD.
2. **COM distance (does the ligand leave the pocket?)**: distance COM(ligand)–COM(pocket) or COM(protein).
   *Why it is strong:* ligand RMSD can increase because of rotation / pose changes; COM distance shows clearly whether it leaves the pocket.
3. **Clustering (dominant binding pose)**: cluster the stable part (e.g. 20–100 ns). Report the top 1–3 clusters + % occupancy + representative pose.
   *Why it is strong:* papers like a "dominant binding mode".
4. **Per-residue H-bond occupancy**: not just "how many H-bonds", but which H-bond lasts longest (%) and which residues are key.
   *Why it is strong:* tells the mechanistic story and guides mutations / assays.
5. **Per-residue energy decomposition (MM/GBSA decomp)**: table of the residues contributing most (top favourable residues).
   *Why it is strong:* links energy ↔ pocket residues ↔ rational drug design.
6. **PCA / essential dynamics (APO vs HOLO)**: compare the main motion subspaces; does HOLO reduce pocket fluctuations?
   *Why it is strong:* shows how the ligand affects global dynamics.
7. **Replicate MD (essential for strong conclusions)**: run 2–3 MD replicas (different seeds) and check that RMSD / Rg / SASA / MM/GBSA agree.
   *Why it is strong:* avoids conclusions based on "one lucky trajectory".

---

## 4. Minimum "evidence stack" to argue that a complex deserves in vitro / in vivo testing

| Layer | Evidence |
|-------|----------|
| 1. Structural stability | stable RMSD / Rg / SASA; no major DSSP changes |
| 2. Interaction stability | ligand RMSD + contacts + mindist + H-bond occupancy |
| 3. Energetics | stable negative MM/GBSA ΔG_bind + (ideally) per-residue decomposition |
| 4. Binding mode | clustering shows 1–2 dominant binding modes |
| 5. Reliability | replicate runs, or comparison with a control ligand / known inhibitor |

With all five layers, the case for a "promising complex" is well rounded.

---

## 5. Important notes for newcomers

- MD + MM/GBSA **does not prove biological efficacy**; it only provides mechanistic / dynamic evidence of stable binding.
- Phrase conclusions as *"supports / is consistent with / likely"* rather than *"confirms"*.
- When comparing ligands, keep **the same protocol** (simulation length, output interval, cut-offs, force field, PBC handling, equilibrated window).
