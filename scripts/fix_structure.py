#!/usr/bin/env python3
"""
fix_structure.py — repair a protein PDB before `gmx pdb2gmx`

Fixes the problems that most often make pdb2gmx / grompp fail or make
energy minimisation blow up:

  * missing residues (gaps in the chain)        → modelled (short gaps) or
                                                  split with TER (long gaps)
  * missing heavy atoms (truncated side chains,
    missing OXT, ...)                           → rebuilt
  * non-standard residues (MSE, HYP, CSO, ...)  → replaced by the parent residue
  * heterogens (waters, ions, ligands, buffer)  → removed (optionally kept)
  * alternate locations                         → first one kept

Missing residues are only known when the file has SEQRES records (files
downloaded from the PDB have them). Docking outputs, AlphaFold models or
structures saved by Chimera/PyMOL often do not: pass the full sequence with
--fasta so the gaps can be found.

Requires PDBFixer + OpenMM:
    conda install -c conda-forge pdbfixer openmm
`--check-only` works without them (pure-Python gap and atom check).

Examples
    python scripts/fix_structure.py 3bik.pdb -o 3bik_fixed.pdb
    python scripts/fix_structure.py docked_protein.pdb -o pro_fixed.pdb --fasta target.fasta
    python scripts/fix_structure.py model.pdb --check-only

Exit codes: 0 = OK, 1 = error, 3 = chain breaks remain (split with TER).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import OrderedDict
from pathlib import Path

THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
STANDARD = set(THREE.values()) | {"HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "ASH", "GLH", "LYN"}
PEPTIDE_BOND_MAX = 2.0  # Å; a C(i)–N(i+1) distance above this is a chain break


# ----------------------------------------------------------------------------
# Lightweight PDB reading (no dependencies) — used for checks and TER splitting
# ----------------------------------------------------------------------------
def read_residues(pdb_text: str):
    """Return OrderedDict {(chain, resseq, icode): {"name": resname, "atoms": {name: (x,y,z)}, "het": bool}}"""
    res = OrderedDict()
    for line in pdb_text.splitlines():
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        alt = line[16]
        if alt not in (" ", "A", "1"):
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21]
        resseq = int(line[22:26])
        icode = line[26]
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        key = (chain, resseq, icode)
        r = res.setdefault(key, {"name": resname, "atoms": {}, "het": rec == "HETATM"})
        r["atoms"].setdefault(name, (x, y, z))
    return res


def dist(a, b):
    return math.sqrt(sum((p - q) ** 2 for p, q in zip(a, b)))


def find_breaks(residues):
    """Chain breaks between consecutive protein residues of the same chain."""
    breaks = []
    keys = [k for k, r in residues.items() if r["name"] in STANDARD or "CA" in r["atoms"]]
    for k1, k2 in zip(keys, keys[1:]):
        if k1[0] != k2[0]:
            continue
        r1, r2 = residues[k1], residues[k2]
        c, n = r1["atoms"].get("C"), r2["atoms"].get("N")
        if c is None or n is None:
            continue
        d = dist(c, n)
        if d > PEPTIDE_BOND_MAX:
            breaks.append({"chain": k1[0], "after": f"{r1['name']}{k1[1]}{k1[2].strip()}",
                           "before": f"{r2['name']}{k2[1]}{k2[2].strip()}", "C_N_distance_A": round(d, 2),
                           "numbering_gap": k2[1] - k1[1] - 1})
    return breaks


BACKBONE = ("N", "CA", "C", "O")


def find_incomplete_backbone(residues):
    out = []
    for (ch, num, ic), r in residues.items():
        if r["name"] in STANDARD:
            miss = [a for a in BACKBONE if a not in r["atoms"]]
            if miss:
                out.append(f"{ch}:{r['name']}{num}{ic.strip()} missing {','.join(miss)}")
    return out


def split_at_breaks(pdb_text: str) -> tuple[str, int]:
    """Insert TER records where a chain break remains so pdb2gmx does not bond across it."""
    residues = read_residues(pdb_text)
    breaks = find_breaks(residues)
    if not breaks:
        return pdb_text, 0
    starts = set()
    for b in breaks:
        # the residue that starts a new segment: parse "RESNUMicode" back from 'before'
        name = b["before"]
        num = "".join(ch for ch in name[3:] if ch.isdigit() or ch == "-")
        icode = name[3 + len(num):] or " "
        starts.add((b["chain"], int(num), icode))
    out, prev_key = [], None
    for line in pdb_text.splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            key = (line[21], int(line[22:26]), line[26])
            if key in starts and prev_key is not None and prev_key != key and out and not out[-1].startswith("TER"):
                out.append("TER")
            prev_key = key
        out.append(line)
    return "\n".join(out) + "\n", len(breaks)


# ----------------------------------------------------------------------------
# FASTA → PDBFixer sequences
# ----------------------------------------------------------------------------
def read_fasta(path: Path):
    records, cur = [], None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            cur = {"header": line[1:].strip(), "seq": ""}
            records.append(cur)
        elif cur is not None:
            cur["seq"] += "".join(ch for ch in line.upper() if ch.isalpha())
    return records


def chain_id_from_header(header: str, known: list[str]):
    """Accept '>A', '>chain A', '>3BIK_A', '>3BIK_1|Chains A, B|...' styles."""
    h = header.replace("|", " ").replace(",", " ").replace("_", " ").split()
    for tok in h[::-1]:
        if tok in known:
            return tok
    if "Chains" in header or "Chain" in header:
        for tok in h:
            if tok in known:
                return tok
    return None


# ----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path, help="input PDB")
    ap.add_argument("-o", "--output", type=Path, help="fixed PDB to write")
    ap.add_argument("--fasta", type=Path, help="full sequence (one record per chain) when the PDB has no SEQRES")
    ap.add_argument("--max-gap", type=int, default=15,
                    help="model internal gaps up to this many residues (default 15); longer gaps are split with TER")
    ap.add_argument("--add-terminal", action="store_true",
                    help="also model missing N-/C-terminal tails (off by default: tails are usually disordered)")
    ap.add_argument("--keep-het", default="", help="comma-separated residue names to keep, e.g. UNL,ZN")
    ap.add_argument("--keep-water", action="store_true", help="keep crystallographic waters")
    ap.add_argument("--add-h", action="store_true", help="add hydrogens at --ph (not needed with pdb2gmx -ignh)")
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--report", type=Path, help="write a JSON report")
    ap.add_argument("--check-only", action="store_true", help="only report problems (no PDBFixer needed)")
    ap.add_argument("--seed", type=int, default=0, help="random seed for rebuilt atoms (reproducible)")
    args = ap.parse_args(argv)

    if not args.input.exists():
        print(f"ERROR: {args.input} not found", file=sys.stderr)
        return 1
    raw = args.input.read_text(errors="ignore")
    before = read_residues(raw)
    report = {
        "input": str(args.input),
        "has_seqres": any(l.startswith("SEQRES") for l in raw.splitlines()),
        "chains": sorted({k[0] for k in before}),
        "breaks_before": find_breaks(before),
        "incomplete_backbone_before": find_incomplete_backbone(before),
        "heterogens_before": sorted({r["name"] for r in before.values() if r["name"] not in STANDARD and r["name"] != "HOH"}),
    }

    print(f"== {args.input.name}: {len(before)} residues, chains {','.join(report['chains'])}, "
          f"SEQRES {'yes' if report['has_seqres'] else 'no'}")
    for b in report["breaks_before"]:
        print(f"   gap  chain {b['chain']}: {b['after']} → {b['before']}  (C–N {b['C_N_distance_A']} Å, "
              f"{b['numbering_gap']} residue numbers skipped)")
    for m in report["incomplete_backbone_before"][:10]:
        print(f"   backbone  {m}")
    if report["heterogens_before"]:
        print(f"   heterogens: {', '.join(report['heterogens_before'])}")

    if args.check_only:
        if args.report:
            args.report.write_text(json.dumps(report, indent=2))
        problems = report["breaks_before"] or report["incomplete_backbone_before"]
        print("== check: " + ("problems found (run without --check-only to fix)" if problems else "no gaps found"))
        return 3 if problems else 0

    if args.output is None:
        print("ERROR: -o/--output is required unless --check-only", file=sys.stderr)
        return 1

    try:
        from pdbfixer import PDBFixer
        from pdbfixer.pdbfixer import Sequence
        from openmm.app import PDBFile
    except ImportError:
        print("ERROR: PDBFixer/OpenMM not installed.\n"
              "       conda install -c conda-forge pdbfixer openmm\n"
              "   or  pip install openmm pdbfixer", file=sys.stderr)
        return 1

    fixer = PDBFixer(filename=str(args.input))
    chains = list(fixer.topology.chains())
    chain_ids = [c.id for c in chains]

    # Sequence from FASTA (for files without SEQRES)
    if args.fasta:
        recs = read_fasta(args.fasta)
        seqs = []
        for i, rec in enumerate(recs):
            cid = chain_id_from_header(rec["header"], chain_ids)
            if cid is None:
                if i < len(chain_ids):
                    cid = chain_ids[i]
                else:
                    print(f"   (FASTA record '{rec['header']}' has no matching chain; ignored)")
                    continue
            bad = [ch for ch in rec["seq"] if ch not in THREE]
            if bad:
                print(f"ERROR: unknown residue letters in FASTA record '{rec['header']}': {sorted(set(bad))}",
                      file=sys.stderr)
                return 1
            seqs.append(Sequence(cid, [THREE[ch] for ch in rec["seq"]]))
            print(f"   FASTA: chain {cid} ← '{rec['header']}' ({len(rec['seq'])} aa)")
        fixer.sequences = seqs
    elif not report["has_seqres"] and report["breaks_before"]:
        print("   WARNING: gaps found but no SEQRES — residue names of the gaps are unknown.\n"
              "            Pass the full sequence with --fasta to model them; otherwise they are split with TER.")

    # 1) non-standard residues
    fixer.findNonstandardResidues()
    report["nonstandard_replaced"] = [f"{r.chain.id}:{r.name}{r.id}→{p}" for r, p in fixer.nonstandardResidues]
    replaced_names = {r.name for r, _ in fixer.nonstandardResidues}
    fixer.replaceNonstandardResidues()

    # 2) heterogens
    keep = {x.strip().upper() for x in args.keep_het.split(",") if x.strip()}
    kept_lines = []
    if keep:
        kept_lines = [l for l in raw.splitlines()
                      if l[:6] in ("ATOM  ", "HETATM") and l[17:20].strip().upper() in keep]
    fixer.removeHeterogens(keepWater=args.keep_water)
    report["heterogens_removed"] = [h for h in report["heterogens_before"]
                                    if h.upper() not in keep and h not in replaced_names]

    # 3) missing residues
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    modelled, skipped, drop = [], [], []
    for (ci, ri), names in sorted(fixer.missingResidues.items()):
        chain = chains[ci]
        nres = len(list(chain.residues()))
        terminal = ri == 0 or ri == nres
        where = "N-terminus" if ri == 0 else ("C-terminus" if ri == nres else f"after residue #{ri}")
        entry = {"chain": chain.id, "where": where, "n": len(names), "residues": "-".join(names)}
        if (terminal and not args.add_terminal) or (not terminal and len(names) > args.max_gap):
            entry["reason"] = "terminal tail (use --add-terminal)" if terminal else f"gap > --max-gap {args.max_gap}"
            skipped.append(entry)
            drop.append((ci, ri))
        else:
            modelled.append(entry)
    for key in drop:  # PDBFixer will not build these
        del fixer.missingResidues[key]
    report["missing_residues_modelled"] = modelled
    report["missing_residues_skipped"] = skipped

    # 4) missing heavy atoms (+ the residues kept above)
    fixer.findMissingAtoms()
    report["residues_with_missing_atoms"] = len(fixer.missingAtoms)
    report["missing_atoms_added"] = sum(len(v) for v in fixer.missingAtoms.values())
    report["missing_terminals_added"] = sum(len(v) for v in fixer.missingTerminals.values())
    fixer.addMissingAtoms(seed=args.seed)
    if args.add_h:
        fixer.addMissingHydrogens(args.ph)

    # 5) write, re-attach kept heterogens, split remaining breaks with TER
    tmp = args.output.with_suffix(".tmp.pdb")
    with open(tmp, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    text = tmp.read_text()
    tmp.unlink()
    lines = [l for l in text.splitlines() if not l.startswith(("END", "CONECT"))]
    if kept_lines:
        lines.append("TER")
        lines += kept_lines
    text = "\n".join(lines) + "\nEND\n"
    text, n_split = split_at_breaks(text)
    args.output.write_text(text)

    after = read_residues(text)
    report["breaks_after"] = find_breaks(after)
    report["incomplete_backbone_after"] = find_incomplete_backbone(after)
    report["output"] = str(args.output)

    print(f"== fixed → {args.output}")
    for x in report["nonstandard_replaced"]:
        print(f"   replaced  {x}")
    if report["heterogens_removed"]:
        print(f"   removed   {', '.join(report['heterogens_removed'])}")
    for m in modelled:
        print(f"   modelled  chain {m['chain']} {m['where']}: {m['n']} residues ({m['residues']})")
    for s in skipped:
        print(f"   SKIPPED   chain {s['chain']} {s['where']}: {s['n']} residues — {s['reason']}")
    print(f"   rebuilt   {report['missing_atoms_added']} heavy atoms in {report['residues_with_missing_atoms']} residues"
          f" (+{report['missing_terminals_added']} terminal atoms)")
    if n_split:
        print(f"   WARNING   {n_split} chain break(s) remain → split with TER (separate chain segments in pdb2gmx).")
        print("             For a long missing loop, consider building it with MODELLER / AlphaFold instead.")
    if report["incomplete_backbone_after"]:
        print(f"   WARNING   incomplete backbone remains: {report['incomplete_backbone_after'][:5]}")

    if args.report:
        args.report.write_text(json.dumps(report, indent=2))
    return 3 if n_split else 0


if __name__ == "__main__":
    sys.exit(main())
