#!/usr/bin/env python3
# plot_compare_metrics.py
#
# Install:
#   pip install numpy matplotlib
#
# Run (PowerShell):
#   python plot_compare_metrics.py --base "C:\path\to\xvg" --outdir "plots" --dpi 600

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Readers: .xvg and simple .dat
# -----------------------------
_XLAB_RE = re.compile(r'@\s*xaxis\s+label\s+"([^"]+)"')
_YLAB_RE = re.compile(r'@\s*yaxis\s+label\s+"([^"]+)"')
_LEG_RE = re.compile(r'@\s*s(\d+)\s+legend\s+"([^"]+)"')


@dataclass
class XYData:
    x: np.ndarray
    ys: List[np.ndarray]  # supports multi-column y (x y1 y2 ...)
    xlab: str
    ylab: str
    legends: List[str]


def read_xvg_or_dat(path: Path) -> XYData:
    xlab: Optional[str] = None
    ylab: Optional[str] = None
    leg_map: Dict[int, str] = {}
    rows: List[List[float]] = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            # XVG headers start with @ or # ; some .dat may start with ; or #
            if line.startswith(("#", "@", ";")):
                if line.startswith("@"):
                    m = _XLAB_RE.search(line)
                    if m:
                        xlab = m.group(1)
                    m = _YLAB_RE.search(line)
                    if m:
                        ylab = m.group(1)
                    m = _LEG_RE.search(line)
                    if m:
                        leg_map[int(m.group(1))] = m.group(2)
                continue

            parts = line.split()
            try:
                nums = [float(p) for p in parts]
            except ValueError:
                continue

            if len(nums) >= 2:
                rows.append(nums)

    if not rows:
        raise ValueError(f"No numeric data found in: {path.name}")

    arr = np.asarray(rows, dtype=float)
    x = arr[:, 0]
    ys = [arr[:, i] for i in range(1, arr.shape[1])]
    legends = [leg_map.get(i, f"series_{i}") for i in range(len(ys))]

    return XYData(
        x=x,
        ys=ys,
        xlab=xlab or "X",
        ylab=ylab or "Y",
        legends=legends,
    )


def convert_time_to_ns_if_needed(x: np.ndarray, xlab: str) -> Tuple[np.ndarray, str]:
    lab = (xlab or "").lower()

    if "ps" in lab:
        return x / 1000.0, xlab.replace("ps", "ns").replace("PS", "ns")
    if "fs" in lab:
        return x / 1e6, xlab.replace("fs", "ns").replace("FS", "ns")

    # heuristic: if it's likely ps (e.g., 100000) and label is Time/generic
    if np.nanmax(x) > 5000 and ("time" in lab or xlab in ("X", "x", "")):
        return x / 1000.0, "Time (ns)"

    return x, xlab or "Time (ns)"


# -----------------------------
# Plotting
# -----------------------------
@dataclass
class SeriesSpec:
    path: Path
    label: str
    y_index: int = 0
    time_to_ns: bool = False


def plot_metric(
    title: str,
    series: List[SeriesSpec],
    out_base: Path,
    dpi: int = 350,
    legend_right: bool = True,
):
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 16,
            "axes.labelsize": 13,
            "legend.fontsize": 11,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
        }
    )

    fig, ax = plt.subplots(figsize=(18, 6))
    chosen_xlab = None
    chosen_ylab = None

    for s in series:
        if not s.path.exists():
            print(f"[WARN] Missing file, skipped: {s.path.name}")
            continue

        d = read_xvg_or_dat(s.path)

        x = d.x
        xlab = d.xlab
        if s.time_to_ns or "time" in (xlab or "").lower():
            x, xlab = convert_time_to_ns_if_needed(x, xlab)

        y_idx = min(max(s.y_index, 0), len(d.ys) - 1)
        y = d.ys[y_idx]

        chosen_xlab = chosen_xlab or xlab
        chosen_ylab = chosen_ylab or d.ylab

        ax.plot(x, y, label=s.label, linewidth=1.3)

    ax.set_title(title)
    ax.set_xlabel(chosen_xlab or "X")
    ax.set_ylabel(chosen_ylab or "Y")
    ax.grid(True, alpha=0.25)

    if legend_right:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            frameon=True,
        )
        fig.tight_layout(rect=[0, 0, 0.84, 1])
    else:
        ax.legend()
        fig.tight_layout()

    out_base.parent.mkdir(parents=True, exist_ok=True)

    png_path = out_base.with_suffix(".png")
    tiff_path = out_base.with_suffix(".tiff")

    fig.savefig(png_path, dpi=dpi)
    fig.savefig(tiff_path, dpi=dpi, format="tiff")

    plt.close(fig)
    print(f"[OK] {png_path}")
    print(f"[OK] {tiff_path}")


# -----------------------------
# Jobs (filenames)
# -----------------------------
def build_jobs(base: Path) -> List[Tuple[str, List[SeriesSpec], Path]]:
    jobs: List[Tuple[str, List[SeriesSpec], Path]] = []

    # RMSD: APO vs HOLO (+ ligand)
    jobs.append(
        (
            "RMSD: APO vs HOLO (+ Ligand)",
            [
                SeriesSpec(base / "rmsd_apo.xvg", "APO (Protein RMSD)", time_to_ns=True),
                SeriesSpec(base / "rmsd_holo.xvg", "HOLO (Protein RMSD)", time_to_ns=True),
                SeriesSpec(base / "rmsd_lig.xvg", "HOLO (Ligand RMSD)", time_to_ns=True),
            ],
            base / "RMSD_APO_vs_HOLO_plus_Ligand",
        )
    )

    # RMSF: APO vs HOLO
    jobs.append(
        (
            "RMSF: APO vs HOLO",
            [
                SeriesSpec(base / "rmsf_apo_at.xvg", "APO (RMSF)"),
                SeriesSpec(base / "rmsf_holo_at.xvg", "HOLO (RMSF)"),
            ],
            base / "RMSF_APO_vs_HOLO",
        )
    )

    # Rg: APO vs HOLO (if rg_apo exists; otherwise HOLO-only)
    rg_apo = base / "rg_apo.xvg"
    rg_series: List[SeriesSpec] = []
    if rg_apo.exists():
        rg_series.append(SeriesSpec(rg_apo, "APO (Rg)", time_to_ns=True))
        out_name = base / "Rg_APO_vs_HOLO"
        title = "Radius of gyration (Rg): APO vs HOLO"
    else:
        out_name = base / "Rg_HOLO_only"
        title = "Radius of gyration (Rg): HOLO"

    rg_series.append(SeriesSpec(base / "rg_holo.xvg", "HOLO (Rg)", time_to_ns=True))
    jobs.append((title, rg_series, out_name))

    # SASA: APO vs HOLO (Protein/Ligand/Complex)
    jobs.append(
        (
            "SASA: APO vs HOLO (Protein/Ligand/Complex)",
            [
                SeriesSpec(base / "sasa_protein_apo.xvg", "APO (Protein SASA)", time_to_ns=True),
                SeriesSpec(base / "sasa_protein_holo.xvg", "HOLO (Protein SASA)", time_to_ns=True),
                SeriesSpec(base / "sasa_lig_holo.xvg", "HOLO (Ligand SASA)", time_to_ns=True),
                SeriesSpec(base / "sasa_complex_holo.xvg", "HOLO (Complex SASA)", time_to_ns=True),
            ],
            base / "SASA_APO_vs_HOLO",
        )
    )

    # BSA: HOLO
    jobs.append(
        (
            "BSA / dSASA (HOLO)",
            [SeriesSpec(base / "bsa_holo_timeseries.xvg", "HOLO (BSA)", time_to_ns=True)],
            base / "BSA_HOLO",
        )
    )

    # H-bonds: Protein-Ligand (HOLO)
    jobs.append(
        (
            "H-bonds (Protein-Ligand, HOLO)",
            [SeriesSpec(base / "hb_prot_lig.xvg", "Protein–Ligand H-bonds", time_to_ns=True)],
            base / "HBonds_HOLO_prot_lig",
        )
    )

    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", type=str, default=".", help="Folder containing .xvg/.dat files")
    ap.add_argument("--outdir", type=str, default="plots", help="Output folder")
    ap.add_argument("--dpi", type=int, default=350, help="DPI for PNG/TIFF")
    args = ap.parse_args()

    base = Path(args.base).resolve()
    outdir = Path(args.outdir).resolve()

    jobs = build_jobs(base)
    for title, series, outbase in jobs:
        outbase = outdir / outbase.name
        plot_metric(title, series, outbase, dpi=args.dpi, legend_right=True)


if __name__ == "__main__":
    main()