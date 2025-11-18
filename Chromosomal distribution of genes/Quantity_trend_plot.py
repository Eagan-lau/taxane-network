#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Adaptive version based on tidy_counts (Family, Hap, Chr, ChrIndex, Count):

Functions included:
1. Convert per-family per-haplotype chromosome counts into percentages 
   (i.e., column-wise normalization → sum = 100% within each Family × Hap).
2. Trend lines: two curves (hap1 and hap2) with an auto-adjusted Y-axis.
3. Heatmaps: one for hap1 and one for hap2, with automatically scaled color ranges.
4. Optional stacked-bar plots: for each chromosome, visualize family composition 
   as 100% column-normalized stacked bars.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ================= Input Configuration =================
# If using Excel input, use: INPUT = "family_chrom_counts.xlsx" with SHEET="tidy_counts"
INPUT = "family_chrom_counts.xlsx"
IS_EXCEL = INPUT.endswith(".xlsx")
SHEET = "tidy_counts"   # Used only when IS_EXCEL=True

OUTDIR = Path("plots_counts_pct_auto")
OUTDIR.mkdir(exist_ok=True)

# ================= Read Input Table =================
if IS_EXCEL:
    df = pd.read_excel(INPUT, sheet_name=SHEET)
else:
    df = pd.read_csv(INPUT)

# Standardize column types
df["Family"] = df["Family"].astype(str)
df["Hap"] = df["Hap"].astype(str)
df["Chr"] = df["Chr"].astype(str)
df["ChrIndex"] = pd.to_numeric(df["ChrIndex"], errors="coerce").astype("Int64")
df["Count"] = pd.to_numeric(df["Count"], errors="coerce").fillna(0).astype(int)

# Helper: convert a Series to within-column percent
def series_to_pct(s: pd.Series) -> pd.Series:
    tot = s.sum()
    return (s / tot * 100.0) if tot > 0 else s * 0.0

# ===========================================================
# 1) Family-level percentage trend curves (hap1 vs hap2)
# ===========================================================
def plot_family_trend_pct_auto(family: str, save_pdf: bool = True):
    sub = df[df["Family"] == family].copy()
    if sub.empty:
        print("[Skip]", family, "has no data")
        return

    # Pivot table: rows = ChrIndex, columns = Hap, values = Count
    pv_count = (
        sub.pivot_table(
            index="ChrIndex",
            columns="Hap",
            values="Count",
            aggfunc="sum"
        )
        .fillna(0)
        .sort_index()
    )

    # Column-wise percentage normalization per hap
    pv_pct = pv_count.copy()
    for col in pv_pct.columns:
        pv_pct[col] = series_to_pct(pv_pct[col])

    fig, ax = plt.subplots(figsize=(8, 3.2))

    # Two trend lines (hap1 and hap2)
    for col in ["hap1", "hap2"]:
        if col in pv_pct.columns:
            ax.plot(
                pv_pct.index,
                pv_pct[col],
                marker="o",
                linewidth=2,
                label=col
            )

    ax.set_title(f"{family} — Chromosome share within family (%, by hap)", fontsize=12)
    ax.set_xlabel("Chromosome index")
    ax.set_ylabel("Percent (%)")
    ax.grid(alpha=0.3, linestyle="--", linewidth=0.8)
    ax.set_xticks(pv_pct.index.tolist())

    # Auto-adjust Y-axis margins
    yvals = pv_pct.values.flatten()
    yvals = yvals[np.isfinite(yvals)]
    if yvals.size:
        ymin, ymax = float(np.min(yvals)), float(np.max(yvals))
        pad = (ymax - ymin) * 0.08 if ymax > ymin else 1.0
        ax.set_ylim(ymin - pad, ymax + pad)

    ax.legend()
    fig.tight_layout()

    fig.savefig(OUTDIR / f"trend_pct_{family}.png", dpi=300)
    if save_pdf:
        fig.savefig(OUTDIR / f"trend_pct_{family}.pdf")
    plt.close(fig)


for fam in sorted(df["Family"].unique()):
    plot_family_trend_pct_auto(fam)

print("All percentage trend plots have been saved to:", OUTDIR)

# ===========================================================
# 2) Percentage heatmap (one for hap1, one for hap2)
# ===========================================================
def plot_heatmap_pct_auto(hap: str, cmap: str = "viridis"):
    sub = df[df["Hap"] == hap].copy()
    if sub.empty:
        print("[Skip]", hap, "has no data")
        return

    # Pivot: rows = Family, columns = ChrIndex, values = Count
    mat = (
        sub.pivot_table(
            index="Family",
            columns="ChrIndex",
            values="Count",
            aggfunc="sum"
        )
        .fillna(0)
        .sort_index(axis=1)
    )

    # Row-wise percentage conversion (share of each chromosome for a given family)
    mat_pct = mat.apply(series_to_pct, axis=1)

    families = mat_pct.index.tolist()
    chridx = mat_pct.columns.tolist()
    arr = mat_pct.values

    fig_w = max(8, len(chridx) * 0.6)
    fig_h = max(4, len(families) * 0.3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Auto-scaled color domain (no vmin/vmax specified)
    im = ax.imshow(arr, aspect="auto", cmap=cmap)

    ax.set_title(f"Chromosome share within family (%, {hap})", fontsize=12)
    ax.set_xticks(range(len(chridx)))
    ax.set_xticklabels(chridx)
    ax.set_yticks(range(len(families)))
    ax.set_yticklabels(families)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Percent (%)")

    fig.tight_layout()
    fig.savefig(OUTDIR / f"heatmap_pct_{hap}.png", dpi=300)
    fig.savefig(OUTDIR / f"heatmap_pct_{hap}.pdf")
    plt.close(fig)

plot_heatmap_pct_auto("hap1")
plot_heatmap_pct_auto("hap2")

print("All percentage heatmaps have been saved to:", OUTDIR)

# ===========================================================
# 3) (Optional) Stacked 100% bars: family composition per chromosome
# ===========================================================
def plot_stacked_by_chr_pct_auto(hap: str):
    sub = df[df["Hap"] == hap].copy()
    if sub.empty:
        print("[Skip]", hap, "has no data")
        return

    # Pivot: rows = ChrIndex, columns = Family, values = Count
    wide = (
        sub.pivot_table(
            index="ChrIndex",
            columns="Family",
            values="Count",
            aggfunc="sum"
        )
        .fillna(0)
        .sort_index()
    )

    # Column-normalization (within each chromosome)
    colsum = wide.sum(axis=1).replace(0, np.nan)
    wide_pct = (wide.T / colsum).T.fillna(0.0) * 100.0

    families = list(wide_pct.columns)
    x = np.arange(len(wide_pct.index))
    bottoms = np.zeros(len(x))

    fig, ax = plt.subplots(figsize=(max(8, len(x) * 0.5), 4))

    # Stacked bars
    for fam in families:
        ax.bar(
            x,
            wide_pct[fam].values,
            bottom=bottoms,
            label=fam
        )
        bottoms += wide_pct[fam].values

    ax.set_title(f"Family composition per chromosome (%, {hap})", fontsize=12)
    ax.set_xlabel("Chromosome index")
    ax.set_ylabel("Percent (%)")
    ax.set_xticks(x)
    ax.set_xticklabels(wide_pct.index.tolist())

    # Auto Y limit (expected ≈100%)
    ymax = bottoms.max() if bottoms.size else 100.0
    ax.set_ylim(0, ymax * 1.05)

    ax.legend(
        ncol=min(4, len(families)),
        fontsize=8,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )
    ax.grid(axis="y", alpha=0.3, linestyle="--", linewidth=0.8)

    fig.tight_layout()
    fig.savefig(OUTDIR / f"stacked_pct_{hap}.png", dpi=300)
    fig.savefig(OUTDIR / f"stacked_pct_{hap}.pdf")
    plt.close(fig)


plot_stacked_by_chr_pct_auto("hap1")
plot_stacked_by_chr_pct_auto("hap2")

print("All stacked percentage bar plots have been saved to:", OUTDIR)
