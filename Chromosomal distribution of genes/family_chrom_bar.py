#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# ================= Configuration =================
XLSX = "family_chrom_counts.xlsx"
SHEET_H1 = "hap1_wide"   # Expected columns: TcChr01a, TcChr02a, ...
SHEET_H2 = "hap2_wide"   # Expected columns: TcChr01b, TcChr02b, ...
OUT_PNG = "family_chrom_bar_pairs.png"
DPI = 300

# If a specific list of five families should be plotted, define them here;
# If empty, the script selects the first five families in the table automatically.
FAMILIES_TO_PLOT = ['2-ODD', 'ACH', 'ACT', 'P450', 'UGT']

# Define chromosome ordering (1..12), pairing hap1(a) with hap2(b).
CHROM_IDX = list(range(1, 13))
CHROM_ORDER_PAIRS = [
    (f"TcChr{str(i).zfill(2)}a", f"TcChr{str(i).zfill(2)}b")
    for i in CHROM_IDX
]
# Flatten the order into a,b,a,b,...
CHROM_ORDER = [x for pair in CHROM_ORDER_PAIRS for x in pair]

# ================= Read and merge data =================
hap1 = pd.read_excel(XLSX, sheet_name=SHEET_H1).set_index('Family')
hap2 = pd.read_excel(XLSX, sheet_name=SHEET_H2).set_index('Family')

# Ensure all expected columns exist; missing ones are filled with zeros.
for col in CHROM_ORDER:
    if col.endswith('a') and col not in hap1.columns:
        hap1[col] = 0
    if col.endswith('b') and col not in hap2.columns:
        hap2[col] = 0

# Reindex individual haplotype tables and ensure integer counts.
hap1 = hap1.reindex(columns=[c for c in CHROM_ORDER if c.endswith('a')]).fillna(0).astype(int)
hap2 = hap2.reindex(columns=[c for c in CHROM_ORDER if c.endswith('b')]).fillna(0).astype(int)

# Combine hap1 and hap2 into one wide table (column order: a,b,a,b,...)
combined = pd.concat([hap1, hap2], axis=1)
combined = combined.reindex(columns=CHROM_ORDER).fillna(0).astype(int)

# Select the five families to be plotted
if FAMILIES_TO_PLOT:
    fam_list = [f for f in FAMILIES_TO_PLOT if f in combined.index][:5]
else:
    fam_list = combined.index.tolist()[:5]

if len(fam_list) < 1:
    raise ValueError("No valid families found for plotting. Check the Excel file or FAMILIES_TO_PLOT.")

# ================= Plotting =================
# Font and vector export settings
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

n_rows = len(fam_list)
fig_h = max(2.2 * n_rows, 4.5)
fig_w = 16

fig, axes = plt.subplots(n_rows, 1, figsize=(fig_w, fig_h), sharex=True)
if n_rows == 1:
    axes = [axes]

# ---- Spacing configuration ----
bar_w = 0.35       # Width of a single bar (a or b)
inner_gap = 0.00   # Gap between sister chromosomes (a vs b); 0 means touching.
gap = 0.60         # Gap between chromosome pairs.

# Compute x positions so that a,b are adjacent and pairs separated by a gap.
x_positions = []
pair_centers = []
base = 0.0

for a_col, b_col in CHROM_ORDER_PAIRS:
    x_a = base
    x_b = base + bar_w + inner_gap
    x_positions.extend([x_a, x_b])
    pair_centers.append((x_a + x_b) / 2.0)

    # Move base for the next chromosome pair
    base += (bar_w + inner_gap) + bar_w + gap

# Colors for hap1 (a) and hap2 (b)
color_a = "#4C78A8"
color_b = "#F58518"

# ---------------- Per-family bar plots ----------------
for ax, fam in zip(axes, fam_list):
    vals = combined.loc[fam, CHROM_ORDER].values.astype(float)

    xs = np.array(x_positions)
    ys_a = vals[0::2]  # Values for ...a
    ys_b = vals[1::2]  # Values for ...b

    # Draw bars for hap1 and hap2
    ax.bar(
        xs[0::2], ys_a,
        width=bar_w,
        label="hap1 (…a)",
        color=color_a,
        edgecolor="black",
        linewidth=0.6
    )
    ax.bar(
        xs[1::2], ys_b,
        width=bar_w,
        label="hap2 (…b)",
        color=color_b,
        edgecolor="black",
        linewidth=0.6
    )

    # Family label on the left of each row
    ax.set_ylabel(
        fam, rotation=0, ha="right", va="center",
        labelpad=24, fontsize=11
    )

    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# ---------------- Bottom x-axis ----------------
major_labels = [f"TcChr{str(i).zfill(2)}" for i in CHROM_IDX]
axes[-1].set_xticks(pair_centers)
axes[-1].set_xticklabels(major_labels, rotation=0, fontsize=10)

# Legend shown on the first subplot
axes[0].legend(ncol=2, frameon=False, loc="upper right")

fig.tight_layout()
plt.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight")
plt.close()

print("Saved figure:", OUT_PNG)
