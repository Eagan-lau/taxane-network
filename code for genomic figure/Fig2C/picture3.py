#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# ===== Fonts  =====
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42

# ------------------------
# Config
# ------------------------
INPUT_XLSX = "new2.xlsx"
SHEET_NAME = 0
OUT_PNG    = "gene_clusters_bar_lines_pies.png"
OUT_PDF    = "gene_clusters_bar_lines_pies.pdf"
FIGSIZE    = (16, 8)
BAR_ALPHA  = 0.8
BAR_EDGE   = "black"

# ======= Pie size controls =======

PIE_FRAC    = 1.40
PIE_MIN_PX  = 38
PIE_MAX_PX  = 150

PIE_GAP_FACTOR = 2.10

# ------------------------
# Helpers
# ------------------------
def normalize_percent_column(series):

    s = series.copy()
    if s.dtype == object:
        s = s.astype(str).str.strip()
        s = s.replace({"nan": np.nan})
        s = s.str.rstrip("%")
        s = pd.to_numeric(s, errors="coerce")
    else:
        s = pd.to_numeric(s, errors="coerce")

    if np.nanmax(s.values) <= 1.5:  
        s = s * 100.0
    return s

def coerce_numeric(s):
    return pd.to_numeric(s, errors="coerce")

def data_dx_to_inches(ax, fig, dx, frac=0.4, min_px=12, max_px=40):
    
    p0 = ax.transData.transform((0, 0))
    p1 = ax.transData.transform((dx, 0))
    dx_pix = abs(p1[0] - p0[0])
    pie_pix = np.clip(dx_pix * frac, min_px, max_px)
    return pie_pix / fig.dpi

# ------------------------
# Load
# ------------------------
df = pd.read_excel(INPUT_XLSX, sheet_name=SHEET_NAME)
df.columns = [str(c).strip() for c in df.columns]

if df.shape[1] < 18:
    raise ValueError("should contain at least 18 columns")


col_gene_label = df.iloc[:, 0].astype(str)        # 1st column as string labels
col_hap1_pct   = normalize_percent_column(df.iloc[:, 1])  # 2nd column
col_hap2_pct   = normalize_percent_column(df.iloc[:, 2])  # 3rd column
col_frequency  = coerce_numeric(df.iloc[:, 13])           # 14th column
pies_vals      = df.iloc[:, 14:18].apply(pd.to_numeric, errors="coerce").fillna(0.0)  # 15–18


mask_valid = col_frequency.notna()
col_gene_label = col_gene_label[mask_valid]
col_hap1_pct   = col_hap1_pct[mask_valid]
col_hap2_pct   = col_hap2_pct[mask_valid]
col_frequency  = col_frequency[mask_valid]
pies_vals      = pies_vals.loc[mask_valid]


x_labels = col_gene_label.astype(str).values
x_pos    = np.arange(len(x_labels), dtype=float)

bar_spacing = 1.0
bar_width   = bar_spacing * 0.8


fig, ax = plt.subplots(figsize=FIGSIZE)


bars = ax.bar(x_pos, col_frequency.values, width=bar_width, alpha=BAR_ALPHA,
              edgecolor=BAR_EDGE, label="Frequency")

ax.set_xlabel("Gene Number")
ax.set_ylabel("Frequency (count)")


positive_freq = col_frequency[col_frequency > 0]
ymin = (positive_freq.min() * 0.8) if len(positive_freq) else 1
ax.set_yscale('log')
ax.set_ylim(bottom=max(ymin, 1e-2))


ax2 = ax.twinx()
line1, = ax2.plot(x_pos, col_hap1_pct.values, marker="o", linewidth=1.8, label="Hap1 (%)")
line2, = ax2.plot(x_pos, col_hap2_pct.values, marker="s", linewidth=1.8, label="Hap2 (%)")


ax2.set_ylim(30, 70)
ax2.set_yticks([30, 40, 50, 60, 70])
ax2.yaxis.set_major_formatter(PercentFormatter(xmax=100))
ax2.set_ylabel("Percentage")


ax.set_xticks(x_pos)
ax.set_xticklabels(x_labels, rotation=0)


ax.grid(axis="y", linestyle="--", alpha=0.35, which="both")
handles1, labels1 = ax.get_legend_handles_labels()
handles2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(handles1 + [line1, line2], labels1 + ["Hap1 (%)", "Hap2 (%)"],
           loc="upper right", frameon=False)


pie_diam_in = data_dx_to_inches(ax, fig, dx=bar_spacing,
                                frac=PIE_FRAC, min_px=PIE_MIN_PX, max_px=PIE_MAX_PX)
pie_labels = ["Hap1-1", "Hap2-2", "Hap1-2", "others"]

for i, (x, freq) in enumerate(zip(x_pos, col_frequency.values)):
    vals = pies_vals.iloc[i].values.astype(float)
    total = vals.sum()
    if total <= 0 or not np.isfinite(freq) or freq <= 0:
        continue
    fracs = vals / total
    y_anchor = freq * PIE_GAP_FACTOR  

    axins = inset_axes(
        ax,
        width=pie_diam_in, height=pie_diam_in,  
        loc="center",
        bbox_to_anchor=(x, y_anchor),
        bbox_transform=ax.transData,
        borderpad=0.0
    )
    axins.pie(fracs, startangle=90)
    axins.set_aspect("equal")
    axins.set_xticks([])
    axins.set_yticks([])
    for spine in axins.spines.values():
        spine.set_visible(False)


legend_ax = inset_axes(ax, width=1.6, height=1.6, loc="upper left", borderpad=1.2)
legend_ax.pie([1, 1, 1, 1], labels=pie_labels, startangle=90)
legend_ax.set_aspect("equal")
legend_ax.set_title("Pie segments", fontsize=10)
legend_ax.set_xticks([]); legend_ax.set_yticks([])
for spine in legend_ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=300)
plt.savefig(OUT_PDF)
plt.show()
