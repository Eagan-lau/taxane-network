#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Read picture.xlsx:
- First row: header
- First column: species (x-axis)
- Columns 3–6: bar plots (percent values)
- Column 2: line plot (percent values)
- Sort x-axis by column 2 descending
- Y-axis range: 0–100
- Font: Arial
- Cleans cells with spaces, '%' etc.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# ===== Font settings =====
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42

# ===== File paths =====
INPUT_XLSX = "picture.xlsx"
OUT_PNG = "picture_plot.png"
OUT_PDF = "picture_plot.pdf"

# ===== Read and clean data =====
df = pd.read_excel(INPUT_XLSX)

# Function to clean numeric cells
def to_float(x):
    if pd.isna(x):
        return np.nan
    x = str(x).replace("%", "").replace(",", "").strip()
    x = x.replace(" ", "")
    try:
        return float(x)
    except ValueError:
        return np.nan

# Apply cleaning to all numeric columns
for c in df.columns[1:]:
    df[c] = df[c].apply(to_float)

# Sort by column 2 (descending)
df = df.sort_values(by=df.columns[1], ascending=False)

# X-axis: species (column 1)
species = df.iloc[:, 0].astype(str)

# Line: column 2
line_data = df.iloc[:, 1]

# Bars: columns 3–6
bar_cols = df.columns[2:6]
bar_data = df.iloc[:, 2:6]

# ===== Plot =====
x = np.arange(len(species))
bar_width = 0.18

fig, ax1 = plt.subplots(figsize=(10, 6))

# Bar plots
for i, col in enumerate(bar_cols):
    ax1.bar(x + (i - 1.5) * bar_width, bar_data[col], width=bar_width, label=col)

# Line plot
ax1.plot(x, line_data, color="black", marker="o", linewidth=2, label=df.columns[1])

# Axis settings
ax1.set_xlabel("Species")
ax1.set_ylabel("Percentage (%)")
ax1.set_xticks(x)
ax1.set_xticklabels(species, rotation=45, ha="right")
ax1.set_ylim(0, 100)

# Legend
ax1.legend(frameon=False)

plt.tight_layout()

# Save
plt.savefig(OUT_PNG, dpi=300)
plt.savefig(OUT_PDF)
plt.show()
