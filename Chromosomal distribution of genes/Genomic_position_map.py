#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Genome-wide gene position visualization for diploid species.
All gene families are plotted in the same figure, one family per row.
For each family row:
- hap1 and hap2 chromosomes are drawn side-by-side.
- A vertical Mb scale bar is shown on the left (optional).
- Gene positions follow the same orientation as the scale bar:
  top represents 0 bp and bottom represents maximum chromosome length.
The script also exports per-family per-chromosome gene counts to Excel.
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse, Arc

# ===================== Configuration =====================
# Input files for hap1
LOC_FILE_H1 = "hap1_gene_loc.txt"
FAI_FILE_H1 = "Taxus_chinensis_var_mairei_hap1_N_20230809.fasta.seqkit.fai"

# Input files for hap2
LOC_FILE_H2 = "hap2_gene_loc.txt"
FAI_FILE_H2 = "Taxus_chinensis_var_mairei_hap2_N_20230809.fasta.seqkit.fai"

# Family annotation file containing both haplotypes
TOP_FILE = "all-90.csv"   # expected columns: protein_ID, family, ...

# Output files
OUTPUT_FIG = "genome_gene_positions.pdf"
OUTPUT_XLSX = "family_chrom_counts.xlsx"

# Figure layout parameters
BASE_FIG_W_PER_BLOCK = 12.0
BASE_FIG_H_PER_ROW = 2.6
TOP_MARGIN = 0.35
BOTTOM_MARGIN = 0.35
LEFT_MARGIN_X = 0.7
RIGHT_MARGIN_X = 1.4

BACKGROUND = "white"

# Chromosome layout
N_CHROM = 12
PAIR_GAP = 1.0         # gap between chromosome pairs
SISTER_GAP = 0.45      # gap between hap1/hap2 within a pair
CHROM_WIDTH = 0.33
ROW_HEIGHT = 1.85
CAP_ASPECT = 1.0
EDGE_COLOR = "#444444"
CHROM_FACE = "#B5B5B5"

# Gene position tick-line style
TICK_LINE_LEN_FACTOR = 0.90
TICK_LINE_WIDTH = 1.6
TICK_CAPSTYLE = "round"

# Gene-count labels above chromosomes
SHOW_ZERO_COUNTS = False
COUNT_FONT_SIZE = 8.2
COUNT_COLOR_MODE = "family"  # "family" / "black" / custom color
COUNT_PAD_ABOVE = 0.10

# Chromosome label styles
LABEL_COLOR = "#111111"
CHR_LABEL_FONTSIZE = 8.2
CHR_LABEL_EXTRA_PAD = 0.12
FAMILY_LABEL_FONTSIZE = 10.0

# Predefined colors for families (optional)
FAMILY_COLOR_PRESET = {}
FALLBACK_COLORS = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756", "#72B7B2",
    "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC", "#FFA600"
]

# Prefix indicating hap1 proteins in TOP_FILE
HAP1_PREFIX = "Tca"

# Scale bar configuration
SHOW_SCALE_BAR = True
SCALE_LEFT_OFFSET = 1.05   # horizontal shift relative to LEFT_MARGIN_X
SCALE_LINE_WIDTH = 2.0
SCALE_TICK_LEN = 0.10
SCALE_TICK_WIDTH = 1.8
SCALE_FONT_SIZE = 9.0
SCALE_COLOR = "#000000"
SCALE_ARROW_SIZE = 6
SCALE_TICKS_MIN = 6
SCALE_TICKS_MAX = 12

DPI = 300
# =========================================================


def read_locations(path):
    """Read gene location file: GeneID, Chr, Start, End, Strand."""
    return pd.read_csv(
        path, sep="\t", header=None, engine="python",
        names=["GeneID", "Chr", "Start", "End", "Strand"],
        dtype={"GeneID": str, "Chr": str, "Start": int, "End": int, "Strand": str}
    )


def read_fai(path):
    """Read FASTA index (.fai) file and extract chromosome lengths."""
    fai = pd.read_csv(
        path, sep="\t", header=None, usecols=[0, 1],
        names=["Chr", "Length"], engine="python"
    )
    fai["Length"] = pd.to_numeric(fai["Length"], errors="coerce").astype("Int64")
    fai = fai.dropna(subset=["Length"]).astype({"Length": int})
    return fai


def read_top_with_hap(path):
    """
    Read family annotation file and infer haplotype from protein_ID prefix.
    Hap assignment rule:
      - If protein_ID begins with HAP1_PREFIX → hap1
      - Otherwise → hap2
    """
    t = pd.read_csv(path)
    t = t[["protein_ID", "family"]].copy()
    t["hap"] = t["protein_ID"].astype(str).apply(
        lambda x: "hap1" if x.startswith(HAP1_PREFIX) else "hap2"
    )
    t = t.rename(columns={"protein_ID": "GeneID", "family": "Family"}).reset_index(drop=True)
    return t


def pick_color_for_family(fam, idx=[0]):
    """Return a stable color for each family, using presets if available."""
    if fam in FAMILY_COLOR_PRESET:
        return FAMILY_COLOR_PRESET[fam]
    c = FALLBACK_COLORS[idx[0] % len(FALLBACK_COLORS)]
    idx[0] += 1
    return c


def draw_capsule(ax, cx, by, h, w, face, edge, cap_aspect=1.0):
    """
    Draw a chromosome as a capsule (rounded top and bottom).
    cx = center x, by = bottom y, h = total height, w = width.
    """
    cap_h = w * cap_aspect
    body_y = by + cap_h / 2.0
    body_h = max(0.0, h - cap_h)

    ax.add_patch(Rectangle((cx - w/2, body_y), w, body_h, lw=0, fc=face, zorder=1))
    bcy = by + cap_h/2.0
    ax.add_patch(Ellipse((cx, bcy), w, cap_h, lw=0, fc=face, zorder=1))
    ax.add_patch(Arc((cx, bcy), w, cap_h, theta1=180, theta2=360, lw=1.0, ec=edge, zorder=2))

    tcy = by + h - cap_h/2.0
    ax.add_patch(Ellipse((cx, tcy), w, cap_h, lw=0, fc=face, zorder=1))
    ax.add_patch(Arc((cx, tcy), w, cap_h, theta1=0, theta2=180, lw=1.0, ec=edge, zorder=2))

    lx, rx = cx - w/2, cx + w/2
    ax.plot([lx, lx], [body_y, body_y + body_h], color=edge, lw=1.0, zorder=2)
    ax.plot([rx, rx], [body_y, body_y + body_h], color=edge, lw=1.0, zorder=2)


def chr_index(chr_name):
    """Extract chromosome numeric index from names like 'TcChr01a'."""
    m = re.search(r"TcChr(\d+)[ab]?$", str(chr_name))
    if not m:
        return None
    return int(m.group(1))


def build_chr_lists_by_index(fai_h1, fai_h2):
    """
    Build ordered chromosome lists (1..N) based on FAI files.
    Returns:
      - list of indices present in either hap
      - mapping index → hap1 chromosome name
      - mapping index → hap2 chromosome name
    """
    map_h1 = {}
    map_h2 = {}
    for c in fai_h1["Chr"].unique():
        idx = chr_index(c)
        if idx is not None:
            map_h1[idx] = c
    for c in fai_h2["Chr"].unique():
        idx = chr_index(c)
        if idx is not None:
            map_h2[idx] = c
    all_idx = sorted(set(map_h1.keys()) | set(map_h2.keys()))
    return all_idx, map_h1, map_h2


def layout_x_pairs(indices, map_h1, map_h2):
    """
    Determine x-locations for chromosome pairs across hap1 and hap2.
    Ensures hap1 and hap2 chromosomes of the same index are drawn side-by-side.
    """
    xs = {}
    x = LEFT_MARGIN_X
    for idx in indices:
        c1 = map_h1.get(idx)
        c2 = map_h2.get(idx)
        if c1:
            xs[("hap1", c1)] = x
            if c2:
                xs[("hap2", c2)] = x + SISTER_GAP
                x += SISTER_GAP
        if (not c1) and c2:
            xs[("hap2", c2)] = x
        x += PAIR_GAP
    return xs


def pick_mb_step(max_bp, ticks_min=6, ticks_max=12):
    """
    Select a suitable Mb step for scale bar ticks so that the number of ticks
    falls into a reasonable range (ticks_min to ticks_max).
    """
    if max_bp <= 0:
        return 100e6

    max_mb = max_bp / 1e6

    candidates = np.array([
        5, 10, 12, 20, 25, 30, 40, 50, 60, 75, 80, 90, 100,
        108, 120, 150, 180, 200, 250, 300, 360, 400, 450, 500,
        600, 750, 800, 900, 1000
    ], dtype=float)

    if max_mb > candidates[-1]:
        scale = 10 ** np.floor(np.log10(max_mb))
        for base in [1, 2, 2.5, 5]:
            candidates = np.append(candidates, base * scale)

    best, best_diff = candidates[0], 1e9
    for step in candidates:
        n = max_mb / step
        if n < ticks_min:
            diff = ticks_min - n
        elif n > ticks_max:
            diff = n - ticks_max
        else:
            diff = abs(n - round(n))
        if diff < best_diff:
            best, best_diff = step, diff

    return best * 1e6


def draw_scale_bar(ax, x_scale, row_bottom, row_height, max_bp_ref,
                   color=SCALE_COLOR):
    """
    Draw a vertical Mb scale bar for each family row.
    The top of the bar corresponds to 0 bp; the bottom corresponds to max_bp_ref.
    Tick labels decrease from top to bottom.
    """
    ax.plot(
        [x_scale, x_scale],
        [row_bottom, row_bottom + row_height],
        color=color, lw=SCALE_LINE_WIDTH, solid_capstyle="butt", zorder=5
    )

    step_bp = pick_mb_step(max_bp_ref, SCALE_TICKS_MIN, SCALE_TICKS_MAX)
    ticks_bp = np.arange(0, max_bp_ref + step_bp * 0.5, step_bp)
    labels = [f"{int(round(bp / 1e6))} Mb" for bp in ticks_bp]

    scale = row_height / max_bp_ref if max_bp_ref > 0 else 0
    for bp, lab in zip(ticks_bp, labels):
        y = row_bottom + (max_bp_ref - bp) * scale
        ax.plot([x_scale - SCALE_TICK_LEN, x_scale],
                [y, y], color=color, lw=SCALE_TICK_WIDTH, zorder=5)
        ax.text(x_scale - SCALE_TICK_LEN - 0.02, y, lab,
                ha="right", va="center",
                color=color, fontsize=SCALE_FONT_SIZE)

    ax.annotate(
        "",
        xy=(x_scale, row_bottom - 0.02),
        xytext=(x_scale, row_bottom + 0.06),
        arrowprops=dict(
            arrowstyle='-|>',
            lw=SCALE_LINE_WIDTH,
            color=color,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=SCALE_ARROW_SIZE
        )
    )


def main():
    # ---------------- Load input data ----------------
    df_loc_h1 = read_locations(LOC_FILE_H1)
    df_loc_h2 = read_locations(LOC_FILE_H2)
    df_fai_h1 = read_fai(FAI_FILE_H1)
    df_fai_h2 = read_fai(FAI_FILE_H2)
    df_top = read_top_with_hap(TOP_FILE)

    # Assign haplotypes based on naming rule
    top_h1 = df_top[df_top["hap"] == "hap1"][["GeneID", "Family"]].copy()
    top_h2 = df_top[df_top["hap"] == "hap2"][["GeneID", "Family"]].copy()

    # Merge with gene location files
    m_h1 = top_h1.merge(df_loc_h1, on="GeneID", how="left")
    m_h2 = top_h2.merge(df_loc_h2, on="GeneID", how="left")

    # Clean and compute midpoints
    def cleanup(m):
        m = m.dropna(subset=["Chr", "Start", "End"]).copy()
        m["Start"] = pd.to_numeric(m["Start"], errors="coerce")
        m["End"] = pd.to_numeric(m["End"], errors="coerce")
        m = m.dropna(subset=["Start", "End"]).copy()
        m["Mid"] = (m["Start"] + m["End"]) / 2.0
        return m

    m_h1 = cleanup(m_h1)
    m_h2 = cleanup(m_h2)

    # Keep only chromosomes present in FAI files
    m_h1 = m_h1[m_h1["Chr"].isin(df_fai_h1["Chr"].unique())].copy()
    m_h2 = m_h2[m_h2["Chr"].isin(df_fai_h2["Chr"].unique())].copy()

    # Determine family list
    families = sorted(
        set(m_h1["Family"].dropna().unique()) |
        set(m_h2["Family"].dropna().unique())
    )
    if not families:
        print("No families found.")
        return

    # Build ordered chromosome index list
    indices, map_h1, map_h2 = build_chr_lists_by_index(df_fai_h1, df_fai_h2)
    if not indices:
        print("No recognizable chromosome IDs found in FAI files.")
        return

    # ---------------- Figure setup ----------------
    fig_w = BASE_FIG_W_PER_BLOCK
    fig_h = TOP_MARGIN + BOTTOM_MARGIN + BASE_FIG_H_PER_ROW * len(families)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor(BACKGROUND)
    ax.set_facecolor(BACKGROUND)

    xs = layout_x_pairs(indices, map_h1, map_h2)

    # Height scales
    max_len_h1 = df_fai_h1["Length"].max() if len(df_fai_h1) else 0
    max_len_h2 = df_fai_h2["Length"].max() if len(df_fai_h2) else 0
    scale_h1 = (ROW_HEIGHT / max_len_h1) if max_len_h1 > 0 else 1.0
    scale_h2 = (ROW_HEIGHT / max_len_h2) if max_len_h2 > 0 else 1.0

    len_h1_map = dict(zip(df_fai_h1["Chr"], df_fai_h1["Length"]))
    len_h2_map = dict(zip(df_fai_h2["Chr"], df_fai_h2["Length"]))
    max_bp_ref = max(max_len_h1, max_len_h2)

    min_y, max_y = 1e9, -1e9

    # Assign colors to families
    used_idx = [0]
    fam2color = {f: pick_color_for_family(f, used_idx) for f in families}

    # Accumulate count table
    count_rows = []

    # Starting baseline
    top_y = TOP_MARGIN + BASE_FIG_H_PER_ROW * len(families)

    # ---------------- Per-family plotting ----------------
    for r, fam in enumerate(families):
        row_bottom = top_y - (r + 1) * BASE_FIG_H_PER_ROW + 0.15

        f_h1 = m_h1[m_h1["Family"] == fam].copy()
        f_h2 = m_h2[m_h2["Family"] == fam].copy()

        # Left-side Mb scale bar
        if SHOW_SCALE_BAR and max_bp_ref > 0:
            x_scale = LEFT_MARGIN_X - SCALE_LEFT_OFFSET
            draw_scale_bar(ax, x_scale, row_bottom, ROW_HEIGHT, max_bp_ref, color=SCALE_COLOR)

        # Family label on the left
        ax.text(
            LEFT_MARGIN_X - 0.55, row_bottom + ROW_HEIGHT / 2,
            str(fam),
            ha="right", va="center",
            fontsize=FAMILY_LABEL_FONTSIZE,
            color=LABEL_COLOR, weight="bold"
        )

        # Gene counts for each hap × chromosome
        cnt_h1 = f_h1.groupby("Chr").size().to_dict()
        cnt_h2 = f_h2.groupby("Chr").size().to_dict()

        # ---- Draw chromosomes (hap1 and hap2) ----
        for idx in indices:
            c1 = map_h1.get(idx)
            c2 = map_h2.get(idx)

            if c1:
                L1 = int(len_h1_map.get(c1, 0))
                h1 = L1 * scale_h1
                x1 = xs[("hap1", c1)]
                draw_capsule(ax, x1, row_bottom, h1, CHROM_WIDTH, CHROM_FACE, EDGE_COLOR, CAP_ASPECT)

                n1 = int(cnt_h1.get(c1, 0))
                count_rows.append({"Family": fam, "Hap": "hap1", "Chr": c1, "ChrIndex": idx, "Count": n1})

                if SHOW_ZERO_COUNTS or n1 > 0:
                    count_color = fam2color[fam] if COUNT_COLOR_MODE == "family" else (COUNT_COLOR_MODE or "black")
                    y_count = row_bottom + h1 + COUNT_PAD_ABOVE
                    ax.text(x1, y_count, str(n1),
                            ha="center", va="bottom",
                            fontsize=COUNT_FONT_SIZE, color=count_color)
                    max_y = max(max_y, y_count)

                if r == 0:
                    y_label = row_bottom + h1 + COUNT_PAD_ABOVE + CHR_LABEL_EXTRA_PAD
                    ax.text(x1, y_label, c1,
                            ha="center", va="bottom",
                            fontsize=CHR_LABEL_FONTSIZE, color=LABEL_COLOR)
                    max_y = max(max_y, y_label)

                min_y = min(min_y, row_bottom)
                max_y = max(max_y, row_bottom + h1)

            if c2:
                L2 = int(len_h2_map.get(c2, 0))
                h2 = L2 * scale_h2
                x2 = xs[("hap2", c2)]
                draw_capsule(ax, x2, row_bottom, h2, CHROM_WIDTH, CHROM_FACE, EDGE_COLOR, CAP_ASPECT)

                n2 = int(cnt_h2.get(c2, 0))
                count_rows.append({"Family": fam, "Hap": "hap2", "Chr": c2, "ChrIndex": idx, "Count": n2})

                if SHOW_ZERO_COUNTS or n2 > 0:
                    count_color = fam2color[fam] if COUNT_COLOR_MODE == "family" else (COUNT_COLOR_MODE or "black")
                    y_count = row_bottom + h2 + COUNT_PAD_ABOVE
                    ax.text(x2, y_count, str(n2),
                            ha="center", va="bottom",
                            fontsize=COUNT_FONT_SIZE, color=count_color)
                    max_y = max(max_y, y_count)

                if r == 0:
                    y_label = row_bottom + h2 + COUNT_PAD_ABOVE + CHR_LABEL_EXTRA_PAD
                    ax.text(x2, y_label, c2,
                            ha="center", va="bottom",
                            fontsize=CHR_LABEL_FONTSIZE, color=LABEL_COLOR)
                    max_y = max(max_y, y_label)

                min_y = min(min_y, row_bottom)
                max_y = max(max_y, row_bottom + h2)

        # ---- Draw gene tick marks (hap1) ----
        half = CHROM_WIDTH * TICK_LINE_LEN_FACTOR / 2.0
        col = fam2color[fam]

        for _, row in f_h1.iterrows():
            c = row["Chr"]
            key = ("hap1", c)
            if key not in xs:
                continue
            x = xs[key]
            L = len_h1_map.get(c, 0)
            h = L * scale_h1
            y = row_bottom + h - float(row["Mid"]) * scale_h1
            ax.plot([x - half, x + half], [y, y],
                    color=col, lw=TICK_LINE_WIDTH, solid_capstyle=TICK_CAPSTYLE, zorder=3)
            min_y = min(min_y, y - 0.1)
            max_y = max(max_y, y + 0.1)

        # ---- Draw gene tick marks (hap2) ----
        for _, row in f_h2.iterrows():
            c = row["Chr"]
            key = ("hap2", c)
            if key not in xs:
                continue
            x = xs[key]
            L = len_h2_map.get(c, 0)
            h = L * scale_h2
            y = row_bottom + h - float(row["Mid"]) * scale_h2
            ax.plot([x - half, x + half], [y, y],
                    color=col, lw=TICK_LINE_WIDTH, solid_capstyle=TICK_CAPSTYLE, zorder=3)
            min_y = min(min_y, y - 0.1)
            max_y = max(max_y, y + 0.1)

    # ---------------- Finalize plot ----------------
    x_min = (LEFT_MARGIN_X - max(SCALE_LEFT_OFFSET, 0.75) - 0.35) if SHOW_SCALE_BAR else (LEFT_MARGIN_X - 0.75)
    x_max = max(xs.values()) + RIGHT_MARGIN_X if xs else LEFT_MARGIN_X + 5
    ax.set_xlim(x_min, x_max)

    pad = 0.35
    ax.set_ylim(max(0, min_y - pad), max_y + pad)

    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    plt.tight_layout()
    plt.savefig(
        OUTPUT_FIG,
        dpi=DPI,
        bbox_inches="tight",
        pad_inches=0.15,
        facecolor=fig.get_facecolor()
    )
    plt.close()
    print("Figure saved:", OUTPUT_FIG)

    # ---------------- Export gene counts to Excel ----------------
    if count_rows:
        df_counts = pd.DataFrame(count_rows)

        hap1_wide = (
            df_counts[df_counts["Hap"] == "hap1"]
            .pivot_table(index="Family", columns="Chr", values="Count", aggfunc="sum")
            .fillna(0).astype(int)
            .sort_index(axis=1, key=lambda c: [chr_index(x) or 999 for x in c])
        )

        hap2_wide = (
            df_counts[df_counts["Hap"] == "hap2"]
            .pivot_table(index="Family", columns="Chr", values="Count", aggfunc="sum")
            .fillna(0).astype(int)
            .sort_index(axis=1, key=lambda c: [chr_index(x) or 999 for x in c])
        )

        with pd.ExcelWriter(OUTPUT_XLSX, engine="xlsxwriter") as writer:
            df_counts.sort_values(["Family", "Hap", "ChrIndex", "Chr"]) \
                .to_excel(writer, index=False, sheet_name="tidy_counts")
            hap1_wide.to_excel(writer, sheet_name="hap1_wide")
            hap2_wide.to_excel(writer, sheet_name="hap2_wide")

        print("Gene counts exported:", OUTPUT_XLSX)
    else:
        print("No gene counts available for export.")


if __name__ == "__main__":
    main()
