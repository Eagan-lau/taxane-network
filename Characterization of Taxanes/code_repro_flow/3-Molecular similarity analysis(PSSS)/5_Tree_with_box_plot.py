import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import re
import mpld3
import numpy as np
from Bio import Phylo
import warnings

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)

# 1. Read the .nwk file
nwk_file = "UPGMA_tree_iTOL.nwk"
tree = Phylo.read(nwk_file, "newick")

# Get the leaf node labels of the tree
raw_leaf_names = [terminal.name for terminal in tree.get_terminals()]
leaf_names = [str(terminal.name).strip("'").strip().lower() for terminal in tree.get_terminals()]
print("Raw tree leaf names:", raw_leaf_names)
print("Tree leaf names:", leaf_names)

# 2. Read and normalize the Excel file
df = pd.read_excel('0_processed_final.xlsx', header=0)
raw_excel_names = df.iloc[:, 0].tolist()
df.iloc[:, 0] = df.iloc[:, 0].astype(str).str.strip().str.lower()
excel_names = df.iloc[:, 0].tolist()

# Check name matching
missing_in_excel = [name for name in leaf_names if name not in excel_names]
missing_in_tree = [name for name in excel_names if name not in leaf_names]
print("Names in tree but not in Excel:", missing_in_excel)
print("Names in Excel but not in tree:", missing_in_tree)

# Ensure the table row order matches the tree leaf nodes
df = df.set_index(df.columns[0]).reindex(leaf_names).reset_index()

# 3. Group into R groups
r_groups = {}
for col in df.columns:
    match = re.match(r'(R\d+)', col)
    if match:
        r_name = match.group(1)
        if r_name not in r_groups:
            r_groups[r_name] = []
        r_groups[r_name].append(col)

r_order = [f'R{i}' for i in range(1, 16)]
r_groups_ordered = {r: r_groups[r] for r in r_order if r in r_groups}

# 4. Read the color mapping
color_df = pd.read_excel('element_color_mapping.xlsx', header=0)
color_dict = dict(zip(color_df.iloc[:, 1], color_df.iloc[:, 3]))
label_dict = dict(zip(color_df.iloc[:, 1], color_df.iloc[:, 0]))
color_dict['**'] = 'none'

# 5. Set shapes
shape_dict = {r: 'rect' for r in r_groups_ordered.keys()}

# 6. Calculate X-axis positions for boxes
n_groups = len(r_groups_ordered)
m = df.shape[0]
max_width = 7.27
max_height = 10.69

total_cols = sum(len(cols) for cols in r_groups_ordered.values())
w = min(0.3, max_width / (total_cols + n_groups * 0.5))
h = min(0.5, max_height / (m + 2))
s_sub = w * 0.1
s_h = w * 0.5
s_v = h * 0.1

x = 0
column_centers = {}
r_label_positions = {}
for r in r_groups_ordered.keys():
    cols = r_groups_ordered[r]
    n = len(cols)
    for j, col in enumerate(cols):
        center_x = x + j * (w + s_sub) + w / 2
        column_centers[col] = center_x
    r_label_positions[r] = x + (n * w + (n - 1) * s_sub) / 2
    x += n * w + (n - 1) * s_sub + s_h

# 7. Calculate the number of color types per column and row
col_color_counts = {col: len(set(df[col].dropna().unique())) for col in column_centers.keys()}
row_color_counts = [len(set(df.iloc[i, 1:].dropna().unique())) for i in range(m)]

# 8. Create the main figure
fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]}, figsize=(11.69, 8.27))

# 9. Draw the tree on the left and assign Y coordinates
Phylo.draw(tree, axes=ax1, do_show=False)
for text in ax1.texts:
    if text.get_text() in raw_leaf_names:
        text.set_visible(False)
ax1.set_title("Phylogenetic Tree")

# Manually assign Y coordinates from the tree structure
y_coords = {}
leaf_count = len(leaf_names)
y_step = max_height / (leaf_count + 1)
for i, terminal in enumerate(tree.get_terminals()):
    y_coords[leaf_names[i]] = -i * y_step
y_centers = [y_coords[name] for name in leaf_names]

# 10. Draw the box plot on the right
ax2.set_xlim(0, x + w)
ax2.set_ylim(min(y_centers) - h, max(y_centers) + h)
for i in range(m):
    for r in r_groups_ordered.keys():
        cols = r_groups_ordered[r]
        for col in cols:
            val = df[col].iloc[i]
            center_x = column_centers[col]
            center_y = y_centers[i]
            if pd.isna(val):
                continue
            else:
                color = color_dict.get(val, 'grey')
                rect = Rectangle((center_x - w / 2, center_y - h / 2), w, h, color=color, ec='black', lw=0.5)
                ax2.add_patch(rect)

# Add R column labels
for r, x_pos in r_label_positions.items():
    ax2.text(x_pos, min(y_centers) - h / 2, r, ha='center', va='center', fontsize=8, fontweight='bold')

# 11. Adjust subplot spacing
plt.subplots_adjust(left=0.1, right=0.85, top=0.85, bottom=0.1)

# 12. Add bar chart at the top
ax2_pos = ax2.get_position()
ax3 = fig.add_axes([ax2_pos.x0, 0.9, ax2_pos.width, 0.08])
ax3.bar(list(column_centers.values()), [col_color_counts[col] for col in column_centers.keys()], 
        width=w, align='center', color='lightblue')
ax3.set_xlim(0, x + w)
ax3.set_ylim(0, max(col_color_counts.values()) + 1)
ax3.set_title("Column Color Types", fontsize=8)
ax3.xaxis.set_visible(False)

# 13. Add bar chart on the right
ax4 = fig.add_axes([0.9, ax2_pos.y0, 0.08, ax2_pos.height])
ax4.barh(y_centers, row_color_counts, height=h, align='center', color='lightgreen')
ax4.set_ylim(min(y_centers) - h, max(y_centers) + h)
ax4.set_xlim(0, max(row_color_counts) + 1)
ax4.set_title("Row Color Types", fontsize=8, rotation=90, ha='center')
ax4.yaxis.set_visible(False)

# 14. Hide axes
ax1.axis('off')
ax2.axis('off')

# 15. Save the main figure as PDF (without legend)
plt.savefig('main_output_with_bars.pdf', format='pdf', bbox_inches='tight', pad_inches=0.1)
print("Main figure with bar charts saved as 'main_output_with_bars.pdf'")

# 16. Save the legend separately as PDF
fig_legend, ax_legend = plt.subplots(figsize=(2, 8))
legend_patches = []
legend_labels = []
for val, color in color_dict.items():
    label = label_dict.get(val, str(val))
    if color == 'none':  # Handle '**', display as an empty box
        patch = Rectangle((0, 0), 0.2, 0.2, color='white', ec='black', lw=0.5)
    else:
        patch = Rectangle((0, 0), 0.2, 0.2, color=color, ec='black', lw=0.5)
    legend_patches.append(patch)
    legend_labels.append(label)
ax_legend.legend(handles=legend_patches, labels=legend_labels, loc='center', title="Legend", fontsize=8)
ax_legend.axis('off')
plt.savefig('legend_output.pdf', format='pdf', bbox_inches='tight')
print("Legend saved as 'legend_output.pdf'")

# 17. HTML output (increase canvas width)
fig_html, (ax1_html, ax2_html) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]}, figsize=(18, 15))  # Increase width from 15 to 18
Phylo.draw(tree, axes=ax1_html, do_show=False, label_func=lambda x: x.name if x.name else '')
for text in ax1_html.texts:
    if text.get_text() in raw_leaf_names:
        text.set_fontsize(6)
ax1_html.set_title("Phylogenetic Tree")

# Use the same Y coordinates
y_centers_html = y_centers

ax2_html.set_xlim(0, x + w)
ax2_html.set_ylim(min(y_centers_html) - h, max(y_centers_html) + h)
for i in range(m):
    for r in r_groups_ordered.keys():
        cols = r_groups_ordered[r]
        for col in cols:
            val = df[col].iloc[i]
            center_x = column_centers[col]
            center_y = y_centers_html[i]
            if pd.isna(val):
                continue
            else:
                color = color_dict.get(val, 'grey')
                rect = Rectangle((center_x - w / 2, center_y - h / 2), w, h, color=color, ec='black', lw=0.5)
                ax2_html.add_patch(rect)

# Add R column labels
for r, x_pos in r_label_positions.items():
    ax2_html.text(x_pos, min(y_centers_html) - h / 2, r, ha='center', va='center', fontsize=8, fontweight='bold')

# HTML bar charts
ax2_html_pos = ax2_html.get_position()
ax3_html = fig_html.add_axes([ax2_html_pos.x0, 0.9, ax2_html_pos.width, 0.08])
ax3_html.bar(list(column_centers.values()), [col_color_counts[col] for col in column_centers.keys()], 
             width=w, align='center', color='lightblue')
ax3_html.set_xlim(0, x + w)
ax3_html.set_ylim(0, max(col_color_counts.values()) + 1)
ax3_html.set_title("Column Color Types", fontsize=8)
ax3_html.xaxis.set_visible(False)

ax4_html = fig_html.add_axes([0.9, ax2_html_pos.y0, 0.08, ax2_html_pos.height])
ax4_html.barh(y_centers_html, row_color_counts, height=h, align='center', color='lightgreen')
ax4_html.set_ylim(min(y_centers_html) - h, max(y_centers_html) + h)
ax4_html.set_xlim(0, max(row_color_counts) + 1)
ax4_html.set_title("Row Color Types", fontsize=8, rotation=90, ha='center')
ax4_html.yaxis.set_visible(False)

# Add legend to HTML
ax2_html.legend(handles=legend_patches, labels=legend_labels, loc='upper left', bbox_to_anchor=(1.1, 1), 
                title="Legend", fontsize=8, frameon=False)

ax1_html.axis('off')
ax2_html.axis('off')

html_str = mpld3.fig_to_html(fig_html)
with open('output_with_bars.html', 'w', encoding='utf-8') as f:
    f.write(html_str)
print("HTML figure with bar charts and legend saved as 'output_with_bars.html'")

plt.close('all')
