import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import pickle

# Read the similarity matrix from the Excel file
excel_file = "global_difference_matrix_difference.xlsx"
df = pd.read_excel(excel_file, index_col=0)

# Obtain the distance matrix
distance_matrix = df.values

# Convert to condensed format
condensed_distance_matrix = squareform(distance_matrix)

# Perform UPGMA clustering analysis
Z = linkage(condensed_distance_matrix, method='average')

# Plot and save as PDF (vertical display)
plt.figure(figsize=(8, 10))
dendrogram(Z, labels=df.index, orientation='left', leaf_rotation=0, leaf_font_size=12)
plt.title("UPGMA Phylogenetic Tree")
plt.ylabel("Ligand ID")
plt.xlabel("Distance")
plt.tight_layout()
plt.savefig("UPGMA_phylogenetic_tree.pdf")

# Save interactive HTML tree
fig = ff.create_dendrogram(distance_matrix, orientation='left', labels=df.index.tolist(),
                           linkagefun=lambda x: linkage(x, 'average'))
fig.update_layout(width=800, height=600)
fig.write_html("UPGMA_phylogenetic_tree.html")

# Save tree data (Pickle)
with open('UPGMA_tree.pkl', 'wb') as f:
    pickle.dump(Z, f)

# Convert to strict iTOL-compatible Newick format
def build_newick(node, leaf_names, parent_dist):
    if node.is_leaf():
        return f"{leaf_names[node.id]}:{max(parent_dist - node.dist, 0):.6f}"
    else:
        left = build_newick(node.get_left(), leaf_names, node.dist)
        right = build_newick(node.get_right(), leaf_names, node.dist)
        return f"({left},{right}):{max(parent_dist - node.dist, 0):.6f}"

def linkage_to_newick(Z, labels):
    tree = to_tree(Z, rd=False)
    return build_newick(tree, labels, tree.dist) + ";"

# Generate and save Newick file
newick_str = linkage_to_newick(Z, df.index.tolist())

with open("UPGMA_tree_iTOL.nwk", "w") as f:
    f.write(newick_str)

print("Files have been successfully generated and saved:")
print("- Vertical PDF: UPGMA_phylogenetic_tree.pdf")
print("- Interactive HTML: UPGMA_phylogenetic_tree.html")
print("- iTOL-compatible Newick file: UPGMA_tree_iTOL.nwk")
print("- Pickle data file: UPGMA_tree.pkl")
