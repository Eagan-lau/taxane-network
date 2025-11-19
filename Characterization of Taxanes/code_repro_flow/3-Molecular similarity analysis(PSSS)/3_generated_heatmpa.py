import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


excel_file = 'global_difference_matrix_difference.xlsx'
df = pd.read_excel(excel_file, index_col=0)


distance_matrix = df.values


condensed_distance_matrix = squareform(distance_matrix)


row_linkage = linkage(condensed_distance_matrix, method='average')
col_linkage = linkage(condensed_distance_matrix, method='average')


sns.set(context="poster", font_scale=0.6)
g = sns.clustermap(
    df,
    row_linkage=row_linkage,
    col_linkage=col_linkage,
    cmap="coolwarm",
    linewidths=0,
    figsize=(15, 10),
    cbar_kws={'label': 'Collection Difference'}
)


plt.savefig('heatmap_UPGMA.pdf', bbox_inches='tight')
plt.close()

print("save 'heatmap_UPGMA.pdf'")
