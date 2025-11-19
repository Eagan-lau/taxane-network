import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from matplotlib.colors import LinearSegmentedColormap

# 1. Read the processed symmetric matrix file
matrix_file = 'symmetric_matrix_avg.csv'
matrix = pd.read_csv(matrix_file, sep='\t', index_col=0)

print("Symmetric Matrix loaded successfully:")
print(matrix.head())

# 2. Perform hierarchical clustering
Z = linkage(matrix, method='average', optimal_ordering=True)

# Save row order
row_order = leaves_list(Z)
pd.DataFrame(row_order).to_csv("row_order.csv", index=False)

# 3. Create custom colormap (0=blue, 0.6=white, 1=red)
cmap = LinearSegmentedColormap.from_list(
    'custom_cmap',
    [(0.0, 'blue'),
     (0.7, 'white'),
     (1.0, 'red')],
    N=256
)

# 4. Plot clustered heatmap
sns.clustermap(matrix, cmap=cmap, row_linkage=Z, col_linkage=Z, vmin=0, vmax=1)

# 5. Show figure
plt.title('Clustered Heatmap of Symmetric Protein Comparisons')
plt.show()
