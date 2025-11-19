import pandas as pd
import numpy as np


num_side_chains = int(input("Enter the number of side chain positions in the molecules of this group:"))


columns_to_parse = [f'R{i}' for i in range(1, num_side_chains + 1)]

# Load data
data_file = '0_processed_final.xlsx'  
df = pd.read_excel(data_file)

# Parse side chains
def parse_side_chains(cell):
    if pd.isna(cell):
        return set()
    return set(cell.split('.'))

parsed_data = df[columns_to_parse].applymap(parse_side_chains)

# Compute global difference matrix
n = len(parsed_data)
global_difference_matrix = np.zeros((n, n))

for pos in range(len(columns_to_parse)):
    position_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            set_i = parsed_data.iloc[i, pos]
            set_j = parsed_data.iloc[j, pos]
            diff_size = len(set_i.symmetric_difference(set_j))
            difference_value = np.ceil(diff_size / 2)
            position_matrix[i, j] = difference_value
    weight = (pos + 1)**2 * 0.001
    global_difference_matrix += weight * position_matrix

# Save global difference matrix to Excel
global_difference_df = pd.DataFrame(global_difference_matrix, index=df['Ligand_ID'], columns=df['Ligand_ID'])
global_difference_df.to_excel('global_difference_matrix_difference.xlsx', index=True)
