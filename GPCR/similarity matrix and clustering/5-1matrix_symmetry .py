import pandas as pd


file_path = 'protein_comparison_matrix_CYP450.txt'  
matrix = pd.read_csv(file_path, sep='\t', index_col=0)


symmetric_matrix_avg = (matrix + matrix.T) / 2


output_path = 'symmetric_matrix_avg.csv'
symmetric_matrix_avg.to_csv(output_path, sep='\t')

print(f"Symmetric matrix saved to: {output_path}")

