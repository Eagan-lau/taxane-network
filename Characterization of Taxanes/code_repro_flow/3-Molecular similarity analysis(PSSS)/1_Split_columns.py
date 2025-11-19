import pandas as pd

# Read the 0.xlsx file
df = pd.read_excel('0.xlsx')

original_columns = df.columns.tolist()

# 1. Insert a blank column between each pair of adjacent columns
new_columns = []
for i, col in enumerate(original_columns):
    new_columns.append(col) 
    if i < len(original_columns) - 1:  
        new_columns.append(col + '_blank')  

# Create a new DataFrame, initializing blank columns as NaN
new_df = pd.DataFrame(index=df.index)
for col in new_columns:
    if '_blank' in col:
        new_df[col] = pd.NA  
    else:
        new_df[col] = df[col] 

# 2. Split the original columns by '.' and fill into the adjacent blank column
for col in original_columns:
    blank_col = col + '_blank'
    if blank_col in new_df.columns:
        # Split the original column, assuming it splits into at most two parts
        split_data = new_df[col].astype(str).str.split('.', expand=True)
        new_df[col] = split_data[0]  
        new_df[blank_col] = split_data[1] if split_data.shape[1] > 1 else pd.NA  # Fill the blank column with the second part (if it exists)

# 3. Remove columns that are entirely blank (NaN) except for the header
final_columns = []
for col in new_df.columns:
    # Check if all cells except the header are NaN
    if new_df[col].iloc[1:].isna().all():  
        continue  
    else:
        final_columns.append(col)  

final_df = new_df[final_columns]

# 4. Save the preliminary processed result to an Excel file
final_df.to_excel('0_processed.xlsx', index=False)

# 5. Read the file again and remove '_blank' from column names
processed_df = pd.read_excel('0_processed.xlsx')
# Replace '_blank' with '' in column names
new_column_names = [col.replace('_blank', '') for col in processed_df.columns]
processed_df.columns = new_column_names

# 6. Save the final result
processed_df.to_excel('0_processed_final.xlsx', index=False)

print("Processing complete, final result saved as '0_processed_final.xlsx'")
