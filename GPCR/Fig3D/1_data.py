import pandas as pd

# Define file paths
input_file = '1.xlsx'           
output_file = 'grouped_statistics.csv'  

# Read Excel file
try:
    data = pd.read_excel(input_file)
    print("Excel file loaded successfully.")
except Exception as e:
    print(f"Error reading Excel file: {e}")
    exit()

# Preview data
print("Data preview:")
print(data.head())

# Check required columns
required_columns = ['cluster', 'length']
for col in required_columns:
    if col not in data.columns:
        print(f"Error: missing required column '{col}'. Please check the Excel file.")
        exit()

# Check data types
print("\nData types:")
print(data.dtypes)

# Ensure 'length' is numeric
if not pd.api.types.is_numeric_dtype(data['length']):
    print("Converting 'length' column to numeric.")
    data['length'] = pd.to_numeric(data['length'], errors='coerce')

# Check missing values
print("\nMissing values summary:")
print(data[required_columns].isnull().sum())

# Remove rows with missing values
data_clean = data.dropna(subset=required_columns)
print(f"\nNumber of rows after cleaning: {data_clean.shape[0]} (original: {data.shape[0]})")

# Group and compute statistics
grouped = data_clean.groupby('cluster')['length']

# Compute statistical metrics
stats = grouped.agg(
    count='count',
    mean_length='mean',
    median_length='median',
    min_length='min',
    max_length='max',
    std_length='std',
    var_length='var'
).reset_index()

# Compute range
stats['range_length'] = stats['max_length'] - stats['min_length']

# Coefficient of variation (optional)
stats['cv_length'] = stats['std_length'] / stats['mean_length']

# Preview statistics
print("\nStatistics preview:")
print(stats.head())

# Save to CSV
try:
    stats.to_csv(output_file, index=False)
    print(f"\nStatistics saved to '{output_file}'.")
except Exception as e:
    print(f"Error saving CSV: {e}")

# Optional: save to Excel
output_excel = 'grouped_statistics.xlsx'
try:
    stats.to_excel(output_excel, index=False)
    print(f"Statistics saved to '{output_excel}'.")
except Exception as e:
    print(f"Error saving Excel file: {e}")
