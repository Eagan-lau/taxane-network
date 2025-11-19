import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set the aesthetic style of the plots
sns.set(style="whitegrid")

# Define file paths
input_file_csv = 'grouped_statistics.csv'      # Input CSV file
output_plot_scatter = 'adjusted_scatter_plot.png'  # Output scatter plot file

# Function to read data from CSV
def read_data(csv_path):
    try:
        data = pd.read_csv(csv_path)
        print(f"Successfully read CSV file: {csv_path}")
    except FileNotFoundError:
        print(f"Error: CSV file '{csv_path}' not found.")
        raise
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        raise
    return data

# Read the data
stats = read_data(input_file_csv)

# Preview the data
print("Data Preview:")
print(stats.head())

# Define required columns
required_columns = ['count', 'cv_length', 'EC(eggNOG)', 'CLEAN']

# Check if required columns exist
for col in required_columns:
    if col not in stats.columns:
        raise ValueError(f"Missing required column: '{col}'. Please check your input file.")

# Ensure 'count' and 'cv_length' are numeric
for col in ['count', 'cv_length']:
    if not pd.api.types.is_numeric_dtype(stats[col]):
        stats[col] = pd.to_numeric(stats[col], errors='coerce')
        print(f"Converted column '{col}' to numeric.")

# Check for missing values
missing_values = stats[required_columns].isnull().sum()
print("\nMissing Values:")
print(missing_values)

# Remove rows with missing values in required columns
stats_clean = stats.dropna(subset=['count', 'cv_length'])
print(f"\nNumber of rows after cleaning: {stats_clean.shape[0]} (Original: {stats.shape[0]})")

# Remove rows with non-positive values in 'count' (since log scale cannot handle them)
stats_clean = stats_clean[(stats_clean['count'] > 0)]
print(f"Number of rows after removing non-positive 'count' values: {stats_clean.shape[0]}")

# Final check
if stats_clean.empty:
    raise ValueError("No data available for plotting after cleaning. Please check your data.")

# Create a new 'Category' column based on 'EC(eggNOG)' and 'CLEAN'
def categorize_row(row):
    ec_known = row['EC(eggNOG)'] == 'known'
    clean_known = row['CLEAN'] == 'known'
    if ec_known and not clean_known:
        return 'EC known only'
    elif not ec_known and clean_known:
        return 'CLEAN known only'
    elif ec_known and clean_known:
        return 'Both known'
    else:
        return 'Others'

stats_clean['Category'] = stats_clean.apply(categorize_row, axis=1)

# Define color palette
palette = {
    'EC known only': 'blue',
    'CLEAN known only': 'green',
    'Both known': 'red',
    'Others': 'gray'
}

# Create a JointGrid with marginal distributions
g = sns.JointGrid(data=stats_clean, x='count', y='cv_length', height=8, ratio=5, marginal_ticks=True)

# Scatter plot with hue based on 'Category'
g.plot_joint(sns.scatterplot, hue='Category', data=stats_clean, palette=palette, s=150, alpha=0.7)


g.ax_joint.set_yscale('linear')
g.ax_joint.set_xscale('log')

# Customize marginal histograms
sns.histplot(
    data=stats_clean, x='count', bins=50, ax=g.ax_marg_x, color='blue', edgecolor='blue', alpha=0.6
)
sns.histplot(
    data=stats_clean, y='cv_length', bins=30, ax=g.ax_marg_y, color='blue', edgecolor='blue', alpha=0.6
)


g.ax_marg_x.set_yscale('log')


g.ax_marg_y.set_xscale('log')

# Customize axes
g.ax_joint.set_title('Adjusted: Group Size vs. CV of Gene Lengths', fontsize=16)
g.ax_joint.set_xlabel('Gene Count per Group (Log Scale)', fontsize=14)
g.ax_joint.set_ylabel('Coefficient of Variation of Gene Lengths (Linear Scale)', fontsize=14)

# Tweak layout
plt.tight_layout()

# Save the adjusted plot
plt.savefig(output_plot_scatter, dpi=300)
print(f"\nAdjusted scatter plot saved as '{output_plot_scatter}'.")

# Show plot
plt.show()
