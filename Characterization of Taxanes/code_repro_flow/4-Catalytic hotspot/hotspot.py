import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, to_hex, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

# Hill number calculation function
def calculate_hill_numbers(counts, q_values):
    total = sum(counts)
    if total == 0:
        return [0] * len(q_values)
    p_i = np.array(counts) / total
    hill_numbers = []
    for q in q_values:
        if q == 1:
            shannon = -sum(p * np.log(p) for p in p_i if p > 0)
            hill = np.exp(shannon)
        else:
            p_i_nonzero = p_i[p_i > 0]
            hill = np.power(np.sum(np.power(p_i_nonzero, q)), 1/(1-q)) if p_i_nonzero.size > 0 else 0
        hill_numbers.append(hill)
    return hill_numbers

# User input for the number of side chains
try:
    num_side_chains = int(input("Enter the number of side chains of this core molecule:"))
    if num_side_chains < 1:
        raise ValueError("The number of side chains must be a positive integer")
except ValueError as e:
    print(f"Input error: {e}")
    exit()

# Automatically generate column names from R1 to R<num_side_chains>
columns_to_parse = [f'R{i}' for i in range(1, num_side_chains + 1)]

# Load data
data_file = '0.xlsx'
try:
    df = pd.read_excel(data_file)
except FileNotFoundError:
    print(f"Error: file {data_file} not found")
    exit()

# Check required columns
if not all(col in df.columns for col in columns_to_parse):
    print(f"Error: data missing required columns {columns_to_parse}")
    exit()

# Function to parse side chains
def parse_side_chains(cell):
    if pd.isna(cell):
        return []
    return cell.split('.')

# Parse cell content for each column
parsed_data = df[columns_to_parse].apply(lambda x: x.map(parse_side_chains))

# Define colors and line styles for R1 to R16 (remove yellow series)
colors_16 = [
    '#FF0000', '#0000FF', '#800080', '#FF00FF', '#4B0082', '#4682B4',
    '#008000', '#000080', '#2F4F4F', '#C71585', '#191970', '#32CD32',
    '#483D8B', '#8B008B', '#006400', '#1E90FF'
]
line_styles = ['-', '--', ':', '-.']
style_map = {}
for i in range(1, 17):  # R1 to R16
    col_name = f'R{i}'
    color = colors_16[i-1] if i <= len(colors_16) else 'black'
    line_style = line_styles[(i-1) % len(line_styles)]
    style_map[col_name] = {'color': color, 'linestyle': line_style}

# Calculate diversity parameters
results = []
all_richness = []
hill_results = {}
q_values = np.linspace(0, 3, 50)

for col in columns_to_parse:
    all_side_chains = [item for sublist in parsed_data[col] for item in sublist]
    total_count = len(all_side_chains)
    
    if total_count == 0:
        richness = 0
        shannon_index = 0
        pielou_evenness = 0
        hill_numbers = [0] * len(q_values)
        hill_q1 = 0
        sidechain_dist = ""
        sidechain_props = ""
        sidechain_types = ""
        simpson_index = 0
        hill_q0 = 0
        hill_q2 = 0
    else:
        side_chain_counts = pd.Series(all_side_chains).value_counts()
        richness = len(side_chain_counts)
        all_richness.append(richness)
        
        p_i = side_chain_counts / total_count
        shannon_index = -sum(p_i * np.log(p_i))
        
        pielou_evenness = shannon_index / np.log(richness) if richness > 1 else 0
        
        # Calculate Hill numbers
        hill_numbers = calculate_hill_numbers(side_chain_counts.values, q_values)
        hill_q1 = calculate_hill_numbers(side_chain_counts.values, [1])[0]
        
        # Extra information
        sidechain_dist = str(side_chain_counts.to_dict())
        sidechain_props = str((side_chain_counts / total_count).to_dict())
        sidechain_types = ",".join(side_chain_counts.index)
        simpson_index = sum(p_i ** 2)
        hill_q0 = richness  # q=0 Hill number = S
        hill_q2 = 1 / simpson_index  # q=2 Hill number
    
    results.append({
        'Column': col,
        'Total_Sidechains': total_count,  # Added: total number of side chains
        'Sidechain_Types': sidechain_types,  # Added: side chain type list
        'Sidechain_Distribution': sidechain_dist,  # Added: side chain distribution
        'Sidechain_Proportions': sidechain_props,  # Added: side chain proportions
        'Richness': richness,
        'Shannon_Index': shannon_index,
        'Pielou_Evenness': pielou_evenness,
        'Simpson_Index': simpson_index,  # Added: Simpson index
        'Hill_q0': hill_q0,  # Added: Hill number at q=0
        'Hill_q1': hill_q1,
        'Hill_q2': hill_q2,  # Added: Hill number at q=2
    })
    hill_results[col] = hill_numbers

# Convert to DataFrame
results_df = pd.DataFrame(results)

# Normalize ^1D and Shannon index
D1_min = 1
D1_max = max(results_df['Hill_q1']) if max(results_df['Hill_q1']) > 1 else 2
results_df['Hill_q1_norm'] = (results_df['Hill_q1'] - D1_min) / (D1_max - D1_min)

H_max = np.log(max(all_richness)) if all_richness else 1
results_df['Shannon_Index_norm'] = results_df['Shannon_Index'] / H_max  # Added: normalized Shannon index

# Normalize richness S
S_min = 1
S_max = max(all_richness) if all_richness else 1
results_df['S_norm'] = (results_df['Richness'] - S_min) / (S_max - S_min)

# Define multi-color gradient for S_norm
colors_s = ['#4B0082', 'blue', 'cyan', 'green', 'yellow']
cmap_s = LinearSegmentedColormap.from_list('custom_s', colors_s, N=256)
norm_s = Normalize(vmin=0, vmax=1)
sm_s = ScalarMappable(cmap=cmap_s, norm=norm_s)

# Define five-color gradient for Pielou_Evenness and Hill_q1_norm
colors_j = ['#00CED1', 'gold', 'orange', '#FF4500', 'red']
color_positions_j = [0.0, 0.25, 0.5, 0.75, 1.0]
color_pos_pairs_j = list(zip(color_positions_j, colors_j))
cmap_j = LinearSegmentedColormap.from_list('custom_j', color_pos_pairs_j, N=256)
norm_j = Normalize(vmin=0, vmax=1)
sm_j = ScalarMappable(cmap=cmap_j, norm=norm_j)

# Generate hexadecimal color codes
results_df['Color_S'] = results_df['S_norm'].apply(lambda x: to_hex(sm_s.to_rgba(x)))
results_df['Color_J'] = results_df['Pielou_Evenness'].apply(lambda x: to_hex(sm_j.to_rgba(x)))
results_df['Color_Hill_q1'] = results_df['Hill_q1_norm'].apply(lambda x: to_hex(sm_j.to_rgba(x)))

# Save to Excel file
results_df.to_excel('side_chain_diversity_analysis.xlsx', index=False)

# Plot color scheme figure and save as PDF
fig = plt.figure(figsize=(3, 6))
ax1 = fig.add_axes([0.1, 0.1, 0.2, 0.8])
ax2 = fig.add_axes([0.4, 0.1, 0.2, 0.8])
ax3 = fig.add_axes([0.7, 0.1, 0.2, 0.8])

cbar_s = fig.colorbar(sm_s, cax=ax1, orientation='vertical')
cbar_s.set_label('S_norm (Richness)', size=12)

cbar_j = fig.colorbar(sm_j, cax=ax2, orientation='vertical')
cbar_j.set_label('Pielou_Evenness', size=12)

cbar_hill = fig.colorbar(sm_j, cax=ax3, orientation='vertical')
cbar_hill.set_label('Hill_q1_norm', size=12)

plt.savefig('color_scheme.pdf', bbox_inches='tight')
plt.close()

# Plot diversity profile (by site)
plt.figure(figsize=(12, 8))
for col, hills in hill_results.items():
    style = style_map.get(col, {'color': 'black', 'linestyle': '-'})
    plt.plot(q_values, hills, label=col, color=style['color'], linestyle=style['linestyle'], linewidth=2, alpha=0.8)
plt.xlabel('q value', fontsize=12)
plt.ylabel('Hill Number (^qD)', fontsize=12)
plt.title('Diversity Profile Across Sites', fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('diversity_profile.pdf', bbox_inches='tight')
plt.close()

print("Analysis complete, results saved to 'side_chain_diversity_analysis.xlsx', color scheme figure saved to 'color_scheme.pdf', diversity profile saved to 'diversity_profile.pdf'")
