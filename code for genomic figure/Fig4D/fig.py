import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import FuncFormatter
import warnings

# Ignore warnings about log(0)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Load data
file_path = 'output_position_stastics.csv'  # Ensure this is the correct file path
data = pd.read_csv(file_path)

# Calculate the changes in different dimensions
data['Delta_RMSD_loss'] = data['Final_RMSD_loss'] - data['Initial_RMSD_loss']
data['Delta_Steric_loss'] = data['Final_Steric_loss'] - data['Initial_Steric_loss']
data['Delta_Angle_Dev_loss'] = data['Final_Angle_Dev_loss'] - data['Initial_Angle_Dev_loss']

# Function to format y-tick labels
def format_func(value, tick_number):
    if value == 0:
        return "0"
    else:
        return "{:.0e}".format(value)  # Scientific notation

# Function to get y-ticks based on data
def get_y_ticks(n):
    n = n[n > 0]  # Exclude zero frequencies
    if len(n) == 0:
        return []
    min_freq = n.min()
    max_freq = n.max()
    min_exp = int(np.floor(np.log10(min_freq)))
    max_exp = int(np.ceil(np.log10(max_freq)))
    y_ticks = [10**i for i in range(min_exp, max_exp+1)]
    return y_ticks

# Create the figure and subplots
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 24))  # Adjust the figure size

# Plot configurations for each subplot
plot_configs = [
    {
        'data': data['Delta_RMSD_loss'],
        'range': (-0.2, 0.1),
        'color': 'green',
        'title': 'Change in RMSD Loss',
        'xlabel': 'Change in RMSD Loss'
    },
    {
        'data': data['Delta_Steric_loss'],
        'range': (-100, 0),
        'color': 'skyblue',
        'title': 'Change in Steric Loss',
        'xlabel': 'Change in Steric Loss'
    },
    {
        'data': data['Delta_Angle_Dev_loss'],
        'range': None,  # Use default range
        'color': 'orange',
        'title': 'Change in Angle Deviation Loss',
        'xlabel': 'Change in Angle Deviation Loss'
    }
]

for i, config in enumerate(plot_configs):
    # Plot histogram
    n, bins, patches = axes[i].hist(
        config['data'],
        bins=50,
        range=config['range'],
        color=config['color'],
        alpha=0.7
    )
    axes[i].set_title(config['title'], fontsize=16)
    axes[i].set_xlabel(config['xlabel'], fontsize=14)
    axes[i].set_ylabel('Frequency (Number of Samples)', fontsize=14)
    axes[i].tick_params(axis='x', labelsize=12)
    axes[i].tick_params(axis='y', labelsize=12)
    axes[i].set_yscale('log')  # Set Y-axis to logarithmic scale

    # Remove minor ticks on the left Y-axis
    axes[i].tick_params(axis='y', which='minor', length=0)

    # Set Y-ticks based on data
    y_ticks = get_y_ticks(n)
    axes[i].set_yticks(y_ticks)
    axes[i].yaxis.set_major_formatter(FuncFormatter(format_func))

    # Add right Y-axis showing log10 of frequency
    axes_right = axes[i].twinx()
    axes_right.set_ylabel('Log Value (log₁₀ Frequency)', fontsize=14)
    axes_right.set_ylim(axes[i].get_ylim())

    # Set right Y-axis ticks to match left Y-axis ticks
    axes_right.set_yticks(y_ticks)
    # Set right Y-axis tick labels to log10 of frequency
    axes_right.set_yticklabels([f"{np.log10(y):.1f}" if y > 0 else "" for y in y_ticks])
    axes_right.tick_params(axis='y', labelsize=12)

    # Remove minor ticks on the right Y-axis
    axes_right.tick_params(axis='y', which='minor', length=0)

# Adjust layout
plt.tight_layout()

# Show plot
plt.show()
