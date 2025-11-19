import matplotlib.pyplot as plt
import pandas as pd

# 载入数据
data = pd.read_csv('output_position_stastics.csv')

# 绘制散点图比较初始值和最终值，使用对数刻度
plt.figure(figsize=(10, 8))
plt.scatter(data.index, data['Initial_Steric_loss'], alpha=0.5, label='Initial Steric Loss', color='blue')
plt.scatter(data.index, data['Final_Steric_loss'], alpha=0.5, label='Final Steric Loss', color='red')
#plt.scatter(data.index, data['Initial_Steric_loss'], alpha=0.5, color='blue')
#plt.scatter(data.index, data['Final_Steric_loss'], alpha=0.5, color='red')
plt.yscale('log')  # 使用对数刻度
plt.legend()
#plt.title('Comparison of Initial and Final Steric Loss Across Proteins')
plt.xlabel('Protein Index', fontsize=20)
plt.ylabel('Steric Loss (log scale)', fontsize=20)

# Adjust spines to remove top and right borders and set linewidth for the remaining spines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2)  # Set the thickness of the left spine
ax.spines['bottom'].set_linewidth(2)  # Set the thickness of the bottom spine

# Adjust tick parameters
ax.tick_params(axis='y', direction='in', length=6, width=2, labelsize=20)  # Y-axis ticks inside and thicker
ax.tick_params(axis='x', direction='out', length=6, width=2, labelsize=20)  # X-axis ticks outside and thicker

# Display the legend
plt.legend(fontsize=20)

# Show the plot
plt.tight_layout()
plt.show()


