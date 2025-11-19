import matplotlib.pyplot as plt
import pandas as pd

# 假设您的数据在一个名为 'output_position_stastics.csv' 的文件中
data = pd.read_csv('output_position_stastics.csv')

# 设置图形大小
plt.figure(figsize=(10, 8))

# 绘制散点图比较初始值和最终值
#plt.scatter(data.index, data['Initial_RMSD_loss'], alpha=0.5, label='Initial RMSD Loss', color='blue')
#plt.scatter(data.index, data['Final_RMSD_loss'], alpha=0.5, label='Final RMSD Loss', color='red')
plt.scatter(data.index, data['Initial_RMSD_loss'], alpha=0.5, color='blue')
plt.scatter(data.index, data['Final_RMSD_loss'], alpha=0.5, color='red')

# 添加图例
plt.legend()

# 设置标题
#plt.title('Comparison of Initial and Final RMSD Loss Across Proteins')

# 设置 X 轴标签
plt.xlabel('Protein Index', fontsize=20)

# 设置 Y 轴标签
plt.ylabel('RMSD Loss', fontsize=20)

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

# 显示图形
plt.tight_layout()
plt.show()

# 计算 RMSD Loss 的变化量
data['Delta_RMSD_loss'] = data['Final_RMSD_loss'] - data['Initial_RMSD_loss']

# 绘制直方图展示变化量的分布
plt.figure(figsize=(10, 8))
plt.hist(data['Delta_RMSD_loss'], bins=50, color='green', alpha=0.7)
# 设置横坐标范围，例如如果您想要显示 -0.2 到 0.2 的变化
plt.xlim(-0.2, 0.1)
#plt.title('Distribution of Change in RMSD Loss')
plt.xlabel('Change in RMSD Loss', fontsize=20)
plt.ylabel('Frequency', fontsize=20)
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

plt.show()
