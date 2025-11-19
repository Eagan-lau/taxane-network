import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# 读取本地CSV文件
file_path = 'AF2_hap_data.csv'  # 替换为你的CSV文件路径
df = pd.read_csv(file_path)

# 查看数据前几行（用于验证数据是否正确读取）
print(df.head())

# 假设第一列是蛋白ID
protein_ids = df.iloc[:, 0]

# 创建自定义的颜色映射：从非常淡的蓝色到深蓝色
colors = [(0.9, 0.9, 1), (0.7, 0.7, 1), (0.4, 0.4, 1), (0.2, 0.2, 0.8), (0, 0, 0.5)]
nodes = [0, 0.2, 0.4, 0.6, 1]  # 设置颜色渐变的节点

# 创建LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("custom_blues", list(zip(nodes, colors)))

# 创建带有颜色渐变的散点图
plt.figure(figsize=(10, 6))

# 使用自定义颜色映射来显示密度
plt.hexbin(df['pLDDT'], df['pTM'], gridsize=50, cmap=cmap, mincnt=1)

# 调整颜色条
cb = plt.colorbar(label='Counts')

# 调整字体和坐标
plt.title('Hexbin Plot of pLDDT vs pTM Scores', fontsize=14)
plt.xlabel('pLDDT Score', fontsize=12)
plt.ylabel('pTM Score', fontsize=12)

# 设置坐标轴刻度字体大小
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

# 添加垂直和水平参考线
plt.axvline(x=70, color='blue', linestyle='--', linewidth=1.5)  # 垂直线
plt.axhline(y=0.5, color='blue', linestyle='--', linewidth=1.5)  # 水平线

# 显示网格线
plt.grid(False)

# 计算四个区域的数据百分比
total_points = len(df)
q1 = df[(df['pLDDT'] > 70) & (df['pTM'] > 0.5)]
q2 = df[(df['pLDDT'] <= 70) & (df['pTM'] > 0.5)]
q3 = df[(df['pLDDT'] <= 70) & (df['pTM'] <= 0.5)]
q4 = df[(df['pLDDT'] > 70) & (df['pTM'] <= 0.5)]

# 计算百分比
q1_percentage = len(q1) / total_points * 100
q2_percentage = len(q2) / total_points * 100
q3_percentage = len(q3) / total_points * 100
q4_percentage = len(q4) / total_points * 100

# 添加图中图到左上方
ax_inset = inset_axes(plt.gca(), width="30%", height="30%", loc='upper left')
ax_inset.bar(['Q1', 'Q2', 'Q3', 'Q4'], 
             [q1_percentage, q2_percentage, q3_percentage, q4_percentage], color=['red', 'green', 'blue', 'orange'])

# 保存每个区域的蛋白ID到本地文本文件
q1_ids = q1.iloc[:, 0]  # 假设第一列是蛋白ID
q2_ids = q2.iloc[:, 0]
q3_ids = q3.iloc[:, 0]
q4_ids = q4.iloc[:, 0]

q1_ids.to_csv('Q1_protein_ids.txt', index=False, header=False)
q2_ids.to_csv('Q2_protein_ids.txt', index=False, header=False)
q3_ids.to_csv('Q3_protein_ids.txt', index=False, header=False)
q4_ids.to_csv('Q4_protein_ids.txt', index=False, header=False)

# 显示图表
plt.show()
