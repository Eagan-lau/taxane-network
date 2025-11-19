import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# 1. 读取 cluster_index.txt，并解析 {name} 和 {num}
cluster_df = pd.read_csv('cluster_index.txt', sep='\t', header=None, names=['name', 'num'])

# 2. 读取 eggnog.xlsx 和 clean.xlsx，并指定列名
eggnog_df = pd.read_excel('eggnog.xlsx')
clean_df = pd.read_excel('clean.xlsx')

expected_columns = ['Column1', 'Column2', 'Column3']
if eggnog_df.columns.tolist() != expected_columns:
    eggnog_df.columns = expected_columns
if clean_df.columns.tolist() != expected_columns:
    clean_df.columns = expected_columns

categories = ['1', '2', '3', '4', '5', '6', '7', '-']
eggnog_percentages = {}
clean_percentages = {}
eggnog_counts_dict = {}
clean_counts_dict = {}
protein_nums = {}
group_names = cluster_df['name'].tolist()

for idx, row in cluster_df.iterrows():
    name = row['name']
    num = row['num']
    eggnog_group_data = eggnog_df[eggnog_df['Column2'] == name]
    clean_group_data = clean_df[clean_df['Column2'] == name]
    eggnog_data = eggnog_group_data['Column3'].astype(str).fillna('-')
    clean_data = clean_group_data['Column3'].astype(str).fillna('-')
    eggnog_counts = {cat: 0 for cat in categories}
    clean_counts = {cat: 0 for cat in categories}
    for value in eggnog_data:
        if value in categories:
            eggnog_counts[value] += 1
        else:
            eggnog_counts['-'] += 1
    for value in clean_data:
        if value in categories:
            clean_counts[value] += 1
        else:
            clean_counts['-'] += 1
    total_entries_eggnog = len(eggnog_data)
    total_entries_clean = len(clean_data)
    if total_entries_eggnog > 0:
        eggnog_percentages[name] = {cat: eggnog_counts[cat] * 100 / total_entries_eggnog for cat in categories}
    else:
        eggnog_percentages[name] = {cat: 0 for cat in categories}
    if total_entries_clean > 0:
        clean_percentages[name] = {cat: clean_counts[cat] * 100 / total_entries_clean for cat in categories}
    else:
        clean_percentages[name] = {cat: 0 for cat in categories}
    eggnog_counts_dict[name] = eggnog_counts
    clean_counts_dict[name] = clean_counts
    protein_nums[name] = float(num)

eggnog_percent_df = pd.DataFrame.from_dict(eggnog_percentages, orient='index')
clean_percent_df = pd.DataFrame.from_dict(clean_percentages, orient='index')
eggnog_percent_df = eggnog_percent_df.loc[group_names]
clean_percent_df = clean_percent_df.loc[group_names]
protein_nums_list = [protein_nums[name] for name in group_names]

N = len(group_names)
group_spacing = 1.2
ind = np.arange(N) * group_spacing
width = 0.7
eggnog_stacks = {}
clean_stacks = {}
for cat in categories:
    eggnog_stacks[cat] = eggnog_percent_df[cat].tolist()
    clean_stacks[cat] = clean_percent_df[cat].tolist()

# 灰色部分的"数量"
eggnog_gray_counts = [eggnog_counts_dict[name]['-'] for name in group_names]
clean_gray_counts = [clean_counts_dict[name]['-'] for name in group_names]

# 自动计算百分比最大绝对值用于Y轴自适应
all_stacks = np.array([eggnog_stacks[cat] + [-x for x in clean_stacks[cat]] for cat in categories])
y_limit = np.ceil(np.max(np.abs(all_stacks.sum(axis=0))) * 1.1)
if y_limit < 10: y_limit = 10   # 防止数据极小溢出

# 颜色字典
colors_dict = {
    '1': '#8dd3c7',
    '2': '#ffffb3',
    '3': '#bebada',
    '4': '#fb8072',
    '5': '#80b1d3',
    '6': '#fdb462',
    '7': '#b3de69',
    '-': 'gray'
}

fig, ax1 = plt.subplots(figsize=(18, 9))

# eggNOG柱形图
bottom_eggnog = np.zeros(N)
for cat in categories:
    values = eggnog_stacks[cat]
    ax1.bar(ind, values, width, bottom=bottom_eggnog,
            color=colors_dict[cat], edgecolor='black', linewidth=1.0)
    bottom_eggnog += values

# CLEAN镜像柱形图，注意bottom_clean要逐步更负
bottom_clean = np.zeros(N)
for cat in categories:
    values = clean_stacks[cat]
    neg_values = -np.array(values)   # array便于计算
    ax1.bar(ind, neg_values, width, bottom=bottom_clean,
            color=colors_dict[cat], edgecolor='none', linewidth=0)
    bottom_clean += neg_values   # 注意，这里逐步“更负”

ax1.set_xticks(ind)
ax1.set_xticklabels(group_names, rotation=45, ha='right')
ax1.set_xlabel('Cluster', fontsize=15)
ax1.set_ylabel('Percentage (%)', fontsize=15)
ax1.set_ylim(-y_limit, y_limit)
ax1.axhline(0, color='black', linewidth=1)

# 新增右侧y轴，两条蓝色折线
ax2 = ax1.twinx()
ax2.plot(ind, eggnog_gray_counts, color='blue', marker='o', linestyle='-', linewidth=2, label="eggNOG '-' Count")
ax2.plot(ind, [-x for x in clean_gray_counts], color='blue', marker='s', linestyle='--', linewidth=2, label="CLEAN '-' Count")
ax2.set_ylabel("Count of '-' (gray) category", fontsize=15)
max_count = max(max(eggnog_gray_counts), max(clean_gray_counts))
ax2.set_ylim(-max_count*1.15, max_count*1.15)

# 图例
legend_elements = [
    Patch(facecolor=colors_dict[cat], edgecolor='black', label=cat) for cat in categories
]
first_legend = ax1.legend(handles=legend_elements, title='Categories', bbox_to_anchor=(1.05, 1), loc='upper left')
edge_legend_elements = [
    Patch(facecolor='white', edgecolor='black', label='eggNOG (Top)'),
    Patch(facecolor='white', edgecolor='none', label='CLEAN (Bottom)')
]
second_legend = ax1.legend(handles=edge_legend_elements, title='Data Source', bbox_to_anchor=(1.05, 0.8), loc='upper left')
ax1.add_artist(first_legend)
ax2.legend(title="Gray Category Count", bbox_to_anchor=(1.05, 0.55), loc='upper left')

plt.tight_layout()
plt.show()
