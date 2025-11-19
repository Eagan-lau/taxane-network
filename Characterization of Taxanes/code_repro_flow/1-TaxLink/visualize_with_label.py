import pandas as pd
import gravis as gv
import networkx as nx

# 加载邻接矩阵
file_path = 'data/linkage_new.csv'
adj_matrix = pd.read_csv(file_path, index_col=0)

# 将邻接矩阵的索引和列转换为整数，确保一致性
adj_matrix.index = adj_matrix.index.astype(int)
adj_matrix.columns = adj_matrix.columns.astype(int)

# 确认邻接矩阵的值分布是否符合预期
print(adj_matrix.describe())

# 创建一个 NetworkX 图
G = nx.Graph()

# 根据邻接矩阵添加边，1 表示有连接，0 表示无连接
for i, row in adj_matrix.iterrows():
    for j, value in row.items():
        # 如果 i 和 j 不同，并且值为 1，则添加边
        if int(i) != int(j) and value == 1.0:
            G.add_edge(i, j)

# 加载 SMILES 数据
smiles_path = 'data/TaxList1.csv'
df_smiles = pd.read_csv(smiles_path)['Isomeric SMILES']

# 确保 SMILES 数据的索引与图节点一致
df_smiles = df_smiles.reset_index(drop=True)

# 打印 SMILES 数据的前几行以确认正确加载
print(df_smiles.head())

# 设置节点的颜色、大小和其他属性
for node in G.nodes:
    G.nodes[node]['color'] = '#007fff'  # 设置节点默认颜色
    G.nodes[node]['size'] = 60  # 设置节点默认大小
    G.nodes[node]['label_size'] = 30  # 设置节点标签大小
    G.nodes[node]['border_size'] = 2  # 设置节点边框宽度
    G.nodes[node]['label'] = node  # 节点标签为节点号

    # 如果节点号超出 SMILES 数据范围，跳过处理
    if int(node) >= len(df_smiles):
        G.nodes[node]['click'] = 'No SMILES data'
    else:
        smiles = str(df_smiles[int(node)])  # 获取 SMILES 数据
        G.nodes[node]['click'] = f'smiles: {smiles}'

# 使用 gravis 进行可视化
fig = gv.d3(G, node_label_data_source='label')

# 确保生成包含 JavaScript 渲染代码的完整 HTML 文件
html_content = fig.to_html()

# 保存 HTML 文件
output_path = 'graph_600_smiles.html'
with open(output_path, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"Graph visualization saved to {output_path}")
