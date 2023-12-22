import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# 读取基因表达数据
expression_data_path = '../datasets/processed/overall.csv'
expression_data_df = pd.read_csv(expression_data_path, index_col=0)

# 读取潜在的10个下游基因和它们的自定义得分
top_genes_path = '../result/top_genes_in_selected_criteria.csv'
top_genes_df = pd.read_csv(top_genes_path)
top_genes = top_genes_df['Gene'].tolist()  # 基因列表
scores = top_genes_df['0'].tolist()    # 对应的得分列表

# 添加NAT10到基因列表和得分列表中
# 假设NAT10的得分为最大值或某个特定值
nat10_score = max(scores)  # 或者是预先定义的得分
top_genes.append('NAT10')
scores.append(nat10_score)

# 创建网络图
G = nx.Graph()

# 添加节点和边
for i, gene1 in enumerate(top_genes):
    for j, gene2 in enumerate(top_genes):
        if i < j:  # 防止添加重复的边
            # 使用基因的得分作为边的权重
            # 这里我们简化为使用得分的平均值
            weight = (scores[i] + scores[j]) / 2
            G.add_edge(gene1, gene2, weight=weight)

# 为了使边的权重在图形中可见，我们将它们正规化并乘以一个因子
weights = [G[u][v]['weight'] for u, v in G.edges()]
max_weight = max(weights)
edge_weights = [(w / max_weight) * 10 for w in weights]  # 根据需要调整因子

pos = nx.spring_layout(G)  # 节点的布局
nx.draw(G, pos, with_labels=True, node_size=700,
        node_color='skyblue', edge_color=edge_weights,
        width=edge_weights, edge_cmap=plt.cm.Blues, font_size=5)  # 将字体大小设置为8

# 如果需要保存图片
plt.savefig('../result/gene_interaction_network_based_on_scores.png', dpi=900)


# 显示图
plt.title('Gene Interaction Network Based on Custom Scores')
plt.show()

