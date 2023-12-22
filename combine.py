import pandas as pd

# 读取数据
combined_scores = pd.read_csv('result/combined_gene_scores.csv', index_col=0)
selected_genes = pd.read_csv('result/selected_genes_common_to_both_criteria.csv', index_col=0)

# 保证基因名列是字符串
selected_genes.index = selected_genes.index.map(str)

# 在selected_genes中找到需要的列，包括logFC, AveExpr, t, P.Value, adj.P.Val, B
# 我们假设combined_scores里面有一个与基因相关的平均得分列，比如叫做'AverageScore'
# 我们将此列和selected_genes中的列合并，以获得完整的数据信息
combined_data = selected_genes.join(combined_scores)

# 按照'AverageScore'排序并选择前10个基因
top_genes = combined_data.sort_values(by='0', ascending=False).head(10)

# 保存结果
output_path = 'result/top_genes_in_selected_criteria.csv'
top_genes.to_csv(output_path)

