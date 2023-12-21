import pandas as pd

# 读取数据
combined_scores = pd.read_csv('result/combined_gene_scores.csv', index_col=0)
selected_genes = pd.read_csv('result/selected_genes_common_to_both_criteria.csv')

# 保证基因名列是字符串
selected_genes['Gene'] = selected_genes['Gene'].astype(str)

# 在combined_gene_scores.csv中查找selected_genes的评分
selected_genes_scores = combined_scores.loc[combined_scores.index.intersection(selected_genes['Gene'])]

# 按评分排序并选择前10个
top_genes = selected_genes_scores.sort_values(by=selected_genes_scores.columns[0], ascending=False).head(10)

# 保存结果
output_path = 'result/top_genes_in_selected_criteria.csv'
top_genes.to_csv(output_path)

print(f"Top genes saved to {output_path}")
