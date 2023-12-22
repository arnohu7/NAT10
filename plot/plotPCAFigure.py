import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# 读取数据
predicted_genes_path = '../result/top_genes_in_selected_criteria.csv'
predicted_genes_df = pd.read_csv(predicted_genes_path)

expression_data_path = '../datasets/processed/overall.csv'
expression_data_df = pd.read_csv(expression_data_path, index_col=0)

# 确保NAT10也包含在数据中
genes_of_interest = predicted_genes_df['Gene'].tolist() + ['NAT10']

# 筛选出相关基因的表达数据
selected_expression_data = expression_data_df.loc[genes_of_interest]

# 标准化数据
scaler = StandardScaler()
scaled_expression_data = scaler.fit_transform(selected_expression_data.T)

# 执行PCA
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_expression_data)

# 将PCA结果转换为DataFrame
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['Condition'] = expression_data_df.columns

# 绘制PCA图
plt.figure(figsize=(10, 8))
sns.scatterplot(x='PC1', y='PC2', hue='Condition', data=pca_df, palette='viridis')
plt.title('PCA of NAT10 and Predicted Downstream Genes')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')

# 保存PCA图
pca_output_path = '../result/pca_with_NAT10.png'
plt.savefig(pca_output_path)

# 如果你想在屏幕上查看这个图
plt.show()
