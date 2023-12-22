import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 读取预测出的下游基因列表
predicted_genes_path = '../result/top_genes_in_selected_criteria.csv'
predicted_genes_df = pd.read_csv(predicted_genes_path)

# 读取基因表达数据文件
expression_data_path = '../datasets/processed/overall.csv'
expression_data_df = pd.read_csv(expression_data_path, index_col=0)

# 假设NAT10在DataFrame中的名称就是'NAT10'
# 添加NAT10到预测的下游基因列表中
genes_of_interest = predicted_genes_df['Gene'].tolist() + ['NAT10']

# 筛选出NAT10及其潜在的10个下游基因在不同条件下的表达水平
selected_expression_data = expression_data_df.loc[genes_of_interest]

# 绘制热图
plt.figure(figsize=(12, 9))  # 调整图形的大小以适合您的数据
sns.heatmap(
    selected_expression_data,
    cmap='viridis',  # 更改颜色映射
    annot=False,     # 是否在热图上显示数值
    cbar=True        # 显示颜色条
)
plt.title('Heatmap of NAT10 and Predicted Downstream Genes')
plt.xlabel('Conditions')
plt.ylabel('Genes')

# 定义热图的保存路径
heatmap_output_path = '../result/heatmap_with_NAT10.png'

# 保存热图
plt.savefig(heatmap_output_path, dpi=300)

plt.show()
