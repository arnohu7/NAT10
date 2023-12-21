import pandas as pd
from sklearn.preprocessing import StandardScaler

# 假设您已经有了这些模型的输出结果，这里我们用三个假设的DataFrame代替
# 请将以下DataFrame替换为您的实际数据
pca_results = pd.read_csv('result/pca_loadings.csv', index_col=0)
svm_results = pd.read_csv('result/svm_feature_importances.csv', index_col=0)
rf_results = pd.read_csv('result/feature_importances.csv', index_col=0)

# 结果标准化
scaler = StandardScaler()
pca_scaled = scaler.fit_transform(pca_results)
svm_scaled = scaler.fit_transform(svm_results)
rf_scaled = scaler.fit_transform(rf_results)

# 将标准化结果转换为DataFrame
pca_df = pd.DataFrame(pca_scaled, index=pca_results.index)
svm_df = pd.DataFrame(svm_scaled, index=svm_results.index)
rf_df = pd.DataFrame(rf_scaled, index=rf_results.index)

# 计算综合评分，这里我们简单地取平均值
combined_score = (pca_df + svm_df + rf_df) / 3

# 将综合评分排序
combined_score_sorted = combined_score.sort_values(by=0, ascending=False)

# 保存综合评分
combined_score_sorted.to_csv('result/combined_gene_scores.csv')

print("Combined gene scores saved to result/combined_gene_scores.csv")
