import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# 读取数据
file_path = './datasets/processed/overall.csv'
data = pd.read_csv(file_path, index_col=0)

# 检查NAT10是否在数据中
if 'NAT10' not in data.index:
    raise ValueError("NAT10 gene is not found in the dataset")

# 分离NAT10的表达数据
y = data.loc['NAT10']
X = data.drop('NAT10', axis=0).transpose()

# 标准化数据
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# 应用PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# 获取PCA的第一个主成分的加载量
loadings = pca.components_[0]

# 创建加载量的DataFrame
loading_scores = pd.Series(loadings, index=X.columns)

# 对加载量进行排序
sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)

# 保存到CSV文件
output_path = './result/pca_loadings.csv'
sorted_loading_scores.to_csv(output_path)

print(f"PCA loadings saved to {output_path}")
