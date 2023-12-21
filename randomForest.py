import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split

# 读取数据
file_path = './datasets/processed/GSE206917.csv'
data = pd.read_csv(file_path, index_col=0)

# 检查NAT10是否在数据中
if 'NAT10' not in data.index:
    raise ValueError("NAT10 gene is not found in the dataset")

# 分离NAT10的表达数据作为目标变量
y = data.loc['NAT10']
X = data.drop('NAT10', axis=0).transpose()

# 分割数据集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 构建随机森林模型
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train, y_train)

# 获取特征重要性
feature_importances = pd.Series(model.feature_importances_, index=X.columns)

# 对特征重要性进行排序
sorted_importances = feature_importances.sort_values(ascending=False)

# 保存到CSV文件
output_path = './result/feature_importances.csv'
sorted_importances.to_csv(output_path)

print(f"Feature importances saved to {output_path}")
