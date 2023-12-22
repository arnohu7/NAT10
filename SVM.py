import pandas as pd
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

# 读取数据
file_path = './datasets/processed/overall.csv'
data = pd.read_csv(file_path, index_col=0)

# 检查NAT10是否在数据中
if 'NAT10' not in data.index:
    raise ValueError("NAT10 gene is not found in the dataset")

# 分离NAT10的表达数据作为目标变量
y = data.loc['NAT10']
X = data.drop('NAT10', axis=0).transpose()

# 分割数据集
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 使用标准化和SVM构建管道
pipeline = make_pipeline(StandardScaler(), SVR(kernel='linear'))

# 训练模型
pipeline.fit(X_train, y_train)

# 获取特征重要性（在线性SVM中，系数可以反映特征重要性）
coefficients = pipeline.named_steps['svr'].coef_

# 创建特征重要性的Series
feature_importances = pd.Series(coefficients.flatten(), index=X.columns)

# 对特征重要性进行排序
sorted_importances = feature_importances.abs().sort_values(ascending=False)

# 保存到CSV文件
output_path = './result/svm_feature_importances.csv'
sorted_importances.to_csv(output_path)

print(f"SVM feature importances saved to {output_path}")
