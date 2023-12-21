import pandas as pd

# 读取CSV文件
df = pd.read_csv('../datasets/raw/GSE82139/GSE82139.csv')

# 删除指定的列
columns_to_remove = ['X.GSM2184200.', 'X.GSM2184201.']
df.drop(columns_to_remove, axis=1, inplace=True)

# 重命名其他列
rename_mapping = {
    'X.GSM2184194.': 'Con1', 'X.GSM2184195.': 'Con2', 'X.GSM2184196.': 'Con3',
    'X.GSM2184197.': 'Con4', 'X.GSM2184198.': 'Con5', 'X.GSM2184199.': 'Con6',
    'X.GSM2184202.': 'RR1', 'X.GSM2184203.': 'RR2', 'X.GSM2184204.': 'RR3',
    'X.GSM2184205.': 'RR4', 'X.GSM2184206.': 'RR5', 'X.GSM2184207.': 'RR6'
}
df.rename(columns=rename_mapping, inplace=True)

# 将第一列重命名为gene_symbol
df.rename(columns={df.columns[0]: "gene_symbol"}, inplace=True)

# 将结果保存到新的CSV文件
df.to_csv('../datasets/processed/GSE82139.csv', index=False)
