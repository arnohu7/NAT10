import pandas as pd

# 加载CSV文件
df1 = pd.read_csv('./datasets/processed/GSE82139.csv')
df2 = pd.read_csv('./datasets/processed/GSE206917.csv')
df3 = pd.read_csv('./datasets/processed/GSE207002.csv')

# 修改GSE206917和GSE207002文件的列名
df2.columns = ['gene_symbol'] + [f'Con{int(col[3:]) + 6}' if 'Con' in col else f'RR{int(col[2:]) + 6}' for col in df2.columns[1:]]
df3.columns = ['gene_symbol'] + [f'Con{int(col[3:]) + 14}' if 'Con' in col else f'RR{int(col[2:]) + 14}' for col in df3.columns[1:]]

# 合并文件并只保留在所有文件中都出现、且没有NA值的基因
overall_df = pd.merge(df1, df2, on='gene_symbol', how='inner')
overall_df = pd.merge(overall_df, df3, on='gene_symbol', how='inner')
overall_df = overall_df.dropna()

# 重新排序列，使Con列在前，RR列在后
con_cols = [col for col in overall_df.columns if 'Con' in col]
rr_cols = [col for col in overall_df.columns if 'RR' in col]
overall_df = overall_df[['gene_symbol'] + con_cols + rr_cols]

# 保存到新的CSV文件
overall_df.to_csv('./datasets/processed/overall.csv', index=False)
