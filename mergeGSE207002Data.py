import pandas as pd
import os
import re

# 设置文件夹路径
raw_folder = 'datasets/raw/GSE207002'
processed_folder = 'datasets/processed'

# 如果processed文件夹不存在，就创建它
if not os.path.exists(processed_folder):
    os.makedirs(processed_folder)

# 初始化空的DataFrame来存储所有基因和它们的表达值
all_genes_df = pd.DataFrame()

# 获取所有以GSM开头的txt文件
files = sorted([f for f in os.listdir(raw_folder) if f.startswith('GSM') and f.endswith('.txt')])

# 读取每个文件并将数据追加到all_genes_df中
for file in files:
    file_path = os.path.join(raw_folder, file)

    # 使用正则表达式提取RRi_j或CONi_j
    match = re.search(r'(RR|Con)(\d+)_(\d+)', file)
    if match:
        group, i, j = match.groups()
        column_name = f'{group}{3 * (int(i) - 1) + int(j)}'
    else:
        continue

    # 读取文件，跳过表头
    df = pd.read_csv(file_path, sep='\t', header=None, names=['gene_symbol', column_name], skiprows=1)
    # 如果all_genes_df为空，直接赋值
    if all_genes_df.empty:
        all_genes_df = df.set_index('gene_symbol')
    else:
        # 如果不为空，将新数据合并到all_genes_df中
        all_genes_df = all_genes_df.merge(df.set_index('gene_symbol'), left_index=True, right_index=True, how='outer')

# 重置索引
all_genes_df.reset_index(inplace=True)

# 保存到CSV文件，NA值留空
all_genes_df.to_csv(os.path.join(processed_folder, 'GSE207002.csv'), index=False, na_rep='')
