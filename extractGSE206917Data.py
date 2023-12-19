import pandas as pd
import os

# 设置数据文件夹路径
raw_data_path = 'datasets/raw/GSE206917'
processed_data_path = 'datasets/processed'

# 确保处理过的数据文件夹存在
os.makedirs(processed_data_path, exist_ok=True)

# 读取组织基因表达数据
tissue_df = pd.read_excel(os.path.join(raw_data_path, 'GSE206917_tissue_gene_expression.xls'))
tissue_gene_symbols = tissue_df['gene_symbol']

# 初始化计数器
rr_count = 1
con_count = 1

# 初始化列名映射
column_mapping = {}

for col in tissue_df.columns:
    if 'read_count' in col:
        if 'Recurrent' in col or 'RR' in col:
            new_col = f'RR{rr_count}'
            rr_count += 1
        else:
            new_col = f'Con{con_count}'
            con_count += 1
        column_mapping[col] = new_col

# 应用列名映射
tissue_df.rename(columns=column_mapping, inplace=True)
read_count_cols = list(column_mapping.values())

# 选择gene_symbol和read_count列
tissue_df = tissue_df[['gene_symbol'] + read_count_cols]
final_df = tissue_df.set_index('gene_symbol')

# 定义一个函数来处理和合并数据
def process_and_merge(file_name, final_df):
    global rr_count, con_count
    df = pd.read_excel(os.path.join(raw_data_path, file_name))
    gene_symbols = df['gene_symbol']
    column_mapping = {}

    for col in df.columns:
        if 'read_count' in col:
            if 'Recurrent' in col or 'RR' in col:
                new_col = f'RR{rr_count}'
                rr_count += 1
            else:
                new_col = f'Con{con_count}'
                con_count += 1
            column_mapping[col] = new_col

    df.rename(columns=column_mapping, inplace=True)
    read_count_cols = list(column_mapping.values())
    df = df[['gene_symbol'] + read_count_cols]
    df = df[df['gene_symbol'].isin(final_df.index)]
    df.set_index('gene_symbol', inplace=True)
    return final_df.join(df[read_count_cols], on='gene_symbol', how='left', rsuffix='_' + file_name.split('_')[1])

# 处理其他两个数据集并合并
final_df = process_and_merge('GSE206917_U87_gene_expression.xls', final_df)
final_df = process_and_merge('GSE206917_U251_gene_expression.xls', final_df)

# 在写入CSV之前对列进行排序
# 先提取出所有Con和RR列
con_cols = [col for col in final_df.columns if col.startswith('Con')]
rr_cols = [col for col in final_df.columns if col.startswith('RR')]

# 对DataFrame列进行重新排序，确保gene_symbol作为索引
final_df = final_df[con_cols + rr_cols]

# 将最终结果写入CSV文件
final_df.to_csv(os.path.join(processed_data_path, 'GSE206917.csv'), index=True)
