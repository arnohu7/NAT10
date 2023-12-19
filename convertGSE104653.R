# 加载必要的库
library(data.table) # 用于读取数据文件
library(hugene10sttranscriptcluster.db) # 用于 GPL6244 平台文件

# 读取数据文件
exprData <- read.table("./datasets/raw/GSE82139/GSE82139_series_matrix.txt", 
                  header = TRUE, sep = "\t", quote = "", fill = TRUE,
                  comment.char = "!")

# 将数据中的行号改为 ID 号，更改完毕后将 ID 号删除
exprData_probe_ids <- exprData[, 1] # 获取 ID 号
rownames(exprData) <- exprData_probe_ids # 将原来的行号设置为 ID 号
exprData <- exprData[, -1] # 将原来的 ID 号删除

# 提取 probe_id 和 symbol 的对应关系
idSymbolMap <- toTable(hugene10sttranscriptclusterSYMBOL)

# 将没有基因名的探针筛掉
idSymbolMap <- idSymbolMap[idSymbolMap$symbol != "", ]
exprData <- exprData[rownames(exprData) %in% idSymbolMap$probe_id,]

# 由于 exprData 每一个 symbol 可能对应多个探针，因此要进行处理，使得每一个 symbol 对应一个探针
# 使用 by 进行分组，同一个 symbol 的探针分到同一组，接着求解每组中每一行的均值，将均值最大的那一行返回
groupByProbeId = by(exprData, idSymbolMap$symbol, function(x) rownames(x)[which.max(rowMeans(x))])
groupByProbeId = as.character(groupByProbeId) # 将分组结果转换为纯字符型向量
exprData = exprData[rownames(exprData) %in% groupByProbeId, ] # 过滤掉基因多余的探针

# 将表达矩阵的行名转换为基因的 symbol
rownames(exprData) = idSymbolMap[match(rownames(exprData), idSymbolMap$probe_id), 2]

# 将数据写入到 csv 文件中
write.csv(exprData, file = "./datasets/raw/GSE82139/GSE82139.csv", row.names = TRUE)
