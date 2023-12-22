# 加载所需的库
library(limma)
library(dplyr)
library(readr)

# 变量定义
dataFile <- "./datasets/processed/overall.csv"
limmaOutFile <- "./result/all.limmaOut.csv"
theNumberOfCon <- 23
theNumberOfRR <- 23
PValueThreshold <- 0.05
significantGenesFile <- "./result/significant_genes_from_limma.csv"
correlationThreshold <- 0.5
highCorrelationGenesFile <- "./result/NAT10_high_correlation_genes.csv"
selectedGenesCommon2BothCriteria <- "./result/selected_genes_common_to_both_criteria.csv"

##################################################
# 进行差异表达分析，并将结果输出为 CSV 文件
##################################################

# 读取数据，header = T 表示 CSV 文件的第一行包含列名，rows.names = 1 表示第一列用作行名
data <- read.csv(dataFile, header = T, row.names = 1)

# 样本信息注释
# 首先生成一个 list 向量，其中包含若干个重复的 Con 和若干个重复的 RR
# 然后将这个向量转换为一个因子（分类变量），其中"Con"和"RR"是两个水平，且这个因子是无序的。
list <- c(rep("Con", theNumberOfCon), rep("RR", theNumberOfRR)) %>% factor(., levels = c("Con", "RR"), ordered = F)
# 接着创建一个设计矩阵，用来进行统计分析，+0 表示不含截距项
list <- model.matrix(~list+0)
colnames(list) <- c("Con", "RR") # 对列名进行改名

# 数据与list进行线性拟合
data.fit <- lmFit(data, list)

# 差异分析
# 首先创建一个对比矩阵，计算 Con 和 RR 的差异
data.matrix <- makeContrasts(Con - RR, levels = list)
# 将对比矩阵应用于 lmFit 的结果，计算指定的差异，简单来说，在 data.fit 基础上进一步精确计算差异
fit <- contrasts.fit(data.fit, data.matrix)
# 这行代码应用了经验贝叶斯（Empirical Bayes）方法来改进 limma 包中的差异表达分析
fit <- eBayes(fit)
# 提取和排序差异表达分析的结果
output <- topTable(fit, n = Inf, adjust = "fdr")

# 导出所有的差异结果
nrDEG = na.omit(output) # 去掉数据中有NA的行或列
write.csv(nrDEG, limmaOutFile)

##################################################
# 选取前几的差异表达基因，输出为另一个文件
##################################################

# 读取 limma 输出的差异表达结果
limma_results <- read.csv(limmaOutFile, header = TRUE, row.names = 1)

# 筛选出调整后的P值小于0.05的基因
significant_genes <- limma_results[limma_results$adj.P.Val < PValueThreshold, ]

# 保存筛选结果到CSV文件
write.csv(significant_genes, significantGenesFile, row.names = TRUE)

##################################################
# 计算 NAT10 与其他基因的相关系数，并输出为文件
##################################################

# 读取 limma 输出的差异表达结果
limma_results <- read.csv(limmaOutFile, header = TRUE, row.names = 1)

# 读取你的基因表达数据
# 假设这个文件是一个表格，行表示样本，列表示基因
expression_data <- read.csv(dataFile, header = TRUE, row.names = 1)

# 提取NAT10的表达数据，确保它是一个向量
NAT10_expression <- expression_data[rownames(expression_data) == "NAT10", ]
NAT10_expression_vector <- as.numeric(NAT10_expression[1, ])

# 计算NAT10与其他所有基因的相关性
correlation_with_NAT10 <- apply(expression_data, 1, function(gene_expression) {
  cor(NAT10_expression_vector, as.numeric(gene_expression), use = "complete.obs")
})

# 将结果转换为数据框
correlation_data <- data.frame(Gene = rownames(expression_data), Correlation = correlation_with_NAT10)

# 筛选正相关性较高且不包含NA的基因
high_correlation_genes <- correlation_data[correlation_data$Correlation > correlationThreshold & !is.na(correlation_data$Correlation), ]

# 根据相关性进行排序，相关性大的在前面
high_correlation_genes_sorted <- high_correlation_genes[order(-high_correlation_genes$Correlation),]

# 保存高相关性基因到CSV文件，忽略包含NA的行
write.csv(high_correlation_genes_sorted, highCorrelationGenesFile, row.names = FALSE)

##################################################
# 找出相关性大，且 P 值低的基因，并输出为文件
##################################################

# 读取显著差异表达的基因
significant_genes <- read.csv(significantGenesFile, header = TRUE, row.names = 1)

# 读取与NAT10高度相关的基因
high_correlation_genes <- read.csv(highCorrelationGenesFile, header = TRUE)

# 将显著性数据和相关性数据合并
combined_genes <- merge(significant_genes, high_correlation_genes, by.x = "row.names", by.y = "Gene")
colnames(combined_genes)[1] <- "Gene"

# 筛选出同时满足两个条件的基因
# 即在significant_genes中且在high_correlation_genes中
common_genes <- intersect(row.names(significant_genes), combined_genes$Gene)

# 创建一个新的数据框，包含这些基因的信息
selected_genes <- combined_genes[combined_genes$Gene %in% common_genes, ]

# 保存筛选结果到CSV文件
write.csv(selected_genes, selectedGenesCommon2BothCriteria, row.names = FALSE)

