# 加载必要的库
library(readr)

# 读取 limma 输出的差异表达结果
limma_results <- read.csv("./result/all.limmaOut.csv", header = TRUE, row.names = 1)

# 读取你的基因表达数据
# 假设这个文件是一个表格，行表示样本，列表示基因
expression_data <- read.csv("./datasets/processed/GSE206917.csv", header = TRUE, row.names = 1)

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
high_correlation_genes <- correlation_data[correlation_data$Correlation > 0.5 & !is.na(correlation_data$Correlation), ]

# 根据相关性进行排序，相关性大的在前面
high_correlation_genes_sorted <- high_correlation_genes[order(-high_correlation_genes$Correlation),]

# 保存高相关性基因到CSV文件，忽略包含NA的行
write.csv(high_correlation_genes_sorted, "./result/NAT10_high_correlation_genes.csv", row.names = FALSE)

# 打印一条消息确认保存完成
print("High correlation genes with NAT10 have been sorted by correlation and saved to NAT10_high_correlation_genes.csv (NA values excluded)")
