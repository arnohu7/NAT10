# 加载所需的库
library(limma)
library(dplyr)

# 读取数据，header = T 表示 CSV 文件的第一行包含列名，rows.names = 1 表示第一列用作行名
data <- read.csv("./datasets/processed/GSE206917.csv", header = T, row.names = 1)

# 样本信息注释
# 首先生成一个 list 向量，其中包含若干个重复的 Con 和若干个重复的 RR
# 然后将这个向量转换为一个因子（分类变量），其中"Con"和"RR"是两个水平，且这个因子是无序的。
list <- c(rep("Con", 8), rep("RR", 8)) %>% factor(., levels = c("Con", "RR"), ordered = F)
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
write.csv(nrDEG, "./result/all.limmaOut.csv")
