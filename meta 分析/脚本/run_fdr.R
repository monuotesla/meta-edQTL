#! /gpfs/hpc/home/chenchao/hanc/miniconda3/envs/brain_qtl/lib/R/bin/Rscript


setwd("/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/metaqtl/nominal_all/metal_results_hg38")
# 1. 设置输入和输出文件名
input_file <- "all_merged_results.txt"
output_file <- "all_significant_FDR0.05.txt"

cat("正在读取文件:", input_file, "...\n")

# 2. 读取数据 (check.names=FALSE 保证列名 P-value 不会被自动改成 P.value)
data <- read.table(input_file, header=TRUE, sep="\t", 
                   stringsAsFactors=FALSE, check.names=FALSE)

# 3. 提取染色体信息
# 假设 MarkerName 格式为 "10_100325140"，我们截取第一个 "_" 之前的部分作为染色体编号
# 如果您的 MarkerName 格式不同，也可以解析 "edSite_File" 列
data$Temp_Chr <- sub("_.*", "", data$MarkerName)

cat("检测到的染色体包括:", unique(data$Temp_Chr), "\n")
cat("正在分染色体进行 FDR 矫正...\n")

# 4. 分染色体进行 FDR (Benjamini-Hochberg) 矫正
# split() 将数据按染色体拆分，lapply() 对每一组单独做 p.adjust
data_list <- split(data, data$Temp_Chr)

data_list_adjusted <- lapply(data_list, function(df) {
  # 确保 P-value 是数值型
  df$`P-value` <- as.numeric(df$`P-value`)
  
  # 计算 q-value (FDR)
  df$q_value <- p.adjust(df$`P-value`, method = "BH")
  
  return(df)
})

# 5. 合并结果
final_data <- do.call(rbind, data_list_adjusted)

# 6. 筛选 q-value < 0.05
sig_data <- subset(final_data, q_value < 0.05)

# 移除临时生成的 Temp_Chr 列（可选）
sig_data$Temp_Chr <- NULL

# 7. 写入结果
cat("筛选前行数:", nrow(data), "\n")
cat("筛选后行数 (q<0.05):", nrow(sig_data), "\n")
cat("正在写入结果到:", output_file, "...\n")

write.table(sig_data, output_file, sep="\t", 
            row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("完成！\n")