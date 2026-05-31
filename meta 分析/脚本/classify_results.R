#! /gpfs/hpc/home/chenchao/hanc/miniconda3/envs/brain_qtl/lib/R/bin/Rscript

setwd("/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/metaqtl/nominal_all/metal_results_hg38")

suppressMessages({
  library(data.table)
  library(ggplot2)
})

# ============================
# 1. 输入输出设置
# ============================
input_file  <- "all_significant_FDR0.05.txt"
output_file <- "all_significant_classified.txt"
fig_file    <- "Direction_Heterogeneity_quadrant.pdf"

message(paste0("[", Sys.time(), "] Reading input file..."))
dt <- fread(input_file, sep = "\t", header = TRUE)

# ============================
# 2. 方向与异质性联合分类
# ============================

# 方向标签（便于画图）
dt[, Direction_Class := fcase(
  Direction %in% c("++", "--"), "Consistent",
  Direction %in% c("+-", "-+"), "Opposite",
  Direction %like% "\\?",       "Ambiguous",
  default = "Other"
)]

# 异质性标签
dt[, Het_Class := fcase(
  HetPVal < 0.05 & HetISq >= 50, "High_Het",
  default = "Low_Het"
)]

# ---- 核心：论文级最终分类 ----
dt[, Final_Class := fcase(
  Direction_Class == "Consistent" & Het_Class == "Low_Het",
    "Consistent_LowHet",
  
  Direction_Class == "Consistent" & Het_Class == "High_Het",
    "Consistent_HighHet",
  
  Direction_Class == "Opposite",
    "Opposite",
  
  Direction_Class == "Ambiguous",
    "Ambiguous",
  
  default = "Other"
)]

# ============================
# 3. 分类统计摘要
# ============================

message("\n========= Final class summary =========")
print(dt[, .N, by = Final_Class][order(-N)])

message("\n========= Direction x Heterogeneity =========")
print(dcast(
  dt,
  Direction_Class ~ Het_Class,
  fun.aggregate = length,
  value.var = "Direction_Class"
))

# ============================
# 4. 保存分类结果
# ============================

message(paste0("[", Sys.time(), "] Writing classified results..."))
fwrite(dt, output_file, sep = "\t", quote = FALSE)

# ============================
# 5. 方向 × 异质性 四象限图
# ============================

# 只画有明确方向的（去掉 Ambiguous）
plot_dt <- dt[Direction_Class %in% c("Consistent", "Opposite")]

# 为了画成“方向轴”，构造一个数值变量
# Consistent = +1, Opposite = -1
plot_dt[, Direction_Num := ifelse(Direction_Class == "Consistent", 1, -1)]

p <- ggplot(
  plot_dt,
  aes(
    x = Direction_Num,
    y = HetISq,
    color = Final_Class
  )
) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6, size = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  scale_x_continuous(
    breaks = c(-1, 1),
    labels = c("Opposite direction", "Consistent direction")
  ) +
  scale_y_continuous(
    name = expression("Heterogeneity ("*I^2*" %)"),
    limits = c(0, 100)
  ) +
  labs(
    x = NULL,
    color = "Meta-QTL class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11)
  )

ggsave(fig_file, p, width = 6.5, height = 5)

message(paste0("[", Sys.time(), "] Done. Figure saved to ", fig_file))
