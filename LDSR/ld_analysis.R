#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(rstatix)

# =========================
# 0️⃣ 输入文件
# =========================
input_file <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/class/metaqtl_significant_final_reclassified.txt"
out_dir <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/毕业论文/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dt <- fread(input_file)

# 更新列名
setnames(dt, (ncol(dt)-2):ncol(dt), c("gene", "genomic_region", "alu"))

# =========================
# 1️⃣ 提取 chr
# =========================
dt[, chr := sub("^([^:_]+)[_:].*$", "\\1", MarkerName)]
dt[, chr := gsub("^chr", "", chr)]
dt[, chr := as.integer(chr)]

# =========================
# 2️⃣ 筛选分类
# =========================
target_classes <- c(
  "Consistent_LowHet",
  "Consistent_HighHet",
  "Opposite_HighHet"
)

dt_sub <- dt[Final_Class %in% target_classes]
dt_sub <- unique(dt_sub, by = c("Final_Class", "rsid"))

# =========================
# 3️⃣ LD 路径函数
# =========================
ld_path <- function(pop, chr) {
  if (pop == "EAS") {
    sprintf("/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/class/LD/1000G_Phase3_EAS_ldsc.allsnps/1000G_Phase3_EAS_allsnps.%d.l2.ldscore.gz", chr)
  } else {
    sprintf("/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/class/LD/1000G_EUR_Phase3_ldsc.allsnps/1000G_EUR_Phase3_allsnps.%d.l2.ldscore.gz", chr)
  }
}

# =========================
# 4️⃣ 提取 LD
# =========================
ld_results <- list()

for (c in target_classes) {

  cat("Processing class:", c, "\n")
  dt_c <- dt_sub[Final_Class == c]

  for (cc in sort(unique(dt_c$chr))) {

    snps_chr <- dt_c[chr == cc, rsid]
    if (length(snps_chr) == 0) next

    eas_file <- ld_path("EAS", cc)
    eur_file <- ld_path("EUR", cc)

    if (!file.exists(eas_file) || !file.exists(eur_file)) next

    ld_eas <- fread(cmd = paste("zcat", eas_file))
    ld_eur <- fread(cmd = paste("zcat", eur_file))

    setnames(ld_eas, c("CHR", "SNP", "BP", "L2"), paste0("EAS_", c("CHR","SNP","BP","L2")))
    setnames(ld_eur, c("CHR", "SNP", "BP", "L2"), paste0("EUR_", c("CHR","SNP","BP","L2")))

    # ⚠️ 仍然是 inner join（与你原逻辑一致）
    ld_m <- merge(
      ld_eas[EAS_SNP %in% snps_chr],
      ld_eur[EUR_SNP %in% snps_chr],
      by.x = "EAS_SNP",
      by.y = "EUR_SNP"
    )

    if (nrow(ld_m) == 0) next

    ld_m[, Final_Class := c]
    ld_m[, chr := cc]

    ld_results[[paste(c, cc, sep = "_")]] <- ld_m
  }
}

ld_dt <- rbindlist(ld_results, use.names = TRUE)

# =========================
# 5️⃣ 计算 ΔLD
# =========================
ld_dt[, delta_L2 := abs(EAS_L2 - EUR_L2)]

# 合并 HetChiSq
ld_dt <- merge(
  ld_dt,
  unique(dt_sub[, .(rsid, HetChiSq)], by = "rsid"),
  by.x = "EAS_SNP",
  by.y = "rsid",
  all.x = TRUE
)

# =========================
# 6️⃣ 定义两组（核心）
# =========================
ld_dt[, Het_Group := ifelse(
  Final_Class == "Consistent_LowHet",
  "LowHet",
  "HighHet"
)]

# =========================
# 7️⃣ 保存中间总表（强烈建议）
# =========================
fwrite(ld_dt,
       file.path(out_dir, "LD_master_table.txt"),
       sep = "\t")

# =========================
# 8️⃣ 分组数据
# =========================
lowhet_dt  <- ld_dt[Het_Group == "LowHet"]
highhet_dt <- ld_dt[Het_Group == "HighHet"]

# =========================
# 9️⃣ 保存原始（未去重）
# =========================
fwrite(lowhet_dt,
       file.path(out_dir, "LowHet_deltaL2_raw.txt"),
       sep = "\t")

fwrite(highhet_dt,
       file.path(out_dir, "HighHet_deltaL2_raw.txt"),
       sep = "\t")

# =========================
# 🔟 去重（按 SNP）
# =========================
lowhet_unique  <- unique(lowhet_dt,  by = "EAS_SNP")
highhet_unique <- unique(highhet_dt, by = "EAS_SNP")

# =========================
# 1️⃣1️⃣ 保存去重版本
# =========================
fwrite(lowhet_unique,
       file.path(out_dir, "LowHet_deltaL2_unique.txt"),
       sep = "\t")

fwrite(highhet_unique,
       file.path(out_dir, "HighHet_deltaL2_unique.txt"),
       sep = "\t")

# =========================
# 1️⃣2️⃣ 统计检验（用去重后的）
# =========================
wilcox_res <- wilcox_test(
  as.data.frame(rbind(lowhet_unique, highhet_unique)),
  delta_L2 ~ Het_Group,
  p.adjust.method = "fdr"
)

# 保存统计结果
fwrite(as.data.frame(wilcox_res),
       file.path(out_dir, "LowHet_vs_HighHet_Wilcox.txt"),
       sep = "\t")

print(wilcox_res)