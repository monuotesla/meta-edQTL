library(qvalue)
library(ggplot2)
library(coloc)
library(dplyr)
library(gridExtra)
library(data.table)

# =======================================================
# 2. Locus 可视化 (参数统一与数据整合)
# =======================================================

# ... (图片保存目录和 top_genes 筛选逻辑保持不变) ...
plot_dir="/gpfs/hpc/home/chenchao/mayinuo/edqtl/coloc/withsQTL/result/figure/"

top_genes=sig_genes_list
if (length(top_genes) > 0) {
    
    for (gene in top_genes) {
        
        # --- 关键整合：数据提取和对齐 (与上方逻辑一致) ---
        # 1. 确定 edSite 中心和区域
        edsite_info <- eqtl_full[eqtl_full$eGene == gene, ][1, ]
        if (is.na(edsite_info$site_chr)) next
        chr_k <- edsite_info$site_chr
        pos_k <- edsite_info$site_pos
        region_start <- pos_k - window
        region_end <- pos_k + window

        # 2. 提取 eQTL/edQTL 数据
        eqtl_reg <- eqtl_full[
            eqtl_full$eGene == gene & eqtl_full$snp_chr == chr_k & eqtl_full$snp_pos >= region_start & eqtl_full$snp_pos <= region_end,
        ]
        
        # 3. 提取 GWAS 数据
        gwas_reg <- gwas[
            gwas$snp_chr == chr_k & gwas$snp_pos >= region_start & gwas$snp_pos <= region_end,
        ]

        # 4. SNP 交集与对齐
        snps_common <- intersect(gwas_reg$eSNP, eqtl_reg$eSNP)
        if (length(snps_common) < 2) next

        gwas_use <- gwas_reg[match(snps_common, gwas_reg$eSNP), ]
        eqtl_use <- eqtl_reg[match(snps_common, eqtl_reg$eSNP), ]
        
        # 5. 重新运行 coloc
        res <- coloc.abf(
            dataset1 = list(pvalues = gwas_use$pval, type = "quant", N = 546, snp = snps_common),
            dataset2 = list(pvalues = eqtl_use$p.value, type = "quant", N = 546, snp = snps_common),
            MAF = gwas_use$maf
        )
            
        snp_info <- res$results # 包含 snp 和 SNP.PP.H4 列
        
        # --- 6. 整理绘图数据 (使用交集结果) ---
        plot_df <- data.frame(
            eSNP = snps_common,
            BP = gwas_use$snp_pos, # 使用 GWAS 的位置信息作为 X 轴
            GWAS_P = gwas_use$pval,
            EQTL_P = eqtl_use$p.value
        )
        
        # 将 SNP.PP.H4 合并回绘图数据, by.y = "snp" 需要确定
        plot_df <- merge(plot_df, snp_info[, c("snp", "SNP.PP.H4")], by.x = "eSNP", by.y = "snp")
        
        # 转换为 Mb 单位
        plot_df$BP_Mb <- plot_df$BP / 1e6
        plot_df$SNP.PP.H4[is.na(plot_df$SNP.PP.H4)] <- 0
        
        # --- 2.4 使用 ggplot2 绘图 ---
        
        # 定义子图 1: GWAS
# --- 1. 识别该区域的 Top SNP ---
# 找到 SNP.PP.H4 最大的那一行
top_snp_row <- plot_df[which.max(plot_df$SNP.PP.H4), ]
top_snp_id <- top_snp_row$eSNP

# 重新对 plot_df 排序，确保 PP.H4 大的点最后画（画在最上层）
plot_df <- plot_df %>% arrange(SNP.PP.H4)

# --- 2. 定义子图 1: GWAS ---
p1 <- ggplot(plot_df, aes(x = BP_Mb, y = -log10(GWAS_P))) +
    # 画背景点
    geom_point(aes(color = SNP.PP.H4, size = SNP.PP.H4), alpha = 0.7) +
    # 专门给 Top SNP 画一个大菱形
    geom_point(data = top_snp_row, aes(x = BP_Mb, y = -log10(GWAS_P)), 
               color = "red", shape = 18, size = 6) + 
    # 标注 Top SNP 的 ID
    geom_text(data = top_snp_row, aes(label = eSNP), 
              vjust = -1.5, color = "black", fontface = "bold", size = 3) +
    scale_color_gradient(low = "grey70", high = "red", limits = c(0, 1), name = "P(Causal)") +
    scale_size_continuous(range = c(1.5, 4), guide = "none") +
    expand_limits(y = max(-log10(plot_df$GWAS_P)) * 1.2) + # 留出文字空间
    labs(title = paste0(gene, " - sQTL Signal"), y = "-log10(eQTL P)", x = NULL) +
    theme_bw() + theme(legend.position = "right", axis.text.x = element_blank())

# --- 3. 定义子图 2: edQTL ---
p2 <- ggplot(plot_df, aes(x = BP_Mb, y = -log10(EQTL_P))) +
    geom_point(aes(color = SNP.PP.H4, size = SNP.PP.H4), alpha = 0.7) +
    # 同样给 edQTL 图也标出 Top SNP
    geom_point(data = top_snp_row, aes(x = BP_Mb, y = -log10(EQTL_P)), 
               color = "red", shape = 18, size = 6) +
    geom_text(data = top_snp_row, aes(label = eSNP), 
              vjust = -1.5, color = "black", fontface = "bold", size = 3) +
    scale_color_gradient(low = "grey70", high = "red", limits = c(0, 1), name = "P(Causal)") +
    scale_size_continuous(range = c(1.5, 4), guide = "none") +
    expand_limits(y = max(-log10(plot_df$EQTL_P)) * 1.2) +
    labs(title = paste0(gene, " - edQTL Signal"), 
         y = "-log10(edQTL P)", x = paste("Position on Chr (Mb)")) +
    theme_bw() + theme(legend.position = "right")


            
        # --- 2.5 拼图并保存 ---
       safe_gene_name <- gsub("[:/]", "_", gene)
        out_file <- paste0(plot_dir, safe_gene_name, "_LocusPlot.pdf")
        
        final_plot <- grid.arrange(p1, p2, ncol = 1, heights = c(1, 1))
        
        # 使用 ggsave 保存拼合的图
        ggsave(filename = out_file, plot = final_plot, width = 8, height = 8)
    }
    
    cat("Done! All plots saved to:", plot_dir, "\n")
    
} else {
    cat("No edSites passed the threshold (PP.H4 > 0.75).\n")
}
