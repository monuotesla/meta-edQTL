############################################################
## 0. Load libraries
############################################################
library(data.table)
library(coloc)
library(ggplot2)
library(gridExtra)
library(dplyr)

############################################################
## 1. Parameters (ADHD)
############################################################
n_case    <- 38691
n_control <- 186843
n_total   <- 225534

window    <- 5e5
N_eqtl    <- 975
case_frac <- n_case / n_total   # ≈ 0.1716

############################################################
## 2. File paths
############################################################
# 显著 coloc 结果（循环依据）
sig_coloc_file <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/coloc/coloc1/ADHD/meta_edsite_coloc_sig_075.summary.txt"

# GWAS（已预处理）
gwas_file <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/coloc/gwas/ADHD_formatted.txt"

# edQTL 全量文件
eqtl_full_file <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/metaqtl/nominal_all/逆方差加权/metal_results_hg38/hg19_all/all_merged_classified.allele_fixed_hg19.txt"

# 输出目录
out_dir <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/coloc/coloc1/ADHD/sig_only_locus/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 3. Read data
############################################################
cat("Loading significant coloc edSites...\n")
sig_coloc <- fread(sig_coloc_file)
sig_edsites <- unique(sig_coloc$edsite)
cat("Number of significant edSites:", length(sig_edsites), "\n")

cat("Loading GWAS...\n")
gwas <- fread(gwas_file)

cat("Loading edQTL full data...\n")
eqtl_full <- fread(eqtl_full_file)

eqtl_full[, eGene := V1]
eqtl_full[, eSNP  := V2]
eqtl_full[, pval  := V11]

############################################################
## 4. Parse coordinates
############################################################
tmp_site <- tstrsplit(eqtl_full$eGene, "_")
eqtl_full[, site_chr := as.integer(tmp_site[[1]])]
eqtl_full[, site_pos := as.integer(tmp_site[[2]])]

tmp_snp <- tstrsplit(eqtl_full$eSNP, "_")
eqtl_full[, snp_chr := as.integer(tmp_snp[[1]])]
eqtl_full[, snp_pos := as.integer(tmp_snp[[2]])]

############################################################
## 5. Function: coloc + SNP table + locus plot
############################################################
run_sig_coloc_with_plot <- function(site_id) {

  cat("Processing:", site_id, "\n")

  edinfo <- eqtl_full[eGene == site_id][1]
  if (is.na(edinfo$site_chr)) return(NULL)

  chr_k <- edinfo$site_chr
  pos_k <- edinfo$site_pos

  eqtl_reg <- eqtl_full[
    eGene == site_id &
    snp_chr == chr_k &
    snp_pos %between% c(pos_k - window, pos_k + window)
  ]
  if (nrow(eqtl_reg) < 2) return(NULL)

  gwas_reg <- gwas[
    chr == chr_k &
    pos %between% c(pos_k - window, pos_k + window)
  ]
  if (nrow(gwas_reg) < 5) return(NULL)

  snps_common <- intersect(gwas_reg$snp, eqtl_reg$eSNP)
  if (length(snps_common) < 2) return(NULL)

  gwas_use <- gwas_reg[match(snps_common, snp)]
  eqtl_use <- eqtl_reg[match(snps_common, eSNP)]

  coloc_res <- coloc.abf(
    dataset1 = list(
      pvalues = gwas_use$pval,
      type    = "cc",
      s       = case_frac,
      N       = n_total,
      snp     = snps_common
    ),
    dataset2 = list(
      pvalues = eqtl_use$pval,
      type    = "quant",
      N       = N_eqtl,
      snp     = snps_common
    ),
    MAF = gwas_use$maf
  )

  ##########################################################
  ## SNP-level table
  ##########################################################
  snp_df <- data.frame(
    SNP    = snps_common,
    chr    = chr_k,
    pos    = gwas_use$pos,
    GWAS_P = gwas_use$pval,
    EQTL_P = eqtl_use$pval,
    PP_H4  = coloc_res$results$SNP.PP.H4
  )

  snp_df$PP_H4[is.na(snp_df$PP_H4)] <- 0
  snp_df$BP_Mb <- snp_df$pos / 1e6
  snp_df <- snp_df[order(-snp_df$PP_H4), ]

  fwrite(
    snp_df,
    file = paste0(out_dir, "/", site_id, "_snp_level.txt"),
    sep = "\t"
  )

  ##########################################################
  ## Locus plot
  ##########################################################
  top_snp <- snp_df[1, ]

  p1 <- ggplot(snp_df, aes(BP_Mb, -log10(GWAS_P))) +
    geom_point(aes(color = PP_H4, size = PP_H4), alpha = 0.7) +
    geom_point(data = top_snp, shape = 18, size = 5, color = "red") +
    geom_text(data = top_snp, aes(label = SNP),
              vjust = -1.3, size = 3) +
    scale_color_gradient(low = "grey70", high = "red") +
    scale_size(range = c(1.5, 4), guide = "none") +
    labs(title = paste(site_id, "GWAS"),
         y = "-log10(P)", x = NULL) +
    theme_bw()

  p2 <- ggplot(snp_df, aes(BP_Mb, -log10(EQTL_P))) +
    geom_point(aes(color = PP_H4, size = PP_H4), alpha = 0.7) +
    geom_point(data = top_snp, shape = 18, size = 5, color = "red") +
    geom_text(data = top_snp, aes(label = SNP),
              vjust = -1.3, size = 3) +
    scale_color_gradient(low = "grey70", high = "red") +
    scale_size(range = c(1.5, 4), guide = "none") +
    labs(title = paste(site_id, "edQTL"),
         y = "-log10(P)", x = "Position (Mb)") +
    theme_bw()

  ggsave(
    filename = paste0(out_dir, "/", site_id, "_locus.pdf"),
    plot = grid.arrange(p1, p2, ncol = 1),
    width = 8, height = 8
  )

  data.frame(
    edsite = site_id,
    nsnps  = coloc_res$summary["nsnps"],
    PP.H4  = coloc_res$summary["PP.H4.abf"]
  )
}

############################################################
## 6. Run for ALL significant edSites
############################################################
res_list <- lapply(sig_edsites, run_sig_coloc_with_plot)
res_list <- res_list[!sapply(res_list, is.null)]

summary_df <- rbindlist(res_list)

############################################################
## 7. Save summary
############################################################
fwrite(
  summary_df,
  file = paste0(out_dir, "/sig_coloc_recomputed.summary.txt"),
  sep = "\t"
)

cat("Finished. All results saved in:\n", out_dir, "\n")



