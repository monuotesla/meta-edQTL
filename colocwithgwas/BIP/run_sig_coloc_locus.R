############################
## 0. Load libraries
############################
library(data.table)
library(coloc)
library(ggplot2)
library(dplyr)
library(gridExtra)

############################
## 1. File paths
############################

# 显著共定位结果（你之前保存的）
sig_file <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/coloc/coloc1/BIP_EUR/meta_edsite_coloc_sig_075.summary.txt"

# GWAS summary（已预处理）
gwas_file <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/coloc/gwas/BIP_formatted_final.txt"

# edQTL 全量文件
eqtl_file <- "/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/metaqtl/nominal_all/逆方差加权/metal_results_hg38/hg19_all/all_merged_classified.allele_fixed_hg19.txt"

# 输出目录
out_dir <- "Sig_Coloc_Locus"
dir.create(out_dir, showWarnings = FALSE)

############################
## 2. Global parameters
############################
window    <- 5e5
N_eqtl    <- 975
N_case    <- 41917
case_frac <- 0.101

############################
## 3. Read data
############################

cat("Loading GWAS...\n")
gwas <- fread(gwas_file)

cat("Loading edQTL full data...\n")
eqtl_full <- fread(eqtl_file)

eqtl_full[, eGene := V1]
eqtl_full[, eSNP  := V2]
eqtl_full[, pval  := V11]

############################
## 4. Parse coordinates
############################
tmp_site <- tstrsplit(eqtl_full$eGene, "_")
eqtl_full[, site_chr := as.integer(tmp_site[[1]])]
eqtl_full[, site_pos := as.integer(tmp_site[[2]])]

tmp_snp <- tstrsplit(eqtl_full$eSNP, "_")
eqtl_full[, snp_chr := as.integer(tmp_snp[[1]])]
eqtl_full[, snp_pos := as.integer(tmp_snp[[2]])]

############################
## 5. Read significant edSites
############################
sig_summary <- fread(sig_file)
sig_edsites <- unique(sig_summary$edsite)

cat("Number of significant edSites:", length(sig_edsites), "\n")
print(sig_edsites)

############################
## 6. Coloc + locus function
############################
run_coloc_and_plot <- function(gene) {

  edinfo <- eqtl_full[eGene == gene][1]
  if (is.na(edinfo$site_chr)) return(NULL)

  chr_k <- edinfo$site_chr
  pos_k <- edinfo$site_pos

  eqtl_reg <- eqtl_full[
    eGene == gene &
    snp_chr == chr_k &
    snp_pos %between% c(pos_k - window, pos_k + window)
  ]

  gwas_reg <- gwas[
    chr == chr_k &
    pos %between% c(pos_k - window, pos_k + window)
  ]

  snps_common <- intersect(gwas_reg$snp, eqtl_reg$eSNP)
  if (length(snps_common) < 5) return(NULL)

  gwas_use <- gwas_reg[match(snps_common, snp)]
  eqtl_use <- eqtl_reg[match(snps_common, eSNP)]

  coloc_res <- coloc.abf(
    dataset1 = list(
      pvalues = gwas_use$pval,
      type = "cc",
      s = case_frac,
      N = N_case,
      snp = snps_common
    ),
    dataset2 = list(
      pvalues = eqtl_use$pval,
      type = "quant",
      N = N_eqtl,
      snp = snps_common
    ),
    MAF = gwas_use$maf
  )

  snp_df <- data.frame(
    edsite = gene,
    SNP    = snps_common,
    BP     = gwas_use$pos,
    GWAS_P = gwas_use$pval,
    EQTL_P = eqtl_use$pval,
    PP_H4  = coloc_res$results$SNP.PP.H4
  )

  snp_df$BP_Mb <- snp_df$BP / 1e6
  snp_df$PP_H4[is.na(snp_df$PP_H4)] <- 0
  snp_df <- snp_df[order(-snp_df$PP_H4), ]

  ################ plot ################
  top_snp <- snp_df[1, ]

  p1 <- ggplot(snp_df, aes(BP_Mb, -log10(GWAS_P))) +
    geom_point(aes(color = PP_H4, size = PP_H4), alpha = 0.7) +
    geom_point(data = top_snp, shape = 18, size = 5, color = "red") +
    geom_text(data = top_snp, aes(label = SNP),
              vjust = -1.3, size = 3) +
    scale_color_gradient(low = "grey70", high = "red") +
    scale_size(range = c(1.5, 4), guide = "none") +
    labs(title = paste(gene, "GWAS"),
         y = "-log10(P)", x = NULL) +
    theme_bw()

  p2 <- ggplot(snp_df, aes(BP_Mb, -log10(EQTL_P))) +
    geom_point(aes(color = PP_H4, size = PP_H4), alpha = 0.7) +
    geom_point(data = top_snp, shape = 18, size = 5, color = "red") +
    geom_text(data = top_snp, aes(label = SNP),
              vjust = -1.3, size = 3) +
    scale_color_gradient(low = "grey70", high = "red") +
    scale_size(range = c(1.5, 4), guide = "none") +
    labs(title = paste(gene, "edQTL"),
         y = "-log10(P)", x = "Position (Mb)") +
    theme_bw()

  ggsave(
    filename = paste0(out_dir, "/", gene, "_locus.pdf"),
    plot = grid.arrange(p1, p2, ncol = 1),
    width = 8, height = 8
  )

  fwrite(
    snp_df,
    file = paste0(out_dir, "/", gene, "_snp_level.txt"),
    sep = "\t"
  )
}

############################
## 7. Run all significant edSites
############################
for (gene in sig_edsites) {
  cat("Processing:", gene, "\n")
  try(run_coloc_and_plot(gene), silent = TRUE)
}

cat("All done. Results saved in:", out_dir, "\n")
