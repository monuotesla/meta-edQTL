#!/usr/bin/bash

#==========#
# Get args #
#==========#
echo "Script started at: $(date '+%Y-%m-%d %H:%M:%S')"
#gwas_file="/gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/GWAS/SCZ/cleaned_SCZgwas_format.coloc"
gwas_file="/gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/GWAS/PD.GWAS_discovery.summaryIDmaf"
gwas_sample_size=4450
gwas_cases=2478
#gwas_cases=22778

# List of eQTL datasets and their respective sample sizes
eqtl_files=("/gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/coloc/eqtl/eqtl.nominalhg19.prenatal.chrall.sig"
"/gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/coloc/eqtl/eqtl.nominalhg19.adult.chrall.sig"
"/gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/coloc/eqtl/eqtl.nominalhg19.aged.chrall.sig")
eqtl_sample_sizes=(107 143 403)  # Example sample sizes

# Create output directory
output_dir="coloc_results"
mkdir -p $output_dir
chmod u+w $output_dir 
# Function for coloc analysis
coloc_analysis() {
  eqtl_file=$1
  eqtl_sample_size=$2
  gwas_file=$3
  gwas_sample_size=$4
  gwas_cases=$5
  tissue=$(basename $eqtl_file | cut -d'.' -f3)

  Rscript - <<EOF
  library(coloc)
  library(data.table)
  eqtl <- fread("${eqtl_file}", head=F)
  gwas <- fread("${gwas_file}", head=F)
  # Rename columns for easier access
  gwas\$rs_id <- gwas\$V10
  gwas\$variant_id <- gwas\$V10
  gwas\$pval_nominal <- gwas\$V8
  gwas\$varbeta <- gwas\$V7
  gwas\$maf <- ifelse(gwas\$V11 > 0.5, 1 - gwas\$V11, gwas\$V11)
  print(head(gwas\$maf))
  eqtl\$variant_id <- eqtl\$V8
  eqtl\$gene_id <- eqtl\$V1
  eqtl\$pval_nominal <- eqtl\$V12
  eqtl\$beta <- eqtl\$V14
  eqtl\$varbeta <- eqtl\$V15
  eqtl\$rs_id <- eqtl\$variant_id
  print(head(eqtl))
  # Merge datasets by SNPs
  input <- merge(eqtl, gwas, by="rs_id", all=FALSE, suffixes=c("_eqtl", "_gwas"))
  print(head(input))
  # Remove duplicate SNPs
  input <- input[!duplicated(input\$rs_id),]
  print(head(input))
  res0 <- matrix(NA, ncol=6, nrow=length(unique(input\$gene_id)))
  res0 <- as.data.frame(res0)
  rownames(res0) <- unique(input\$gene_id)

  for (i in unique(input\$gene_id)) {
    result <- coloc.abf(dataset1 = list(snp = input[input\$gene_id == i,]\$variant_id_gwas,
                                    pvalues = input[input\$gene_id == i,]\$pval_nominal_gwas,
                                    type = "cc",
                                    N = ${gwas_sample_size},
                                    s = ${gwas_cases} / ${gwas_sample_size}),
                    dataset2 = list(snp = input[input\$gene_id == i,]\$variant_id_eqtl,
                                    pvalues = input[input\$gene_id == i,]\$pval_nominal_eqtl,
                                    type = "quant",
                                    N = ${eqtl_sample_size}),
                    MAF = input[input\$gene_id == i,]\$maf)
  res0[i, ] <- as.numeric(result\$summary)
  }
  colnames(res0) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
  write.table(res0, file="Edu.$tissue.eqtl.pp.coloc", sep="\t", row.names=T, quote=F, col.names=NA)
EOF
}

# Export the function to use in parallel
export -f coloc_analysis

# Run coloc analysis in paral
parallel --tmpdir /gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/coloc/abf/tmp coloc_analysis ::: "${eqtl_files[@]}" ::: "${eqtl_sample_sizes[@]}" ::: "${gwas_file}" ::: "${gwas_sample_size}" ::: "${gwas_cases}"

echo "Script completed at: $(date '+%Y-%m-%d %H:%M:%S')"
