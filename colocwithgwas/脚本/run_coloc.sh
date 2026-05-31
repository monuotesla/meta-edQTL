#!/bin/bash
#SBATCH --job-name=coloc_analysis
#SBATCH --output=coloc_analysis.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12  # 使用12个CPU核
#SBATCH --mem=64G           # 分配64G内存

sh /gpfs/hpc/home/chenchao/mayinuo/edqtl/coloc/coloc_meta.sh # 运行R脚本
