#!/bin/bash

# --- 1. 【已修改】使用带有 "COLUMN_X" 表头的新文件 ---
EAS_FILE="/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/metaqtl/edqtl_nominal_EAS_final_cols.txt"
EUR_FILE="/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/metaqtl/edqtl_nominal_EUR_final_cols.txt"

# --- 2. 设置你的样本量 (N) ---
N_EAS=546
N_EUR=429  #

# --- 3. 设置工作目录 ---
BASE_DIR="/gpfs/hpc/home/chenchao/mayinuo/edqtl/meta/metaqtl"
TEMP_DIR="$BASE_DIR/metal_temp"
OUTPUT_DIR="$BASE_DIR/metal_results_hg38"

# --- 4. 创建目录 (如果它们还不存在) ---
mkdir -p "$TEMP_DIR"
mkdir -p "$OUTPUT_DIR"

echo "=== Meta-analysis V5.0 (Literal Column Header) Start ==="
echo "EAS File: $EAS_FILE (N=$N_EAS)"
echo "EUR File: $EUR_FILE (N=$N_EUR)"
echo "Out Dir:  $OUTPUT_DIR"
echo "-----------------------------------"

# --- 5. 生成唯一的 edSite 列表 ---
SITES_LIST_FILE="$TEMP_DIR/unique_edsites_hg38.txt"
echo "Generating unique edSite list..."
# 【已修改】awk 'NR > 1 {print $1}' -> 跳过我们刚添加的表头行 (NR > 1)
awk -F'\t' 'NR > 1 {print $1}' "$EAS_FILE" "$EUR_FILE" | sort -u > "$SITES_LIST_FILE"
SITE_COUNT=$(wc -l < "$SITES_LIST_FILE")
echo "Found $SITE_COUNT unique edSites to process."

CURRENT_SITE=0

# --- 6. 开始循环 ---
while read -r site; do
    if [ -z "$site" ]; then
        continue
    fi

    CURRENT_SITE=$((CURRENT_SITE + 1))
    echo "Processing site $CURRENT_SITE / $SITE_COUNT : $site"

    # --- 6a. 定义此循环的临时文件 ---
    EAS_TMP="$TEMP_DIR/${site}_eas.tmp"
    EUR_TMP="$TEMP_DIR/${site}_eur.tmp"
    METAL_SCRIPT="$TEMP_DIR/${site}_metal_script.txt"

    # --- 6b. 【已修改】awk 提取数据，并*保留表头* ---
    #     (我们把表头行 (NR==1) 复制过来，然后提取匹配的数据行)
    (head -n 1 "$EAS_FILE"; awk -F'\t' -v site="$site" '$1 == site' "$EAS_FILE") > "$EAS_TMP"
    (head -n 1 "$EUR_FILE"; awk -F'\t' -v site="$site" '$1 == site' "$EUR_FILE") > "$EUR_TMP"

    # --- 6c. 检查文件是否为空 (检查行数是否 > 1, 因为现在有表头) ---
    eas_exists=0
    eur_exists=0
    if [ $(wc -l < "$EAS_TMP") -gt 1 ]; then
        eas_exists=1
    fi
    if [ $(wc -l < "$EUR_TMP") -gt 1 ]; then
        eur_exists=1
    fi

    # 如果两个文件都只有表头 (或只有一个有数据)，则跳过
    if [ "$eas_exists" -eq 0 ] || [ "$eur_exists" -eq 0 ]; then
        echo "   -> Skipping (data missing in one or both populations)"
        rm -f "$EAS_TMP" "$EUR_TMP"
        continue
    fi

    # --- 6d. 【已修改】生成 METAL 脚本 (使用 "COLUMN_X" 表头名) ---
    cat > "$METAL_SCRIPT" << EOF
#
# METAL SCRIPT FOR: $site
#

# --- 选择方案：逆方差加权 (BETA/SE) ---
SCHEME STDERR

# --- 【已修改】使用我们新添加的 "COLUMN_X" 表头名 ---
MARKER     COLUMN_17
ALLELE     COLUMN_18 COLUMN_19
EFFECT     COLUMN_14
STDERR     COLUMN_13
PVALUE     COLUMN_12
FREQ       COLUMN_20

# --- 打开频率 QC 开关 ---
AVERAGEFREQ ON
MINMAXFREQ ON

# --- 处理 EAS 文件 ---
DEFAULT    $N_EAS
PROCESS    $EAS_TMP

# --- 处理 EUR 文件 ---
DEFAULT    $N_EUR
PROCESS    $EUR_TMP

# --- 定义输出文件 (无引号) ---
OUTFILE    $OUTPUT_DIR/meta_${site}_ .tbl

# --- 
ANALYZE HETEROGENEITY
QUIT
EOF

    # --- 6e. 运行 METAL ---
    metal "$METAL_SCRIPT"

done < "$SITES_LIST_FILE"

echo "-----------------------------------"
echo "All METAL jobs complete."
echo "Results are in: $OUTPUT_DIR"