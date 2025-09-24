#!/bin/bash
# 文件名: 10_create_mge_annotation_script.sh
# 用途: 创建一个经过优化的、支持多样本的BacMet注释脚本并执行。
# 版本: 2.0 (与ARG注释脚本结构保持一致)

# --- 1. 用户配置区 (将重要前缀和参数提到前面) ---

# 定义需要处理的样本列表 (应与ARG注释脚本的列表保持一致)
SAMPLES=(
    "A1A"
    "A2A"
    "A3A"
    "A1B"
)

# 定义数据库路径 (只需修改这里)
BACMET_DB_PATH="$HOME/databases/bacmet_combined_db"

# 定义DIAMOND比对的核心参数 (方便统一调整)
THREADS=16         # 使用的CPU线程数
MAX_TARGET_SEQS=1  # 每个查询序列报告的最佳匹配数


# --- 2. 脚本生成区 ---

# 确保scripts目录存在
mkdir -p ./scripts

# 使用 cat 生成功能脚本
cat > ./scripts/run_bacmet_annotation.sh << 'EOF'
#!/bin/bash
# 如果任何命令失败，则立即退出脚本
set -e

# --- A. 配置区 (这些值将由外部脚本动态注入) ---
SAMPLES=( "A1A" "A2A" "A3A" "A1B" )
BACMET_DB_PATH="$HOME/databases/bacmet_combined_db"
THREADS=16
MAX_TARGET_SEQS=1

# --- B. 环境准备 ---
BASE_OUTPUT_DIR="./analysis_results"
LOG_DIR="./logs"

mkdir -p "${BASE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# 检查数据库文件是否存在 (只需检查一次)
echo "--- 检查 BacMet (合并版) 数据库文件 ---"
if [ -f "${BACMET_DB_PATH}.dmnd" ]; then
    echo "成功: 找到数据库文件: ${BACMET_DB_PATH}.dmnd"
else
    echo "错误: 未找到数据库文件: ${BACMET_DB_PATH}.dmnd"
    echo "请确认数据库构建脚本 (08) 已成功运行。"
    exit 1
fi
echo ""

# --- C. 循环处理每个样本 ---
for SAMPLE in "${SAMPLES[@]}"; do
    echo "----------------------------------------------------"
    echo "--- 开始处理样本: ${SAMPLE} ---"
    
    # 动态定义该样本的文件路径
    QUERY_FAA="${BASE_OUTPUT_DIR}/${SAMPLE}_predicted_proteins.faa"
    OUTPUT_BACMET_HITS="${BASE_OUTPUT_DIR}/${SAMPLE}_bacmet_hits.tsv"
    LOG_FILE="${LOG_DIR}/${SAMPLE}_diamond_bacmet.log"

    # 1. 检查输入文件
    echo "检查输入蛋白文件: ${QUERY_FAA}"
    if [ ! -s "${QUERY_FAA}" ]; then
        echo "警告: 未找到输入文件或文件为空。将跳过此样本。"
        continue
    fi
    
    # 注意: Prodigal脚本已包含清理步骤，此处无需重复清理。
    # 如果输入文件来源不确定，可以像ARG脚本一样加入清理步骤。

    # 2. 执行Diamond比对
    echo "正在运行 DIAMOND 进行 BacMet 注释..."
    diamond blastp \
        --db "${BACMET_DB_PATH}.dmnd" \
        --query "${QUERY_FAA}" \
        --out "${OUTPUT_BACMET_HITS}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --threads ${THREADS} \
        --max-target-seqs ${MAX_TARGET_SEQS} \
        --sensitive \
        > "${LOG_FILE}" 2>&1

    # 3. 检查并报告结果
    echo "检查输出结果..."
    if [ -s "${OUTPUT_BACMET_HITS}" ]; then
        MATCH_COUNT=$(wc -l < "${OUTPUT_BACMET_HITS}")
        echo "成功: 为样本 ${SAMPLE} 找到 ${MATCH_COUNT} 个潜在的重金属抗性基因匹配。"
        echo "结果文件: ${OUTPUT_BACMET_HITS}"
    else
        echo "警告: 未在样本 ${SAMPLE} 中找到任何匹配项。"
        echo "详情请查看日志: ${LOG_FILE}"
    fi
done

echo "----------------------------------------------------"
echo "所有样本的 BacMet 注释流程已全部完成。"
EOF


# --- 3. 动态注入配置到生成的脚本中 ---

# 使用 sed 安全地替换脚本中的配置项
SAMPLES_STRING=$(printf "\"%s\" " "${SAMPLES[@]}")
sed -i "s|SAMPLES=(.*)|SAMPLES=( ${SAMPLES_STRING})|" ./scripts/run_bacmet_annotation.sh
sed -i "s|BACMET_DB_PATH=.*|BACMET_DB_PATH=\"${BACMET_DB_PATH}\"|" ./scripts/run_bacmet_annotation.sh
sed -i "s|THREADS=.*|THREADS=${THREADS}|" ./scripts/run_bacmet_annotation.sh
sed -i "s|MAX_TARGET_SEQS=.*|MAX_TARGET_SEQS=${MAX_TARGET_SEQS}|" ./scripts/run_bacmet_annotation.sh


# --- 4. 执行脚本 ---

# 添加执行权限
chmod +x ./scripts/run_bacmet_annotation.sh

# 执行MGE注释脚本
echo "正在执行生成的 BacMet 注释脚本 ./scripts/run_bacmet_annotation.sh..."
./scripts/run_bacmet_annotation.sh