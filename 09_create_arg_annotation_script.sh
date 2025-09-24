#!/bin/bash
# 文件名: 09_create_arg_annotation_script.sh
# 用途: 创建一个经过优化的、支持多样本的ARG注释脚本并执行。
# 版本: 2.0 (集成多样本循环和顶部参数配置)

# --- 1. 用户配置区 (将重要前缀和参数提到前面) ---

# 定义需要处理的样本列表
SAMPLES=(
    "A1A"
    "A2A"
    "A3A"
    "A1B"
)

# 定义数据库路径 (只需修改这里)
CARD_DB_PATH="$HOME/databases/card_protein_homolog_db"

# 定义DIAMOND比对的核心参数 (方便统一调整)
THREADS=74         # 使用的CPU线程数
EVALUE="1e-5"      # E-value 阈值
QUERY_COVER=50     # 查询序列覆盖度阈值 (%)
IDENTITY=40        # 一致性阈值 (%)


# --- 2. 脚本生成区 ---

# 确保scripts目录存在
mkdir -p ./scripts

# 使用 cat 生成功能脚本
cat > ./scripts/run_arg_annotation.sh << 'EOF'
#!/bin/bash
# 如果任何命令失败，则立即退出脚本
set -e

# --- A. 配置区 (这些值将由外部脚本动态注入) ---
SAMPLES=( "A1A" "A2A" "A3A" "A1B" )
CARD_DB_PATH="$HOME/databases/card_protein_homolog_db"
THREADS=74
EVALUE="1e-5"
QUERY_COVER=50
IDENTITY=40

# --- B. 环境准备 ---
BASE_OUTPUT_DIR="./analysis_results"
LOG_DIR="./logs"

mkdir -p "${BASE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# 检查数据库文件是否存在 (只需检查一次)
echo "--- 检查 CARD 数据库文件 ---"
if [ -f "${CARD_DB_PATH}.dmnd" ]; then
    echo "成功: 找到数据库文件: ${CARD_DB_PATH}.dmnd"
else
    echo "错误: 未找到数据库文件: ${CARD_DB_PATH}.dmnd"
    echo "请确认数据库构建步骤已成功完成。"
    exit 1
fi
echo ""

# --- C. 循环处理每个样本 ---
for SAMPLE in "${SAMPLES[@]}"; do
    echo "----------------------------------------------------"
    echo "--- 开始处理样本: ${SAMPLE} ---"
    
    # 动态定义该样本的文件路径
    QUERY_FAA="${BASE_OUTPUT_DIR}/${SAMPLE}_predicted_proteins.faa"
    OUTPUT_ARG_M8="${BASE_OUTPUT_DIR}/${SAMPLE}_card_hits.m8"
    LOG_FILE="${LOG_DIR}/${SAMPLE}_diamond_card.log"

    # 1. 检查并清理输入文件
    echo "检查输入蛋白文件: ${QUERY_FAA}"
    if [ ! -s "${QUERY_FAA}" ]; then
        echo "警告: 未找到输入文件或文件为空。将跳过此样本。"
        continue
    fi
    
    # 使用临时文件进行清理，避免污染原始预测结果
    CLEAN_QUERY_FAA=$(mktemp)
    awk '/^>/ {print; next} {gsub(/[^A-Z*]/, "X"); print}' "${QUERY_FAA}" > "${CLEAN_QUERY_FAA}"

    # 2. 执行Diamond比对
    echo "正在运行 DIAMOND 进行 ARG 注释..."
    diamond blastp \
        --db "${CARD_DB_PATH}.dmnd" \
        --query "${CLEAN_QUERY_FAA}" \
        --out "${OUTPUT_ARG_M8}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp \
        --threads ${THREADS} \
        --evalue ${EVALUE} \
        --query-cover ${QUERY_COVER} \
        --id ${IDENTITY} \
        --max-target-seqs 5 \
        --ultra-sensitive \
        > "${LOG_FILE}" 2>&1
    
    # 清理临时文件
    rm "${CLEAN_QUERY_FAA}"

    # 3. 检查并报告结果
    echo "检查输出结果..."
    if [ -s "${OUTPUT_ARG_M8}" ]; then
        MATCH_COUNT=$(wc -l < "${OUTPUT_ARG_M8}")
        echo "成功: 为样本 ${SAMPLE} 找到 ${MATCH_COUNT} 个潜在的 ARG 匹配。"
        echo "结果文件: ${OUTPUT_ARG_M8}"
    else
        echo "警告: 未在样本 ${SAMPLE} 中找到任何 ARG 匹配。"
        echo "这可能是正常的生物学结果，或参数过于严格。详情请查看日志: ${LOG_FILE}"
    fi
done

echo "----------------------------------------------------"
echo "所有样本的 ARG 注释流程已全部完成。"
EOF


# --- 3. 动态注入配置到生成的脚本中 ---

# 使用 sed 安全地替换脚本中的配置项
SAMPLES_STRING=$(printf "\"%s\" " "${SAMPLES[@]}")
sed -i "s|SAMPLES=(.*)|SAMPLES=( ${SAMPLES_STRING})|" ./scripts/run_arg_annotation.sh
sed -i "s|CARD_DB_PATH=.*|CARD_DB_PATH=\"${CARD_DB_PATH}\"|" ./scripts/run_arg_annotation.sh
sed -i "s|THREADS=.*|THREADS=${THREADS}|" ./scripts/run_arg_annotation.sh
sed -i "s|EVALUE=.*|EVALUE=\"${EVALUE}\"|" ./scripts/run_arg_annotation.sh
sed -i "s|QUERY_COVER=.*|QUERY_COVER=${QUERY_COVER}|" ./scripts/run_arg_annotation.sh
sed -i "s|IDENTITY=.*|IDENTITY=${IDENTITY}|" ./scripts/run_arg_annotation.sh


# --- 4. 执行脚本 ---

# 添加执行权限
chmod +x ./scripts/run_arg_annotation.sh

# 执行ARG注释脚本
echo "正在执行生成的 ARG 注释脚本 ./scripts/run_arg_annotation.sh..."
./scripts/run_arg_annotation.sh