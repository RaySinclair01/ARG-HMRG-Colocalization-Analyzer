#!/bin/bash
# 文件名: create_prodigal_script.sh
# 用途: 创建一个支持并行处理的、经过优化的基因预测脚本并执行。
# 版本: 3.0 (集成并行化处理与并发任务控制)

# 定义样本列表（前缀）
SAMPLES=(
    "A1A"
    "A2A"
    "A3A"
    "A1B"
    # 将来可以添加更多样本，脚本会自动按设定的并行数处理
)

# 确保scripts目录存在
mkdir -p ./scripts

# --- 生成新的、支持并行的基因预测脚本 ---
cat > ./scripts/run_prodigal.sh << 'EOF'
#!/bin/bash
# 如果任何命令在子shell中失败，则立即退出该子shell
set -e

# --- 1. 定义变量和环境 ---

# 定义样本列表 (此列表将由外部脚本动态注入)
SAMPLES=(
    "A1A" "A2A" "A3A" "A1B"
)

# --- 并行化配置 ---
# 设置最大并行任务数。对于4个样本，设置为4即可全部同时运行。
# 如果您有20个样本，但只想同时运行8个，可将此值改为8。
MAX_JOBS=4

# 定义基础目录
BASE_INPUT_DIR="./input_data"
BASE_OUTPUT_DIR="./analysis_results"
LOG_DIR="./logs"

# 确保输出目录和日志目录存在
mkdir -p "${BASE_OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"


# --- 2. 循环启动并行任务 ---

echo "开始并行基因预测流程，最大并行任务数: ${MAX_JOBS}"

for SAMPLE in "${SAMPLES[@]}"; do
    # 将每个样本的处理逻辑放入一个子shell '( ... )' 并在后台 '&' 运行
    (
        echo "--- [启动] 正在处理样本: ${SAMPLE} ---"
        
        # 定义完整的文件路径
        INPUT_FASTA="${BASE_INPUT_DIR}/${SAMPLE}_contigs.fasta"
        OUTPUT_GFF="${BASE_OUTPUT_DIR}/${SAMPLE}_predicted_genes.gff"
        OUTPUT_FAA="${BASE_OUTPUT_DIR}/${SAMPLE}_predicted_proteins.faa"
        OUTPUT_GENES="${BASE_OUTPUT_DIR}/${SAMPLE}_predicted_genes.fna"
        LOG_FILE="${LOG_DIR}/${SAMPLE}_prodigal.log"
        
        # 检查输入文件是否存在且非空
        if [ ! -s "${INPUT_FASTA}" ]; then
            echo "--- [警告] 输入文件 ${INPUT_FASTA} 不存在或为空 (样本: ${SAMPLE})。跳过此样本。 ---"
            exit 0 # 正常退出此子shell，不影响其他任务
        fi
        
        # 执行基因预测
        prodigal -i "${INPUT_FASTA}" \
                 -o "${OUTPUT_GFF}" \
                 -a "${OUTPUT_FAA}" \
                 -d "${OUTPUT_GENES}" \
                 -p meta -f gff -q \
                 > "${LOG_FILE}" 2>&1
        
        # 检查结果并进行清理
        if [ -s "${OUTPUT_FAA}" ]; then
            # 关键修复：清理蛋白质FASTA文件中的非法字符
            sed -i '/^>/! s/;//g' "${OUTPUT_FAA}"
            
            GENE_COUNT=$(grep -c ">" "${OUTPUT_FAA}")
            echo "--- [成功] 样本 ${SAMPLE} 预测到 ${GENE_COUNT} 个基因。 ---"
        else
            echo "--- [错误] Prodigal 未能为样本 ${SAMPLE} 生成输出。请检查日志: ${LOG_FILE} ---"
        fi
    ) & # & 将子shell置于后台执行

    # --- 并发任务数量控制 ---
    # 当后台运行的任务数量达到上限时，暂停并等待任一任务完成后再继续
    # bash 4.3+ 支持 wait -n, 兼容性更好
    if [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; then
        wait -n
    fi

done

# --- 3. 等待所有任务完成 ---

echo ""
echo "所有样本任务已启动。正在等待剩余的后台任务全部完成..."
# wait 命令会阻塞在这里，直到所有后台任务都结束
wait

echo "----------------------------------------------------"
echo "所有样本的并行基因预测流程已全部完成。"
EOF

# 使用sed将最新的SAMPLES数组注入到生成的脚本中
SAMPLES_STRING=$(printf "\"%s\" " "${SAMPLES[@]}")
sed -i "s|\"A1A\" \"A2A\" \"A3A\" \"A1B\"|${SAMPLES_STRING}|" ./scripts/run_prodigal.sh


# 添加执行权限
chmod +x ./scripts/run_prodigal.sh

# 执行基因预测脚本
echo "正在执行生成的并行脚本 ./scripts/run_prodigal.sh..."
./scripts/run_prodigal.sh