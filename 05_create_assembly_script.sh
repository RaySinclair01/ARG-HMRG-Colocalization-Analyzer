#!/bin/bash
# 文件名: run_batch_assembly.sh
# 用途: 批量循环处理双端宏基因组数据，并使用 MEGAHIT 进行组装

# ========================== 用户配置区 ==========================
# 在下面的括号中，以空格分隔，添加或删除您需要处理的样本ID前缀。
# 例如: SAMPLES=("A1A" "A1B" "A2A" "A3A")
SAMPLES=(
    "A3A"
)
# ===============================================================

# --- 脚本主体 ---

# 首先确保所有需要的基础目录都存在
echo "正在创建所需目录: ./assembly, ./logs, ./input_data"
mkdir -p ./assembly
mkdir -p ./logs
mkdir -p ./input_data

# 遍历样本列表中的每一个样本ID
for SAMPLE_ID in "${SAMPLES[@]}"; do
    
    echo "=========================================================="
    echo "            开始处理样本: ${SAMPLE_ID}"
    echo "=========================================================="

    # --- 为当前样本定义变量 ---
    # 定义输出目录
    OUTPUT_DIR="./assembly/${SAMPLE_ID}_megahit_out"
    # 定义质控后的输入文件路径
    QC_R1="./qc_data/${SAMPLE_ID}_1.fq.clean.gz"
    QC_R2="./qc_data/${SAMPLE_ID}_2.fq.clean.gz"
    # 定义日志文件路径
    LOG_FILE="./logs/megahit_${SAMPLE_ID}.log"

    # 检查输入文件是否存在
    if [ ! -f "${QC_R1}" ] || [ ! -f "${QC_R2}" ]; then
        echo "❌ 错误: 找不到样本 ${SAMPLE_ID} 的输入文件!"
        echo "检查路径: ${QC_R1} 和 ${QC_R2}"
        continue # 跳过当前样本，继续处理下一个
    fi

    # 检查并移除已存在的输出目录 (避免MEGAHIT因目录已存在而报错)
    if [ -d "${OUTPUT_DIR}" ]; then
        echo "发现已存在的输出目录 ${OUTPUT_DIR}，正在移除..."
        rm -rf "${OUTPUT_DIR}"
    fi

    # --- 执行双端数据组装 ---
    echo "正在使用 MEGAHIT 对样本 ${SAMPLE_ID} 进行组装..."
    megahit -1 ${QC_R1} -2 ${QC_R2} \
            --num-cpu-threads 74 \
            --memory 0.9 \
            --min-contig-len 500 \
            --out-dir "${OUTPUT_DIR}" \
            --out-prefix L1HII0100002-TS_T1 \
            > ${LOG_FILE} 2>&1

    # --- 检查组装结果 ---
    FINAL_CONTIGS="${OUTPUT_DIR}/L1HII0100002-TS_T1.contigs.fa"
    if [ -f "${FINAL_CONTIGS}" ]; then
        echo "✅ 样本 ${SAMPLE_ID} 组装成功！"
        echo "Contigs 文件位于: ${FINAL_CONTIGS}"
        
        # 复制最终组装结果到统一的输入目录，并重命名
        cp "${FINAL_CONTIGS}" "./input_data/${SAMPLE_ID}_contigs.fasta"
        echo "最终结果已复制到 ./input_data/${SAMPLE_ID}_contigs.fasta"
    else
        echo "❌ 错误: 样本 ${SAMPLE_ID} 组装失败，未找到最终的 contigs 文件。"
        echo "请查看日志了解详情: ${LOG_FILE}"
    fi

    echo "样本 ${SAMPLE_ID} 处理完毕。"
    echo "" # 添加空行，方便区分不同样本的输出信息

done

echo "=========================================================="
echo "🎉 所有样本已处理完毕！"
echo "=========================================================="