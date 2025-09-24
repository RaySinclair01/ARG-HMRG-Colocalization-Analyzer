#!/bin/bash
# 文件名: create_quast_script.sh
# 用途: 创建质量评估脚本并执行，支持多个样本

# 定义样本列表（前缀，如A3A_）
SAMPLES=(
    "A1A" "A2A" "A3A"
    # 可以添加更多样本，例如: "A1A" "A1B" "A2A"
)

# 确保scripts目录存在
mkdir -p ./scripts

cat > ./scripts/run_quast.sh << EOF
#!/bin/bash

# 定义样本列表
SAMPLES=(
    "${SAMPLES[@]}"
)

# 创建必要的目录
mkdir -p ./assembly
mkdir -p ./logs

# 循环处理每个样本
for SAMPLE in "\${SAMPLES[@]}"; do
    # 定义输入文件，使用前缀加下划线
    INPUT_FASTA="./input_data/\${SAMPLE}_contigs.fasta"
    
    # 检查输入文件
    if [ ! -f "\$INPUT_FASTA" ]; then
        echo "错误：输入文件 \$INPUT_FASTA 不存在 (样本: \$SAMPLE)"
        echo "请确保组装步骤已成功完成并生成了\${SAMPLE}_contigs文件"
        continue  # 跳过当前样本，继续下一个
    fi
    
    # 定义输出目录和日志文件，使用前缀加下划线
    OUTPUT_DIR="./assembly/\${SAMPLE}_quast_report"
    LOG_FILE="./logs/\${SAMPLE}_quast.log"
    
    # 删除已存在的QUAST报告目录(如果存在)
    if [ -d "\$OUTPUT_DIR" ]; then
        echo "发现已存在的QUAST报告目录 (\$OUTPUT_DIR)，正在移除... (样本: \$SAMPLE)"
        rm -rf "\$OUTPUT_DIR"
    fi
    
    # 使用QUAST评估组装质量
    echo "评估组装质量 (样本: \$SAMPLE)..."
    quast.py "\${INPUT_FASTA}" \\
            -o "\$OUTPUT_DIR" \\
            --threads 74 \\
            > "\$LOG_FILE" 2>&1
    
    # 检查QUAST是否成功运行
    if [ -f "\$OUTPUT_DIR/report.txt" ]; then
        echo "组装质量评估完成 (样本: \$SAMPLE)，结果如下："
        cat "\$OUTPUT_DIR/report.txt"
    else
        echo "错误：QUAST评估似乎未成功完成 (样本: \$SAMPLE)"
        echo "请检查日志文件: \$LOG_FILE"
    fi
done

echo "所有样本的QUAST评估已完成。"
EOF

# 添加执行权限
chmod +x ./scripts/run_quast.sh

# 执行质量评估脚本
./scripts/run_quast.sh