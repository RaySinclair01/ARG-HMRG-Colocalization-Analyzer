#!/bin/bash
# 文件名: 08_create_database_script.sh
# 用途: 创建数据库准备脚本并执行。
# 更新: 自动下载并合并 BacMet 的“已验证”和“预测”数据库。

cat > ./scripts/prepare_databases.sh << 'EOF'
#!/bin/bash

# 如果任何命令失败，则立即退出
set -e

# --- 设置环境 ---
echo "创建数据库目录并进入..."
mkdir -p ~/databases
cd ~/databases

# --- 1. 处理CARD数据库 (此部分不变) ---
echo "--- 开始处理 CARD 数据库 ---"
# 检查是否已下载，避免重复下载
if [ ! -f "protein_fasta_protein_homolog_model.fasta" ]; then
    echo "下载并解压CARD数据库..."
    wget https://card.mcmaster.ca/latest/data -O card_data.tar.bz2
    tar -xjf card_data.tar.bz2
    rm card_data.tar.bz2
else
    echo "CARD 文件已存在，跳过下载。"
fi

echo "构建CARD Diamond数据库..."
diamond makedb --in protein_fasta_protein_homolog_model.fasta -d card_protein_homolog_db
echo "CARD 数据库准备完成。"
echo ""


# --- 2. 处理BacMet数据库 (新逻辑) ---
echo "--- 开始处理 BacMet 数据库 ---"
BACMET_EXP_FASTA="BacMet2_EXP_database.fasta"
BACMET_PRED_FASTA_GZ="BacMet2_predicted_database.fasta.gz"
BACMET_PRED_FASTA="BacMet2_predicted_database.fasta"
BACMET_COMBINED_FASTA="BacMet_combined.fasta"

# 下载 BacMet EXP (已验证) 数据库
echo "下载 BacMet (已验证) 数据库..."
wget http://bacmet.biomedicine.gu.se/download/${BACMET_EXP_FASTA}

# 下载 BacMet Predicted (预测) 数据库
echo "下载 BacMet (预测) 数据库..."
wget http://bacmet.biomedicine.gu.se/download/${BACMET_PRED_FASTA_GZ}

# 解压预测数据库
echo "解压 BacMet (预测) 数据库..."
gunzip -f ${BACMET_PRED_FASTA_GZ} # -f 强制解压，即使文件已存在

# 合并两个数据库
echo "合并已验证和预测的 BacMet 数据库..."
cat ${BACMET_EXP_FASTA} ${BACMET_PRED_FASTA} > ${BACMET_COMBINED_FASTA}

# 构建合并后的 BacMet Diamond 数据库
echo "构建合并后的 BacMet Diamond 数据库..."
# 注意：数据库名称已更新为 bacmet_combined_db
diamond makedb --in ${BACMET_COMBINED_FASTA} -d bacmet_combined_db

# 清理中间文件，节省空间
echo "清理中间文件..."
rm ${BACMET_EXP_FASTA} ${BACMET_PRED_FASTA}

echo "BacMet 数据库准备完成。"
echo ""


# --- 完成 ---
echo "所有数据库准备完成！"
echo "最终文件列表:"
ls -lah # 使用 -h 参数使文件大小更易读

# 返回原始工作目录
cd - > /dev/null
EOF

# 添加执行权限
chmod +x ./scripts/prepare_databases.sh

# 执行数据库准备脚本
./scripts/prepare_databases.sh