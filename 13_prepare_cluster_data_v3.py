#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
共定位基因簇绘图数据格式化脚本 (最终版)

功能:
1. 读取由11号脚本生成的、已完全注释的详细共定位报告。
2. **不进行任何筛选**，将所有共定位contigs上的基因信息，
   严格按照用户指定的格式进行整理。
3. 输出的文件可以直接用于绘图软件或后续的筛选脚本。

使用方法:
1. 确保已安装 pandas 库: pip install pandas
2. 将此脚本与 'colocalization_analysis_final_v3' 文件夹放在同一目录下。
3. 在终端中运行此脚本: python 13_format_data_for_plotting.py
"""

import os
import pandas as pd
from typing import List

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
INPUT_FOLDER = "colocalization_analysis_final_v3"
OUTPUT_FOLDER = "cluster_plotting_data_final" # 新建一个文件夹存放整理好的数据

# 文件后缀定义
DETAILS_SUFFIX = "_colocalization_full_details.tsv"


def process_sample(sample_name: str):
    """处理单个样本，将其详细报告转换为绘图格式。"""
    print(f"--- 开始处理样本: {sample_name} ---")

    # 定义输入文件的路径
    details_file = os.path.join(INPUT_FOLDER, f"{sample_name}{DETAILS_SUFFIX}")

    if not os.path.exists(details_file):
        print(f"警告: 样本 {sample_name} 缺少详细报告文件，已跳过。")
        return

    try:
        # --- 1. 读取已经整合好的详细报告 ---
        # 这一步已经包含了所有需要的信息，无需再读取GFF或DIAMOND文件
        df_details = pd.read_csv(details_file, sep='\t')
        
        if df_details.empty:
            print(f"信息: 样本 {sample_name} 的详细报告为空，已跳过。")
            return

        # --- 2. 严格按照指定格式创建新的DataFrame ---
        # 'ID' 列使用 protein_id
        # 'source' 列使用 contig_id
        output_df = pd.DataFrame({
            'ID': df_details['protein_id'],
            'source': df_details['contig_id'],
            'start': df_details['start'],
            'end': df_details['end'],
            'strand': df_details['strand'],
            # 我们额外保留这两列，它们对于后续绘图标注颜色和名称至关重要
            'gene_type': df_details['gene_type'],
            'gene_name': df_details['gene_name']
        })

        # --- 3. 保存整理好的绘图数据 ---
        output_file = os.path.join(OUTPUT_FOLDER, f"{sample_name}_plot_data.tsv")
        output_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"成功: 已将用于绘图的完整基因簇数据保存至: {output_file}")

    except Exception as e:
        print(f"错误: 处理样本 {sample_name} 时发生严重错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主执行函数"""
    print("=" * 60)
    print("--- 整理用于基因簇可视化的数据 ---")
    print("=" * 60)

    if not os.path.isdir(INPUT_FOLDER):
        print(f"错误: 输入文件夹 '{INPUT_FOLDER}' 不存在！")
        return

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    samples: List[str] = sorted(
        [f.replace(DETAILS_SUFFIX, "") for f in os.listdir(INPUT_FOLDER) if f.endswith(DETAILS_SUFFIX)]
    )
    
    if not samples:
        print(f"错误: 在 '{INPUT_FOLDER}' 中没有找到任何详细报告文件。")
        return
        
    print(f"成功检测到 {len(samples)} 个样本: {', '.join(samples)}\n")

    for sample in samples:
        process_sample(sample)
        print("-" * 50)

    print("所有样本数据整理完毕！")
    print(f"所有格式化后的文件均已保存在 '{OUTPUT_FOLDER}' 文件夹中。")
    print("=" * 60)


if __name__ == "__main__":
    main()