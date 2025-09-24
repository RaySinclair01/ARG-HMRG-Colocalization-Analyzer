#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
共定位基因簇绘图数据格式化脚本 (增加去重功能)

功能:
1. 读取由14号筛选脚本生成的Top N contigs的绘图数据。
2. 将数据格式严格转换为示例文件 'gene_cluster_gap.tsv' 的格式。
3. **新增：在保存前去除完全重复的行，确保数据唯一性。**
4. 为每个样本生成一个可以直接上传进行可视化的TSV文件。

使用方法:
1. 确保已安装 pandas 库: pip install pandas
2. 将此脚本与 'top_10_contigs_for_plotting' 文件夹放在同一目录下。
3. 在终端中运行此脚本: python 15_format_for_external_plotter.py
"""

import os
import pandas as pd
from typing import List

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
INPUT_FOLDER = "top_10_contigs_for_plotting"
OUTPUT_FOLDER = "final_data_for_visualization"

# 文件后缀定义
INPUT_SUFFIX = "_plot_data.tsv"


def process_file_for_plotting(sample_name: str):
    """处理单个样本文件，将其转换为最终的绘图格式。"""
    print(f"--- 正在处理样本: {sample_name} ---")
    
    input_file = None
    for fname in os.listdir(INPUT_FOLDER):
        if fname.startswith(sample_name) and fname.endswith("_plot_data.tsv"):
            input_file = os.path.join(INPUT_FOLDER, fname)
            break
            
    if not input_file:
        print(f"警告: 未找到样本 {sample_name} 的绘图数据文件，已跳过。")
        return

    try:
        df = pd.read_csv(input_file, sep='\t')
        
        if df.empty:
            print(f"信息: 样本 {sample_name} 的绘图数据文件为空，已跳过。")
            return

        df_formatted = df[['ID', 'source', 'start', 'end', 'strand']].copy()
        
        unique_contigs = df_formatted['source'].unique()
        contig_map = {contig_id: f"Contig_{i+1}" for i, contig_id in enumerate(unique_contigs)}
        df_formatted['source'] = df_formatted['source'].map(contig_map)
        
        # --- *** 新增的去重处理 *** ---
        # 记录去重前的行数
        rows_before = len(df_formatted)
        # 根据所有列去除完全重复的行
        df_formatted.drop_duplicates(inplace=True)
        # 记录去重后的行数
        rows_after = len(df_formatted)
        
        if rows_before > rows_after:
            print(f"  > 执行了去重处理，移除了 {rows_before - rows_after} 个重复行。")
        else:
            print("  > 数据检查完毕，未发现重复行。")

        # --------------------------------

        output_file = os.path.join(OUTPUT_FOLDER, f"{sample_name}_final_plot_data.tsv")
        df_formatted.to_csv(output_file, sep='\t', index=False)
        
        print(f"成功: 已将最终绘图数据保存至: {output_file}")
        print(f"此文件包含 {len(unique_contigs)} 个基因簇 (contigs) 的 {rows_after} 条基因记录。")

    except Exception as e:
        print(f"错误: 处理样本 {sample_name} 时发生严重错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主执行函数"""
    print("=" * 60)
    print("--- 格式化数据以适配外部基因簇绘图工具 (带去重) ---")
    print("=" * 60)

    if not os.path.isdir(INPUT_FOLDER):
        print(f"错误: 输入文件夹 '{INPUT_FOLDER}' 不存在！请确保已运行之前的筛选脚本。")
        return

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    samples: List[str] = sorted(list(set(
        fname.split('_')[0] for fname in os.listdir(INPUT_FOLDER)
        if fname.endswith(INPUT_SUFFIX)
    )))
    
    if not samples:
        print(f"错误: 在 '{INPUT_FOLDER}' 中没有找到任何绘图数据文件。")
        return
        
    print(f"成功检测到 {len(samples)} 个样本: {', '.join(samples)}\n")

    for sample in samples:
        process_file_for_plotting(sample)
        print("-" * 50)

    print("所有样本数据格式化完毕！")
    print(f"所有最终用于绘图的文件均已保存在 '{OUTPUT_FOLDER}' 文件夹中。")
    print("=" * 60)


if __name__ == "__main__":
    main()