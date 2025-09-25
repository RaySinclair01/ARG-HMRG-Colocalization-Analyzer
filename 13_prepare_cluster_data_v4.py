#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
共定位基因簇绘图数据格式化脚本 (v2.0 - 增加智能去重)

功能:
1. 读取由11号脚本生成的、已完全注释的详细共定位报告。
2. **新增：智能合并**具有双重注释 (ARG/HMRG) 的基因，避免绘图时重复。
3. 将所有共定位contigs上的基因信息，严格按照用户指定的格式进行整理。
4. 输出的文件可以直接用于绘图软件或后续的筛选脚本。

使用方法:
1. 确保已安装 pandas 库: pip install pandas
2. 将此脚本与 'colocalization_analysis_final_annotated' 文件夹放在同一目录下。
3. 在终端中运行此脚本: python 13_format_data_for_plotting.py
"""

import os
import pandas as pd
from typing import List

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
INPUT_FOLDER = "colocalization_analysis_final"
OUTPUT_FOLDER = "cluster_plotting_data_final"

# 文件后缀定义
DETAILS_SUFFIX = "_colocalization_full_details_annotated.tsv"


def process_sample(sample_name: str):
    """处理单个样本，将其详细报告转换为带智能去重的绘图格式。"""
    print(f"--- 开始处理样本: {sample_name} ---")

    details_file = os.path.join(INPUT_FOLDER, f"{sample_name}{DETAILS_SUFFIX}")

    if not os.path.exists(details_file):
        print(f"警告: 样本 {sample_name} 缺少详细报告文件，已跳过。")
        return

    try:
        df_details = pd.read_csv(details_file, sep='\t')
        
        if df_details.empty:
            print(f"信息: 样本 {sample_name} 的详细报告为空，已跳过。")
            return

        # --- *** 新增的智能去重逻辑 *** ---
        # 按protein_id分组，并将多行的信息合并到一行
        def aggregate_duplicates(group):
            # 如果只有一行，直接返回
            if len(group) == 1:
                return group.iloc[0]
            # 如果有多行（即双重注释）
            else:
                # 以第一行为基础创建新行
                aggregated_row = group.iloc[0].copy()
                # 合并 gene_type
                aggregated_row['gene_type'] = '/'.join(sorted(group['gene_type'].astype(str).unique()))
                # 合并 gene_name
                aggregated_row['gene_name'] = ' / '.join(group['gene_name'].dropna().astype(str).unique())
                return aggregated_row

        # 应用聚合函数
        # 我们只对有注释的基因进行此操作
        df_annotated = df_details[df_details['gene_type'] != 'Other'].copy()
        df_other = df_details[df_details['gene_type'] == 'Other'].copy()

        df_deduplicated_annotated = df_annotated.groupby('protein_id').apply(aggregate_duplicates).reset_index(drop=True)
        
        # 将处理过的已注释基因与未注释的基因重新合并
        df_final_deduplicated = pd.concat([df_deduplicated_annotated, df_other], ignore_index=True)
        
        print(f"  > 智能去重完成，将 {len(df_details)} 行记录合并为 {len(df_final_deduplicated)} 行。")
        # --- 结束新增逻辑 ---

        # --- 严格按照指定格式创建新的DataFrame ---
        output_df = pd.DataFrame({
            'ID': df_final_deduplicated['protein_id'],
            'source': df_final_deduplicated['contig_id'],
            'start': df_final_deduplicated['start'],
            'end': df_final_deduplicated['end'],
            'strand': df_final_deduplicated['strand'],
            'gene_type': df_final_deduplicated['gene_type'],
            'gene_name': df_final_deduplicated['gene_name']
        })
        
        # 重新排序
        output_df = output_df.sort_values(by=['source', 'start'])

        # --- 保存整理好的绘图数据 ---
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
    print("--- 整理用于基因簇可视化的数据 (带智能去重) ---")
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