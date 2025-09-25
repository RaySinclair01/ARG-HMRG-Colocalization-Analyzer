#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
高丰度共定位基因簇数据筛选脚本 (按Contig密度筛选 - 增加智能去重功能)

功能:
1. 读取完整的共定位详细报告。
2. 计算每个contig的“共定位密度分数”。
3. 筛选出排名最高的N个contigs。
4. **新增：智能合并双重注释** (ARG/HMRG) 的基因，避免重复行。
5. 将最终的、干净的数据输出，用于可视化。
"""

import os
import pandas as pd
from typing import List

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
DETAILS_FOLDER = "colocalization_analysis_final"
OUTPUT_FOLDER = "top_10_contigs_for_plotting"

# 文件后缀定义
DETAILS_SUFFIX = "_colocalization_full_details_annotated.tsv"

# *** 关键控制参数 ***
NUM_TOP_CONTIGS_TO_PLOT = 10


def process_sample(sample_name: str):
    """处理单个样本，筛选出密度最高的contigs并进行智能去重。"""
    print(f"--- 开始处理样本: {sample_name} ---")

    details_file = os.path.join(DETAILS_FOLDER, f"{sample_name}{DETAILS_SUFFIX}")

    if not os.path.exists(details_file):
        print(f"警告: 样本 {sample_name} 缺少详细报告文件，已跳过。")
        return

    try:
        df_details = pd.read_csv(details_file, sep='\t')
        if df_details.empty:
            print(f"信息: 样本 {sample_name} 的详细报告为空，已跳过。")
            return

        # --- 2. 计算每个contig的共定位密度分数 ---
        gene_counts = df_details.groupby(['contig_id', 'gene_type']).size().unstack(fill_value=0)
        
        if 'ARG' not in gene_counts.columns: gene_counts['ARG'] = 0
        if 'HMRG' not in gene_counts.columns: gene_counts['HMRG'] = 0
        
        gene_counts['density_score'] = gene_counts['ARG'] * gene_counts['HMRG']
        top_contigs_df = gene_counts.sort_values(by='density_score', ascending=False)
        top_contigs_df = top_contigs_df[top_contigs_df['density_score'] > 0]
        
        if top_contigs_df.empty:
            print(f"信息: 样本 {sample_name} 中没有找到包含ARG和HMRG的contigs。")
            return

        target_contigs = set(top_contigs_df.head(NUM_TOP_CONTIGS_TO_PLOT).index)
        
        print(f"根据共定位密度，筛选出Top {len(target_contigs)} 个contigs进行绘图。")

        # --- 4. 筛选并进行智能去重 ---
        df_plot_data = df_details[df_details['contig_id'].isin(target_contigs)].copy()

        # *** 新增的智能去重逻辑 ***
        # 按protein_id分组，并将多行的信息合并到一行
        def aggregate_duplicates(group):
            # 如果只有一行，直接返回
            if len(group) == 1:
                return group
            # 如果有多行（即双重注释）
            else:
                # 以第一行为基础
                first_row = group.iloc[0].copy()
                # 合并 gene_type 和 gene_name
                first_row['gene_type'] = '/'.join(sorted(group['gene_type'].unique()))
                first_row['gene_name'] = ' / '.join(group['gene_name'].dropna().astype(str).unique())
                # 返回合并后的单行
                return first_row.to_frame().T

        # 应用聚合函数
        df_deduplicated = df_plot_data.groupby('protein_id').apply(aggregate_duplicates).reset_index(drop=True)
        
        print(f"  > 智能去重完成，将 {len(df_plot_data)} 行记录合并为 {len(df_deduplicated)} 行。")
        # *** 结束新增逻辑 ***

        
        # 按照您指定的格式创建新的DataFrame
        output_df = pd.DataFrame({
            'ID': df_deduplicated['protein_id'],
            'source': df_deduplicated['contig_id'],
            'start': df_deduplicated['start'],
            'end': df_deduplicated['end'],
            'strand': df_deduplicated['strand'],
            'gene_type': df_deduplicated['gene_type'],
            'gene_name': df_deduplicated['gene_name']
        })
        
        output_df = output_df.sort_values(by=['source', 'start'])

        # --- 5. 保存整理好的绘图数据 ---
        output_file = os.path.join(OUTPUT_FOLDER, f"{sample_name}_top_{NUM_TOP_CONTIGS_TO_PLOT}_contigs_plot_data.tsv")
        output_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"成功: 已将用于绘图的Top {len(target_contigs)} contigs数据保存至: {output_file}")

    except Exception as e:
        print(f"错误: 处理样本 {sample_name} 时发生严重错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主执行函数"""
    print("=" * 60)
    print(f"--- 筛选Top {NUM_TOP_CONTIGS_TO_PLOT} 个高密度共定位Contig (带智能去重) ---")
    print("=" * 60)

    if not os.path.isdir(DETAILS_FOLDER):
        print(f"错误: 必需的输入文件夹 '{DETAILS_FOLDER}' 不存在！")
        return

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    samples: List[str] = sorted(
        [f.replace(DETAILS_SUFFIX, "") for f in os.listdir(DETAILS_FOLDER) if f.endswith(DETAILS_SUFFIX)]
    )
    
    if not samples:
        print(f"错误: 在 '{DETAILS_FOLDER}' 中没有找到任何详细报告文件。")
        return
        
    print(f"成功检测到 {len(samples)} 个样本: {', '.join(samples)}\n")

    for sample in samples:
        process_sample(sample)
        print("-" * 50)

    print("所有样本筛选完毕！")
    print(f"所有用于绘图的精华数据文件均已保存在 '{OUTPUT_FOLDER}' 文件夹中。")
    print("=" * 60)


if __name__ == "__main__":
    main()