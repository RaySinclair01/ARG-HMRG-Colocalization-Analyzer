#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
共定位基因对丰度分析脚本 (版本 6.0 - 最终版)

功能:
1. 读取由最终版11号脚本生成的、已完全注释的详细共定位报告。
2. 直接利用报告中已有的基因名进行分析，不再进行任何重新解析。
3. 计算每个 "ARG名称 - HMRG名称" 配对共同出现的 unique contig 数量 (丰度)。
4. 为每个样本生成一份按丰度排序的、清晰的共现基因对报告。

使用方法:
1. 确保已安装 pandas 库: pip install pandas
2. 将此脚本与 'colocalization_analysis_final_v3' 文件夹放在同一目录下。
3. 在终端中运行此脚本: python 12_analyze_abundance_final.py
"""

import os
import pandas as pd
from typing import List
from collections import Counter
import itertools

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
INPUT_FOLDER = "colocalization_analysis_final_v3"
OUTPUT_FOLDER = "abundance_analysis_final"
DETAILS_SUFFIX = "_colocalization_full_details.tsv"


def process_details_file(sample_name: str):
    """
    处理单个样本的详细报告文件，计算基因对丰度。
    """
    details_file = os.path.join(INPUT_FOLDER, f"{sample_name}{DETAILS_SUFFIX}")

    print(f"--- 正在分析样本: {sample_name} ---")

    if not os.path.exists(details_file):
        print(f"信息: 未找到样本 {sample_name} 的详细报告文件，已跳过。")
        return

    try:
        df = pd.read_csv(details_file, sep='\t')
        
        # 筛选出有明确基因名和类型的行 (ARG或HMRG)
        df_annotated = df[df['gene_type'].isin(['ARG', 'HMRG'])].dropna(subset=['gene_name']).copy()

        # --- 核心逻辑：计算共现丰度 ---
        # 1. 按contig分组，找出每个contig上的所有ARG和HMRG名称
        args_per_contig = df_annotated[df_annotated['gene_type'] == 'ARG'].groupby('contig_id')['gene_name'].unique().apply(list)
        hmrgs_per_contig = df_annotated[df_annotated['gene_type'] == 'HMRG'].groupby('contig_id')['gene_name'].unique().apply(list)

        # 2. 对每个contig，创建所有可能的 ARG-HMRG 配对
        pair_counter = Counter()
        for contig_id, arg_list in args_per_contig.items():
            if contig_id in hmrgs_per_contig.index:
                hmrg_list = hmrgs_per_contig[contig_id]
                # 使用 itertools.product 创建笛卡尔积，得到所有配对
                all_pairs_on_contig = set(itertools.product(arg_list, hmrg_list))
                # 对每个独特的配对，计数器加1
                pair_counter.update(all_pairs_on_contig)

        if not pair_counter:
            print(f"信息: 样本 {sample_name} 中没有有效的共现基因对。")
            return

        # 3. 将计数结果转换为DataFrame
        abundance_data = [
            {'arg_gene_name': pair[0], 'hmrg_gene_name': pair[1], 'shared_contig_count': count}
            for pair, count in pair_counter.items()
        ]
        abundance_df = pd.DataFrame(abundance_data)
        
        # 按丰度（共享contig数量）降序排序
        abundance_df = abundance_df.sort_values(by='shared_contig_count', ascending=False)

        # --- 4. 保存丰度报告 ---
        output_file = os.path.join(OUTPUT_FOLDER, f"{sample_name}_abundance_report.tsv")
        abundance_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"成功: 为样本 {sample_name} 生成了最终丰度报告。")
        print(f"报告已保存至: {output_file}")
        print("丰度最高的Top 5基因对:")
        print(abundance_df.head(5).to_string(index=False))

    except Exception as e:
        print(f"错误: 处理样本 {sample_name} 时发生严重错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主执行函数"""
    print("=" * 60)
    print("--- 共定位基因对丰度分析 (最终版) ---")
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
        
    print(f"成功检测到 {len(samples)} 个样本的详细报告: {', '.join(samples)}\n")

    for sample in samples:
        process_details_file(sample)
        print("-" * 60)

    print("所有样本分析完毕！")
    print(f"所有最终丰度报告均已保存在 '{OUTPUT_FOLDER}' 文件夹中。")
    print("=" * 60)


if __name__ == "__main__":
    main()