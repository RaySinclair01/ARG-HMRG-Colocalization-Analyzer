#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
基因丰度综合分析脚本 (版本 7.0)

功能:
1. **保留已有功能**: 计算 "ARG-HMRG" 共定位基因对的丰度。
2. **新增功能**: 单独统计每个样本中 HMRG 的丰度。
3. 使用预先生成的注释地图文件，确保所有输出都包含正确的生物学基因名。
4. 将两类结果分别输出到不同的文件夹中。

使用方法:
1. 确保已安装 pandas 库: pip install pandas
2. 将此脚本与 'colocalization_analysis_...', 'analysis_results' 文件夹
   和 'bacmet_annotation_map.tsv' 文件放在同一目录下。
3. 在终端中运行此脚本: python 12_analyze_abundance_combined.py
"""

import os
import pandas as pd
from typing import List
from collections import Counter
import itertools

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
# 输入文件夹
COLOCALIZATION_FOLDER = "colocalization_analysis_final"
DIAMOND_FOLDER = "analysis_results"
HMRG_MAP_FILE = "bacmet_annotation_map.tsv" # HMRG "翻译字典"

# 输出文件夹
PAIR_ABUNDANCE_OUTPUT_FOLDER = "abundance_analysis_pairs" # 基因对丰度
HMRG_ABUNDANCE_OUTPUT_FOLDER = "abundance_analysis_hmrg"   # HMRG 单独丰度

# 文件后缀定义
DETAILS_SUFFIX = "_colocalization_full_details_annotated.tsv"
BACMET_SUFFIX = "_bacmet_hits.tsv"


# ==============================================================================
# --- 功能 1: 计算共定位基因对丰度 (保留的旧功能) ---
# ==============================================================================
def process_pair_abundance(sample_name: str):
    """
    处理单个样本的详细报告文件，计算共定位基因对丰度。
    """
    details_file = os.path.join(COLOCALIZATION_FOLDER, f"{sample_name}{DETAILS_SUFFIX}")
    print(f"--- 任务1: 正在分析样本 '{sample_name}' 的共定位基因对丰度 ---")

    if not os.path.exists(details_file):
        print(f"  > 信息: 未找到详细报告文件，已跳过此任务。")
        return

    try:
        df = pd.read_csv(details_file, sep='\t')
        df_annotated = df[df['gene_type'].isin(['ARG', 'HMRG'])].dropna(subset=['gene_name']).copy()

        args_per_contig = df_annotated[df_annotated['gene_type'] == 'ARG'].groupby('contig_id')['gene_name'].unique().apply(list)
        hmrgs_per_contig = df_annotated[df_annotated['gene_type'] == 'HMRG'].groupby('contig_id')['gene_name'].unique().apply(list)

        pair_counter = Counter()
        for contig_id, arg_list in args_per_contig.items():
            if contig_id in hmrgs_per_contig.index:
                hmrg_list = hmrgs_per_contig[contig_id]
                all_pairs_on_contig = set(itertools.product(arg_list, hmrg_list))
                pair_counter.update(all_pairs_on_contig)

        if not pair_counter:
            print(f"  > 信息: 样本 {sample_name} 中没有有效的共现基因对。")
            return

        abundance_data = [{'arg_gene_name': p[0], 'hmrg_gene_name': p[1], 'shared_contig_count': c} for p, c in pair_counter.items()]
        abundance_df = pd.DataFrame(abundance_data).sort_values(by='shared_contig_count', ascending=False)

        output_file = os.path.join(PAIR_ABUNDANCE_OUTPUT_FOLDER, f"{sample_name}_pair_abundance.tsv")
        abundance_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"  > 成功: 基因对丰度报告已保存至: {output_file}")
    except Exception as e:
        print(f"  > 错误: 处理基因对丰度时发生错误: {e}")


# ==============================================================================
# --- 功能 2: 计算HMRG丰度 (新增功能) ---
# ==============================================================================
def parse_hmrg_accession(sseqid: str) -> str:
    """从BacMet的sseqid中提取登录号 (e.g., WP_...)"""
    try:
        parts = str(sseqid).strip().split('|')
        db_codes = ['ref', 'gb', 'emb', 'sp', 'tr']
        for code in db_codes:
            if code in parts: return parts[parts.index(code) + 1].split('.')[0]
        if len(parts) > 3: return parts[3]
        return "unknown_accession"
    except: return "parse_error"

def process_hmrg_abundance(sample_name: str, hmrg_map: dict):
    """
    处理单个样本的原始比对文件，计算HMRG丰度。
    """
    bacmet_file = os.path.join(DIAMOND_FOLDER, f"{sample_name}{BACMET_SUFFIX}")
    print(f"--- 任务2: 正在分析样本 '{sample_name}' 的HMRG丰度 ---")

    if not os.path.exists(bacmet_file):
        print(f"  > 信息: 未找到HMRG比对文件，已跳过此任务。")
        return

    try:
        # 1. 读取原始HMRG比对结果
        bacmet_cols = ['protein_id', 'sseqid']
        df_bacmet = pd.read_csv(bacmet_file, sep='\t', header=None, usecols=[0, 1], names=bacmet_cols, dtype=str)
        
        # 2. 解析登录号并翻译为生物学基因名
        df_bacmet['accession'] = df_bacmet['sseqid'].apply(parse_hmrg_accession)
        df_bacmet['gene_name'] = df_bacmet['accession'].map(hmrg_map).fillna(df_bacmet['accession'])

        # 3. 统计每个基因名的出现次数 (丰度)
        abundance_counts = df_bacmet['gene_name'].value_counts()
        
        # 4. 转换为DataFrame并排序
        abundance_df = abundance_counts.reset_index()
        abundance_df.columns = ['hmrg_gene_name', 'abundance_count']
        
        # 5. 保存丰度报告
        output_file = os.path.join(HMRG_ABUNDANCE_OUTPUT_FOLDER, f"{sample_name}_hmrg_abundance.tsv")
        abundance_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"  > 成功: HMRG丰度报告已保存至: {output_file}")
    except Exception as e:
        print(f"  > 错误: 处理HMRG丰度时发生错误: {e}")


# ==============================================================================
# --- 主执行流程 ---
# ==============================================================================
def main():
    """主执行函数"""
    print("=" * 60)
    print("--- 基因丰度综合分析流程 ---")
    print("=" * 60)

    # 检查必需的文件夹和文件
    if not os.path.isdir(COLOCALIZATION_FOLDER) or not os.path.isdir(DIAMOND_FOLDER):
        print("错误: 必需的输入文件夹 'colocalization_analysis_...' 或 'analysis_results' 不存在！")
        return
    if not os.path.exists(HMRG_MAP_FILE):
        print(f"致命错误: HMRG注释文件 '{HMRG_MAP_FILE}' 未找到！")
        return
        
    # 加载HMRG注释地图
    hmrg_map = pd.read_csv(HMRG_MAP_FILE, sep='\t', index_col='accession', header=0, dtype=str).squeeze("columns").to_dict()
    print(f"成功加载 {len(hmrg_map)} 条HMRG注释。\n")

    # 创建输出文件夹
    os.makedirs(PAIR_ABUNDANCE_OUTPUT_FOLDER, exist_ok=True)
    os.makedirs(HMRG_ABUNDANCE_OUTPUT_FOLDER, exist_ok=True)

    # 自动识别所有样本
    samples: List[str] = sorted(list(set(
        f.replace(DETAILS_SUFFIX, "") for f in os.listdir(COLOCALIZATION_FOLDER) if f.endswith(DETAILS_SUFFIX)
    )))
    
    print(f"成功检测到 {len(samples)} 个样本: {', '.join(samples)}\n")

    # 为每个样本执行两个任务
    for sample in samples:
        process_pair_abundance(sample)
        process_hmrg_abundance(sample, hmrg_map)
        print("-" * 60)

    print("所有样本分析完毕！")
    print(f"共定位基因对丰度报告保存在: '{PAIR_ABUNDANCE_OUTPUT_FOLDER}'")
    print(f"HMRG单独丰度报告保存在: '{HMRG_ABUNDANCE_OUTPUT_FOLDER}'")
    print("=" * 60)

if __name__ == "__main__":
    main()