#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ARG 与 HMRG 共定位分析与数据整合脚本 (版本 6.0 - 最终版)

功能:
1. 整合 GFF (位置), CARD (ARG注释), BacMet (HMRG注释) 三大来源的数据。
2. **精确地**根据用户指定的最终规则从sseqid中解析出ARG和HMRG的基因名。
   - ARG 名称: 取'|'分隔的最后一部分 (例如: ceoB)。
   - HMRG 名称: 取'|'分隔的第三部分 (例如: ref, gb, emb)。
3. 识别共定位contigs并生成一个包含所有正确信息的终极详细报告。

使用方法:
1. 确保已安装 pandas 库: pip install pandas
2. 将此脚本与 'analysis_results' 文件夹放在同一目录下。
3. 在终端中运行此脚本: python 11_run_colocalization_analysis_final.py
"""

import os
import pandas as pd
from typing import List

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
INPUT_FOLDER = "analysis_results"
OUTPUT_FOLDER = "colocalization_analysis_final_v3" # 使用新文件夹以避免混淆

# 文件后缀定义
CARD_SUFFIX = "_card_hits.m8"
BACMET_SUFFIX = "_bacmet_hits.tsv"
GFF_SUFFIX = "_predicted_genes.gff"


def parse_gff_attributes(attr_string: str) -> str:
    """从GFF的attributes列中解析出蛋白质ID的数字部分。"""
    for item in attr_string.split(';'):
        if item.startswith('ID='):
            return item.split('_')[-1]
    return ""

def parse_card_gene_name(sseqid: str) -> str:
    """
    从CARD比对结果的sseqid中解析出基因名称 (例如: ceoB)。
    规则：取'|'分隔的最后一部分。
    """
    try:
        return sseqid.strip().split('|')[-1]
    except Exception:
        return "card_parse_error"

def parse_bacmet_gene_name(sseqid: str) -> str:
    """
    *** 最终修正版 HMRG 解析函数 ***
    从BacMet比对结果的sseqid中解析出基因来源类型 (例如: ref, gb, emb)。
    规则：取'|'分隔的第三部分。
    """
    try:
        parts = sseqid.strip().split('|')
        # 基因来源类型是第三个字段 (索引为2)
        if len(parts) > 2:
            return parts[2]
        return "unknown_format"
    except Exception:
        return "bacmet_parse_error"


def analyze_sample(sample_name: str):
    """处理单个样本，整合所有数据源。"""
    print(f"--- 开始处理样本: {sample_name} ---")

    # 定义所有文件的路径
    card_file = os.path.join(INPUT_FOLDER, f"{sample_name}{CARD_SUFFIX}")
    bacmet_file = os.path.join(INPUT_FOLDER, f"{sample_name}{BACMET_SUFFIX}")
    gff_file = os.path.join(INPUT_FOLDER, f"{sample_name}{GFF_SUFFIX}")

    # 检查文件是否存在
    for f in [card_file, bacmet_file, gff_file]:
        if not os.path.exists(f):
            print(f"警告: 缺少文件 {f}，跳过样本 {sample_name}。")
            return

    try:
        # --- 1. 读取并处理CARD (ARG) 比对结果 ---
        card_cols = ['protein_id', 'sseqid', 'pident', 'evalue', 'bitscore']
        df_card = pd.read_csv(card_file, sep='\t', header=None, usecols=[0, 1, 2, 10, 11], names=card_cols, dtype=str)
        df_card['gene_name'] = df_card['sseqid'].apply(parse_card_gene_name)
        df_card['gene_type'] = 'ARG'
        df_card = df_card.sort_values('bitscore', ascending=False).drop_duplicates('protein_id')

        # --- 2. 读取并处理BacMet (HMRG) 比对结果 ---
        bacmet_cols = ['protein_id', 'sseqid', 'pident', 'evalue', 'bitscore']
        df_bacmet = pd.read_csv(bacmet_file, sep='\t', header=None, usecols=[0, 1, 2, 10, 11], names=bacmet_cols, dtype=str)
        df_bacmet['gene_name'] = df_bacmet['sseqid'].apply(parse_bacmet_gene_name) # 使用最终修正的HMRG解析函数
        df_bacmet['gene_type'] = 'HMRG'
        df_bacmet = df_bacmet.sort_values('bitscore', ascending=False).drop_duplicates('protein_id')

        # 合并所有注释
        df_annotations = pd.concat([df_card, df_bacmet])
        df_annotations['contig_id'] = df_annotations['protein_id'].str.rsplit('_', n=1).str[0]

        # --- 3. 找出共定位的Contigs ---
        contig_summary = df_annotations.groupby('contig_id')['gene_type'].unique().apply(set)
        colocalized_contigs = set(contig_summary[contig_summary.apply(lambda x: 'ARG' in x and 'HMRG' in x)].index)

        if not colocalized_contigs:
            print(f"信息: 样本 {sample_name} 中未发现共定位的contigs。")
            return
        
        print(f"找到 {len(colocalized_contigs)} 个共定位的contigs。正在整合详细信息...")

        # --- 4. 读取GFF文件，并与注释信息合并 ---
        gff_cols = ['contig_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df_gff = pd.read_csv(gff_file, sep='\t', comment='#', header=None, names=gff_cols)
        df_gff = df_gff[df_gff['type'] == 'CDS'].copy()
        
        df_gff['protein_id'] = df_gff['contig_id'] + '_' + df_gff['attributes'].apply(parse_gff_attributes)
        
        df_target_gff = df_gff[df_gff['contig_id'].isin(colocalized_contigs)].copy()
        
        # 合并注释信息
        df_final = pd.merge(df_target_gff, df_annotations, on='protein_id', how='left')

        df_final['gene_type'] = df_final['gene_type'].fillna('Other')
        df_final = df_final.sort_values(by=['contig_id_x', 'start'])
        
        output_cols = {
            'contig_id_x': 'contig_id', 'protein_id': 'protein_id', 'start': 'start',
            'end': 'end', 'strand': 'strand', 'gene_type': 'gene_type',
            'gene_name': 'gene_name', 'pident': 'pident', 'evalue': 'evalue',
            'bitscore': 'bitscore'
        }
        df_final = df_final[list(output_cols.keys())].rename(columns=output_cols)

        # --- 5. 保存最终报告 ---
        output_file = os.path.join(OUTPUT_FOLDER, f"{sample_name}_colocalization_full_details.tsv")
        df_final.to_csv(output_file, sep='\t', index=False, float_format='%.2e')
        
        print(f"成功: 最终详细报告已保存至: {output_file}")

    except Exception as e:
        print(f"错误: 处理样本 {sample_name} 时发生严重错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主执行函数"""
    print("=" * 60)
    print("--- ARG/HMRG 共定位数据整合流程 (v6.0 Definitive) ---")
    print("=" * 60)

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    
    samples: List[str] = sorted(list(set(
        f.replace(CARD_SUFFIX, "") for f in os.listdir(INPUT_FOLDER) if f.endswith(CARD_SUFFIX)
    )))
    
    if not samples:
        print(f"错误: 在 '{INPUT_FOLDER}' 中没有找到任何样本的比对结果文件。")
        return
        
    print(f"成功检测到 {len(samples)} 个样本: {', '.join(samples)}\n")

    for sample in samples:
        analyze_sample(sample)
        print("-" * 50)

    print("所有样本处理完毕！")
    print(f"所有详细报告均已保存在 '{OUTPUT_FOLDER}' 文件夹中。")
    print("=" * 60)


if __name__ == "__main__":
    main()