#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ARG 与 HMRG 共定位分析与最终注释脚本 (版本 8.0 - 最终版)

功能:
1. 加载预先生成的、正确的HMRG注释映射文件 (bacmet_annotation_map.tsv)。
2. 读取GFF、CARD和BacMet的原始文件。
3. **正确地**解析ARG的生物学名称和HMRG的**登录号**。
4. **使用映射文件将HMRG登录号翻译为生物学基因名**。
5. 整合所有信息，生成包含**双重正确生物学名称**的最终详细报告。
"""

import os
import pandas as pd
from typing import List

# --- 用户配置区 ---
INPUT_FOLDER = "analysis_results"
OUTPUT_FOLDER = "colocalization_analysis_final" # Final output
HMRG_MAP_FILE = "bacmet_annotation_map.tsv" # The essential mapping file

CARD_SUFFIX = "_card_hits.m8"
BACMET_SUFFIX = "_bacmet_hits.tsv"
GFF_SUFFIX = "_predicted_genes.gff"


def parse_gff_attributes(attr_string: str) -> str:
    """从GFF属性中解析出基因编号。"""
    try: return attr_string.split('ID=')[1].split(';')[0].split('_')[-1]
    except: return ""

def parse_card_gene_name(sseqid: str) -> str:
    """从CARD的sseqid中提取基因名。"""
    try: return sseqid.strip().split('|')[-1]
    except: return "unknown_arg"

def parse_hmrg_accession(sseqid: str) -> str:
    """从BacMet的sseqid中提取登录号 (e.g., Q5FAM9 or WP_...)."""
    try:
        parts = str(sseqid).strip().split('|')
        db_codes = ['ref', 'gb', 'emb', 'sp', 'tr']
        
        # 检查 NCBI 格式
        for code in db_codes:
            if code in parts:
                idx = parts.index(code)
                if len(parts) > idx + 1:
                    return parts[idx + 1].split('.')[0] # 返回不带版本号的
        
        # 备用方案，处理 BacMet 内部格式 (e.g., >BAC...|abeM|tr|Q5FAM9|...)
        if len(parts) > 3:
            return parts[3]

        return "unknown_accession"
    except: return "parse_error"


def analyze_sample(sample_name: str, hmrg_map: dict):
    """处理单个样本，进行完整的注释和整合。"""
    print(f"--- 正在处理样本: {sample_name} ---")
    
    card_file = os.path.join(INPUT_FOLDER, f"{sample_name}{CARD_SUFFIX}")
    bacmet_file = os.path.join(INPUT_FOLDER, f"{sample_name}{BACMET_SUFFIX}")
    gff_file = os.path.join(INPUT_FOLDER, f"{sample_name}{GFF_SUFFIX}")

    if not all(os.path.exists(f) for f in [card_file, bacmet_file, gff_file]):
        print(f"警告: 样本 {sample_name} 缺少必需文件，已跳过。")
        return

    try:
        # 1. 读取并注释 CARD (ARG) 数据
        card_cols = ['protein_id', 'sseqid', 'pident', 'evalue', 'bitscore']
        df_card = pd.read_csv(card_file, sep='\t', header=None, usecols=[0, 1, 2, 10, 11], names=card_cols, dtype=str).sort_values('bitscore', ascending=False).drop_duplicates('protein_id')
        df_card['gene_name'] = df_card['sseqid'].apply(parse_card_gene_name)
        df_card['gene_type'] = 'ARG'
        
        # 2. 读取、注释并**翻译** BacMet (HMRG) 数据
        bacmet_cols = ['protein_id', 'sseqid', 'pident', 'evalue', 'bitscore']
        df_bacmet = pd.read_csv(bacmet_file, sep='\t', header=None, usecols=[0, 1, 2, 10, 11], names=bacmet_cols, dtype=str).sort_values('bitscore', ascending=False).drop_duplicates('protein_id')
        df_bacmet['accession'] = df_bacmet['sseqid'].apply(parse_hmrg_accession)
        
        # *** 关键翻译步骤 ***
        # 使用映射文件将登录号翻译为生物学基因名
        df_bacmet['gene_name'] = df_bacmet['accession'].map(hmrg_map).fillna(df_bacmet['accession'])
        df_bacmet['gene_type'] = 'HMRG'
        
        # 3. 合并所有注释并识别共定位contigs
        df_annotations = pd.concat([df_card, df_bacmet])
        df_annotations['contig_id'] = df_annotations['protein_id'].str.rsplit('_', n=1).str[0]
        
        contig_summary = df_annotations.groupby('contig_id')['gene_type'].unique().apply(set)
        colocalized_contigs = set(contig_summary[contig_summary.apply(lambda x: 'ARG' in x and 'HMRG' in x)].index)

        if not colocalized_contigs:
            print(f"信息: 样本 {sample_name} 中没有共定位的contigs。")
            return

        # 4. 读取GFF并与完全注释的数据合并
        gff_cols = ['contig_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df_gff = pd.read_csv(gff_file, sep='\t', comment='#', header=None, names=gff_cols)
        df_gff = df_gff[df_gff['type'] == 'CDS'].copy()
        df_gff['protein_id'] = df_gff['contig_id'] + '_' + df_gff['attributes'].apply(parse_gff_attributes)
        
        df_target_gff = df_gff[df_gff['contig_id'].isin(colocalized_contigs)]
        
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

        # 5. 保存报告
        output_file = os.path.join(OUTPUT_FOLDER, f"{sample_name}_colocalization_full_details_annotated.tsv")
        df_final.to_csv(output_file, sep='\t', index=False, float_format='%.2e')
        
        print(f"成功: 为样本 {sample_name} 生成了最终注释的详细报告。")

    except Exception as e:
        print(f"错误: 处理样本 {sample_name} 时发生严重错误: {e}")
        import traceback
        traceback.print_exc()

def main():
    """主执行函数"""
    print("=" * 60)
    print("--- 共定位分析与最终注释流程 (v8.0) ---")
    print("=" * 60)
    
    if not os.path.exists(HMRG_MAP_FILE):
        print(f"致命错误: HMRG注释文件 '{HMRG_MAP_FILE}' 未找到！")
        print("请确保 `create_hmrg_annotation_map.py` 脚本已成功运行。")
        return
        
    # 读取映射文件，将登录号作为索引
    hmrg_map = pd.read_csv(HMRG_MAP_FILE, sep='\t', index_col='accession', header=0, dtype=str).squeeze("columns").to_dict()
    print(f"成功加载 {len(hmrg_map)} 条HMRG注释。")

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    samples = sorted([f.replace(CARD_SUFFIX, "") for f in os.listdir(INPUT_FOLDER) if f.endswith(CARD_SUFFIX)])
    print(f"成功检测到 {len(samples)} 个样本: {', '.join(samples)}\n")

    for sample in samples:
        analyze_sample(sample, hmrg_map)
        print("-" * 50)

    print("所有样本分析完毕！")
    print(f"所有最终详细报告均已保存在 '{OUTPUT_FOLDER}' 文件夹中。")

if __name__ == "__main__":
    main()