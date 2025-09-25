#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HMRG 注释地图创建脚本 (v5.0 - 最终版)

功能:
1. 读取包含多种复杂格式头信息的 BacMet FASTA 文件。
2. **使用全新的、能处理所有已知格式的智能解析器**。
3. 按以下优先级提取基因名:
   - 优先提取末尾方括号 [...] 内的名称 (如 'acrF/envD')。
   - 其次提取 GN= 字段后的名称 (如 'abeM')。
   - 最后使用备用方案。
4. 创建一个完美的、可用于后续所有分析的“翻译字典”。
"""

import os
import re

# ==============================================================================
# --- 用户配置区 ---
# ==============================================================================
BACMET_FASTA_PATH = "BacMet_combined.fasta"
OUTPUT_MAP_FILE = "bacmet_annotation_map.tsv"
# --- 结束配置 ---


def process_header(header: str, bracket_pattern: re.Pattern, gn_pattern: re.Pattern) -> tuple:
    """
    智能处理单条FASTA头信息，返回 (登录号, 基因名)
    """
    accession = None
    gene_name = None
    header = header.strip().strip('>')
    
    # --- 步骤 1: 提取基因名 ---
    
    # 策略 1 (最高优先级): 寻找末尾的方括号 [...]
    bracket_match = bracket_pattern.search(header)
    if bracket_match:
        gene_name = bracket_match.group(1)
    
    # 策略 2 (第二优先级): 寻找 GN= 字段
    elif gn_pattern:
        gn_match = gn_pattern.search(header)
        if gn_match:
            gene_name = gn_match.group(1)

    # --- 步骤 2: 提取登录号 ---
    parts = header.split('|')
    
    # 处理 NCBI 格式 (e.g., >gi|...|ref|WP_...|)
    db_codes = ['ref', 'gb', 'emb', 'sp', 'tr']
    for code in db_codes:
        if code in parts:
            try:
                idx = parts.index(code)
                if len(parts) > idx + 1:
                    accession = parts[idx + 1]
                    break
            except ValueError:
                continue
    
    # 如果不是 NCBI 格式，则处理 BacMet 内部格式 (e.g., >BAC0001|abeM|...)
    if not accession and len(parts) > 3:
        accession = parts[3] # 第四个字段是登录号，如 Q5FAM9
        # 如果此时基因名仍然未知，使用第二个字段作为备用
        if not gene_name:
            gene_name = parts[1]

    # --- 步骤 3: 最终检查 ---
    # 如果基因名仍然未知，但登录号已知，则用登录号（不带版本号）作为最终的备用名
    if not gene_name and accession:
        gene_name = accession.split('.')[0]
        
    return accession, gene_name


def main():
    """主执行函数"""
    print("--- 正在创建 HMRG 注释地图 (v5.0 Final) ---")
    print(f"读取 FASTA 文件: {BACMET_FASTA_PATH}")

    if not os.path.exists(BACMET_FASTA_PATH):
        print(f"\n致命错误: FASTA 文件未在 '{BACMET_FASTA_PATH}' 找到！")
        return

    annotation_map = {}
    # 编译正则表达式以提高效率
    bracket_pattern = re.compile(r"\[([^\]]+)\]\s*$") # 匹配末尾的 [...]
    gn_pattern = re.compile(r"GN=([\w\(\)\-\_/\.]+)") # 匹配 GN=...

    with open(BACMET_FASTA_PATH, "r", encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith(">"):
                accession, gene_name = process_header(line, bracket_pattern, gn_pattern)
                
                if accession and gene_name:
                    # 我们需要两个版本的登录号作为键：带版本号的和不带的
                    annotation_map[accession] = gene_name
                    annotation_map[accession.split('.')[0]] = gene_name

    with open(OUTPUT_MAP_FILE, "w", encoding='utf-8') as out_f:
        out_f.write("accession\tgene_name\n")
        # 只写入去重后的映射关系，以不带版本号的登录号为准
        unique_map = {acc.split('.')[0]: name for acc, name in annotation_map.items()}
        for acc, name in unique_map.items():
            out_f.write(f"{acc}\t{name}\n")

    print(f"\n成功！已创建包含 {len(unique_map)} 条记录的注释地图。")
    print(f"文件已保存为: {OUTPUT_MAP_FILE}")

if __name__ == "__main__":
    main()