#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文件名: homology_clustering_TPMC-S.py
# 描述: 基于物理位置聚类的结果，对 TPMC-S 防御基因进行序列同源性聚类 (CD-HIT)

import pandas as pd
import os
import subprocess
from Bio import SeqIO
import sys
import warnings
import multiprocessing as mp
from Bio import BiopythonDeprecationWarning

# 忽略 Biopython 的注释解析警告
warnings.simplefilter('ignore', BiopythonDeprecationWarning)

def main():
    print("开始 TPMC-S 防御基因同源聚类分析 (Algorithm: CD-HIT, Threshold: 40%)...")
    
    # === 路径配置 ===
    base_dir = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang"
    
    # 依赖上一步(物理聚类)生成的真实坐标 csv 文件
    input_csv = os.path.join(base_dir, "02_physical_clustering/TPMC-S_20260420/real_position_clusters_with_context_TPMC-S.csv")
    
    # FASTA文件所在总目录 (保留原始路径)
    genome_dir = "/data5/zhanghaohong/projects/TPMC_BGC/prodigal/TPMC-S"
    
    # 输出目录 (如果是首次运行会自己创建)
    output_dir = os.path.join(base_dir, "03_homology_clustering/TPMC-S_20260420")
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(input_csv):
        print(f"错误: 找不到输入文件 {input_csv}")
        return

    # 读取CSV文件
    print(f"读取物理聚类结果: {input_csv}")
    df = pd.read_csv(input_csv)
    threads = 64
    # 1. 提取基因序列
    all_proteins_file = extract_sequences(df, genome_dir, output_dir, threads=threads)
    
    # 2. 运行 CD-HIT 聚类 (设为 40% 相似度, 即 0.4)
    if all_proteins_file:
        cdhit_clstr_file = run_cdhit(all_proteins_file, output_dir, similarity=0.4, threads=threads)
        
        # 3. 分析聚类结果 (解析 .clstr 文件)
        clusters = analyze_clusters(cdhit_clstr_file, df, output_dir)
        
        # 4. 提取代表序列
        extract_representatives(clusters, all_proteins_file, output_dir)
        
        print("\n所有分析步骤圆满完成!")

def _process_fasta_worker(args):
    """多进程执行的辅助函数：解析单个FASTA文件并返回匹配的序列字符串列表"""
    prefix, faa_path, required_genes = args
    matched_seqs = []
    try:
        # 直接读取文件
        for record in SeqIO.parse(faa_path, "fasta-pearson"):
            if record.id in required_genes:
                # 转换成 FASTA 字符串返回，避免在进程间传输复杂的 SeqRecord 对象
                fasta_str = f">{record.id} source_prefix={prefix} mapped_id={record.id}\n{str(record.seq)}\n"
                matched_seqs.append(fasta_str)
    except Exception as e:
        print(f"读取 {faa_path} 出错: {e}")
    return matched_seqs

def extract_sequences(df, genome_dir, output_dir, threads=64):
    print("正在搜寻并提取簇内蛋白质序列...")
    all_proteins_file = os.path.join(output_dir, "TPMC-S_cluster_proteins.faa")
    
    # 获取唯一的基因列表
    required_genes = set(df['gene_name'].astype(str))
    print(f"共有 {len(required_genes)} 个靶标基因需要提取。")
    
    # [修复核心]: 从 Contig ID 中提取出原始的 Bin 前缀
    file_prefixes = set(df['genome_id'].astype(str).apply(lambda x: x.split('_k141_')[0]).unique())
    
    # 扫描目录寻找对应的 FASTA 文件
    protein_files = {}
    for root, dirs, files in os.walk(genome_dir):
        for file in files:
            # 增加常见的蛋白序列后缀
            if file.endswith(('.fasta', '.faa', '.fa', '.pep')):
                file_path = os.path.join(root, file)
                for prefix in file_prefixes:
                    # 改为子串包含匹配，以防止中间有其他后缀干扰
                    if prefix in file:
                        protein_files[prefix] = file_path
                        break
                        
    print(f"在目录中定位到 {len(protein_files)} 个相关的源 FASTA 文件。")
    
    found_count = 0
    
    # 准备多进程参数
    worker_args = [(prefix, faa_path, required_genes) for prefix, faa_path in protein_files.items()]
    print(f"启动 {threads} 个进程进行并行提取...")
    with open(all_proteins_file, 'w') as out_f:
        with mp.Pool(processes=threads) as pool:
            # imap_unordered 会在任务完成时立刻生成结果
            for matched_seqs in pool.imap_unordered(_process_fasta_worker, worker_args):
                for seq_str in matched_seqs:
                    out_f.write(seq_str)
                    found_count += 1
                
    print(f"成功提取了 {found_count} 个蛋白质序列，保存在 {all_proteins_file}")
    
    if found_count == 0:
        print("警告：未能匹配提取到任何序列，请检查文件名前缀或路径！")
        return None
    return all_proteins_file

def run_cdhit(input_fasta, output_dir, similarity=0.4, threads=32):
    print(f"\n运行CD-HIT聚类 (相似度阈值: {similarity})...")
    
    if os.path.getsize(input_fasta) == 0:
        raise ValueError("输入序列文件为空！")
        
    output_prefix = os.path.join(output_dir, "cdhit_clusters")
    
    # 词长选择建议：0.7 -> 5, 0.6 -> 4, <0.5 -> 2
    word_len = 3 if similarity >= 0.5 else 2
    
    # 大多数 Conda 环境 (如您的 torch 环境) 安装 CD-HIT 后，系统命令直接为 cd-hit
    cdhit_exe = "/home/xiongjiayun/data5/anti_virus/Rfseq_defense_system/find_newsystem/cdhit/cd-hit"

    cmd = [
        cdhit_exe,
        "-i", input_fasta,
        "-o", output_prefix,
        "-c", str(similarity),  # 相似度阈值
        "-n", str(word_len),    # 词长
        "-d", "0",              # 描述行长度 (0=完整不截断)
        "-T", str(threads),     # 线程数
        "-M", "0"               # 内存限制 (0=无限制)
    ]
    
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        print(f"\nCRITICAL ERROR: 找不到系统命令 '{cdhit_exe}'")
        print("请在您的 conda 环境中安装 CD-HIT，例如运行命令：")
        print("conda install -c bioconda cd-hit")
        sys.exit(1)
    
    clstr_file = output_prefix + ".clstr"
    if not os.path.exists(clstr_file):
        raise FileNotFoundError(f"CD-HIT 未能生成结果文件: {clstr_file}")
        
    print(f"CD-HIT聚类完成，结果保存在 {clstr_file}")
    return clstr_file

def analyze_clusters(clstr_file, original_df, output_dir):
    print("\n解析 CD-HIT 聚类结果...")
    
    clusters = {}
    cluster_rows = []
    
    current_cluster_id = 0
    current_genes = []
    
    def save_cluster(cid, genes_list):
        if not genes_list: return
        real_cid = cid + 1 # 设为 1-based
        clusters[real_cid] = genes_list
        for gene in genes_list:
            cluster_rows.append({'homology_cluster_id': real_cid, 'gene_name': gene})

    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                if current_genes:
                    save_cluster(current_cluster_id, current_genes)
                    current_cluster_id += 1
                current_genes = []
            else:
                try:
                    parts = line.split('>')
                    if len(parts) > 1:
                        content = parts[1]
                        # 提取并丢弃可能存在的附加描述字段(source_prefix等)
                        gene_id_full = content.split('...')[0].replace('*', '').strip()
                        gene_id = gene_id_full.split()[0]
                        current_genes.append(gene_id)
                except Exception:
                    continue
                    
        if current_genes:
            save_cluster(current_cluster_id, current_genes)


    cluster_df = pd.DataFrame(cluster_rows)
    print(f"共识别出 {len(clusters)} 个独立同源家族 (CD-HIT 簇)。")
    
    # 合并原始组分信息
    meta_df = original_df[['gene_name', 'gene_type', 'genome_id', 'cluster_id']].drop_duplicates() 
    meta_df = meta_df.rename(columns={'cluster_id': 'physical_cluster_id'})
    merged_df = pd.merge(cluster_df, meta_df, on='gene_name', how='left')
    
    # 保存详细映射表
    detail_path = os.path.join(output_dir, "TPMC-S_homology_clusters_detail.csv")
    merged_df.to_csv(detail_path, index=False)
    
    # 分类统计逻辑
    stats = merged_df.groupby(['homology_cluster_id', 'gene_type']).size().unstack(fill_value=0)
    if 'candidate' not in stats.columns: stats['candidate'] = 0
    if 'known_defense' not in stats.columns: stats['known_defense'] = 0
    
    stats['total_genes'] = stats.sum(axis=1)
    stats = stats.sort_values('total_genes', ascending=False).reset_index()
    
    def classify_homology(row):
        if row['known_defense'] == 0 and row['candidate'] > 0:
            return "Potential_Novel_Family"
        elif row['known_defense'] > 0 and row['candidate'] > 0:
            return "Expanded_Known_Family"
        else:
            return "Known_Family_Only"
            
    stats['homology_type'] = stats.apply(classify_homology, axis=1)
    
    stats_path = os.path.join(output_dir, "TPMC-S_homology_cluster_statistics.csv")
    stats.to_csv(stats_path, index=False)
    
    print("\n同源簇(家族)统计排名前10:")
    print(stats.head(10).to_string(index=False))
    
    return clusters

def extract_representatives(clusters, all_proteins_file, output_dir):
    print("\n提取代表序列...")
    seq_dict = SeqIO.to_dict(SeqIO.parse(all_proteins_file, "fasta-pearson"))
    
    representative_seqs = {}
    for cid, genes in clusters.items():
        if genes:
            representative_seqs[cid] = genes[0] 
            
    rep_file_path = os.path.join(output_dir, "TPMC-S_homology_cluster_representatives.faa")
    
    count = 0
    with open(rep_file_path, 'w') as out_f:
        for cid in sorted(clusters.keys()):
            gene_id = representative_seqs.get(cid)
            if gene_id and gene_id in seq_dict:
                record = seq_dict[gene_id]
                record.id = f"HomologyCluster_{cid}|{record.id}"
                SeqIO.write(record, out_f, "fasta")
                count += 1
                
    print(f"已成功提取 {count} 条家族代表序列至 {rep_file_path}")

if __name__ == "__main__":
    main()