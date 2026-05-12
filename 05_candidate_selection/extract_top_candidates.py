import pandas as pd
import os
import ast
import re

# ================= 配置区域 =================
BASE_DIR = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang"
MIN_OCCURRENCE = 2  # 最低跨基因组重现次数阈值

def parse_fasta(fasta_path):
    """读取 FASTA 文件，返回 {header: sequence} 字典"""
    sequences = {}
    current_header = None
    current_seq = []
    
    if not os.path.exists(fasta_path):
        return {}

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                # 只取第一个空格前的ID，去除 >
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = "".join(current_seq)
    return sequences

def process_dataset(env_name):
    print(f"\n{'='*50}")
    print(f"🚀 开始处理数据集: {env_name}")
    print(f"{'='*50}")
    
    # 动态路径配置
    arch_freq_file = os.path.join(BASE_DIR, f"04_system_analysis/{env_name}_20260420/{env_name}_system_architecture_frequency.csv")
    homology_detail_file = os.path.join(BASE_DIR, f"03_homology_clustering/{env_name}_20260420/{env_name}_homology_clusters_detail.csv")
    # 原蛋白序列文件 (这里假设您从聚类前的全体蛋白库提取，因为代表序列文件名可能不同)
    all_proteins_file = os.path.join(BASE_DIR, f"03_homology_clustering/{env_name}_20260420/{env_name}_cluster_proteins.faa")
    
    output_dir = os.path.join(BASE_DIR, f"06_candidate_selection/{env_name}_20260420")
    os.makedirs(output_dir, exist_ok=True)

    # 1. 检查文件
    if not os.path.exists(arch_freq_file) or not os.path.exists(homology_detail_file):
        print(f"❌ 警告: 找不到 {env_name} 的输入文件，跳过...")
        return

    # 2. 读取数据
    df_arch = pd.read_csv(arch_freq_file)
    df_detail = pd.read_csv(homology_detail_file)
    
    # 3. 筛选高频且含新基因的组合
    mask = (df_arch['system_frequency'] >= MIN_OCCURRENCE) & (df_arch['has_novel_family'] == True)
    df_candidates = df_arch[mask].copy().sort_values(by='system_frequency', ascending=False)
    
    if df_candidates.empty:
        print(f"⚠️ {env_name} 中未找到符合条件的候选系统。")
        return

    candidate_summary_path = os.path.join(output_dir, f"{env_name}_candidate_systems_summary.csv")
    df_candidates.to_csv(candidate_summary_path, index=False)
    print(f"✅ 已保存候选架构摘要: 发现 {len(df_candidates)} 种高频架构")

    # 4. 收集需要提取的 target_cluster_ids
    target_cluster_ids = set()
    for _, row in df_candidates.iterrows():
        try:
            comp_tuple = ast.literal_eval(row['composition'])
            for cid in comp_tuple:
                target_cluster_ids.add(int(cid))
        except:
            pass

    # 建立 反查字典: Homology_ID -> 代表性 gene_name (每个簇挑第一个基因作为代表以提取序列)
    target_genes_to_extract = set()
    id_to_gene = {}
    for cluster_id in target_cluster_ids:
        genes_in_cluster = df_detail[df_detail['homology_cluster_id'] == cluster_id]['gene_name'].tolist()
        if genes_in_cluster:
            rep_gene = genes_in_cluster[0] # 取该家族的第一个基因作为代表序列
            target_genes_to_extract.add(rep_gene)
            id_to_gene[rep_gene] = cluster_id

    # 5. 提取序列
    print(f"正在从全体蛋白库中捞取 {len(target_cluster_ids)} 个家族的代表序列...")
    all_seqs = parse_fasta(all_proteins_file)
    
    output_fasta_path = os.path.join(output_dir, f"{env_name}_candidate_system_representatives.faa")
    extracted_count = 0
    
    with open(output_fasta_path, 'w') as out_f:
        for gene_name, seq in all_seqs.items():
            if gene_name in target_genes_to_extract:
                fam_id = id_to_gene[gene_name]
                # 重命名 Header，方便后续识别
                out_f.write(f">Fam_{fam_id}|{gene_name}\n{seq}\n")
                extracted_count += 1
                
    print(f"✅ 已保存代表序列 (共 {extracted_count} 条) 以供 AlphaFold/HHpred 分析")

    # 6. 生成详细上下文样本
    print("正在追溯原始物理簇上下文 (Context Sampling)...")
    
    # 获取物理簇信息 (利用 cluster_id 列，如果之前的表叫 physical_cluster_id，这里需要对应。我们假设 df_detail 中包含 phys_cluster_id)
    # 因为您的 df_detail 没有 physical_cluster_id，我们需要读 phy_csv
    phy_csv = os.path.join(BASE_DIR, f"02_physical_clustering/{env_name}_20260420/real_position_clusters_with_context_{env_name}.csv")
    df_phy = pd.read_csv(phy_csv)
    
    # 把 homology_id 映射进去
    gene_to_fam = dict(zip(df_detail['gene_name'], df_detail['homology_cluster_id']))
    df_phy['homology_cluster_id'] = df_phy['gene_name'].map(gene_to_fam).fillna(0).astype(int)

    samples = []
    for _, row in df_candidates.iterrows():
        comp_str = row['composition']
        target_fams = set(ast.literal_eval(comp_str))
        
        # 找出包含这些家族的基因组中的物理簇
        subset = df_phy[df_phy['homology_cluster_id'].isin(target_fams)]
        candidate_phy_ids = subset['cluster_id'].unique()
        
        found_count = 0
        for phy_id in candidate_phy_ids:
            if found_count >= 5: break 
            
            cluster_genes = df_phy[df_phy['cluster_id'] == phy_id].sort_values('start_position')
            present_fams = set(cluster_genes['homology_cluster_id'].astype(int))
            
            if target_fams.issubset(present_fams):
                genome_id = cluster_genes.iloc[0]['genome_id']
                gene_list_str = []
                for _, g in cluster_genes.iterrows():
                    fam = int(g['homology_cluster_id'])
                    is_target = " [★Target]" if fam in target_fams else ""
                    gene_list_str.append(f"{g['gene_name']}(Fam_{fam}{is_target})")
                
                samples.append({
                    'Architecture': comp_str,
                    'Frequency': row['system_frequency'],
                    'Novel_Families': row['novel_families'],
                    'Genome': genome_id,
                    'Cluster_ID': phy_id,
                    'Gene_Context': " --> ".join(gene_list_str)
                })
                found_count += 1

    df_samples = pd.DataFrame(samples)
    samples_path = os.path.join(output_dir, f"{env_name}_candidate_system_context_samples.csv")
    df_samples.to_csv(samples_path, index=False)
    print(f"✅ 已保存物理上下文样本表: {samples_path}")

if __name__ == "__main__":
    for env in ["TPMC-A", "TPMC-S"]:
        process_dataset(env)
    print("\n🎉 全部提取完成！现在您可以拿着 `.faa` 去跑 AlphaFold 或 HHpred 啦！")