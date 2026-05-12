import pandas as pd
import os
import re
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import multiprocessing as mp

def extract_position_from_header(header):
    """从FASTA头部提取位置信息"""
    # 匹配GTDB/Prodigal默认格式：序列ID # 开始位置 # 结束位置 # 方向 # 其他信息
    match = re.search(r'# (\d+) # (\d+) # [+-]?[01] #', header)
    if match:
        start = int(match.group(1))
        end = int(match.group(2))
        return (start, end)
    
    # 尝试其他格式
    patterns = [
        r'location=(\d+)\.\.(\d+)',
        r'\[location=(\d+)\.\.(\d+)\]',
        r'(\d+)\.\.(\d+)'
    ]
    
    for pattern in patterns:
        match = re.search(pattern, header)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            return (start, end)
    
    return None

def _process_single_fasta(args):
    """多进程处理单个文件的辅助函数"""
    prefix, protein_file = args
    genes_in_file = []
    try:
        for record in SeqIO.parse(protein_file, "fasta"):
            position_info = extract_position_from_header(record.description)
            
            # 尝试从ID中提取contig信息
            gene_id = record.id
            match = re.search(r'(.+)_\d+$', gene_id)
            contig_id = match.group(1) if match else gene_id
            
            if position_info:
                start_pos, end_pos = position_info
                genes_in_file.append({
                    'gene_id': record.id,
                    'contig_id': contig_id, 
                    'start': start_pos,
                    'end': end_pos,
                    'record_id': record.id,
                    'record_desc': record.description,
                    'record_seq': str(record.seq) # 转为字符串避免多进程序列化报错
                })
        if genes_in_file:
            genes_in_file.sort(key=lambda x: x['start'])
        return prefix, protein_file, genes_in_file, None
    except Exception as e:
        return prefix, protein_file, [], str(e)


def find_all_gene_positions(genome_dir, target_file_prefixes, threads=64):
    """
    查找指定前缀文件中的【所有】基因的位置信息 (多核并行版本)
    """
    # 查找蛋白质文件 (支持 .faa 或 .fasta)
    protein_files = {}
    for root, dirs, files in os.walk(genome_dir):
        for file in files:
            if file.endswith('.fasta') or file.endswith('.faa'):
                file_path = os.path.join(root, file)
                for prefix in target_file_prefixes:
                    # 确保前缀匹配文件名本体，避免跨目录误匹配
                    if file.startswith(prefix):
                        protein_files[prefix] = file_path
                        break 
    
    print(f"找到 {len(protein_files)} 个涉及的基因组文件")
    
    genome_data = {}
    worker_args = [(prefix, target_path) for prefix, target_path in protein_files.items()]
    
    print(f"启动 {threads} 个进程并发解析基因组序列...")
    with mp.Pool(processes=threads) as pool:
        for prefix, protein_file, genes_in_file, err_msg in pool.imap_unordered(_process_single_fasta, worker_args):
            if err_msg:
                print(f"处理文件 {protein_file} 时出错: {err_msg}")
            
            if genes_in_file:
                genome_data[prefix] = {
                    'genes': genes_in_file,
                    'file_path': protein_file
                }
            
    return genome_data

def analyze_defense_genes_with_real_positions(candidate_csv_path, defense_pool_csv_path, genome_dir, output_dir, max_gap=1500):
    """
    使用基因的实际位置分析防御基因并识别基因簇。
    逻辑：
    1. 仅计算 binary_label='defense' 的基因。
    2. 如果两个防御基因间距 <= max_gap，则连起来（忽略中间的非防御基因）。
    3. 结果簇必须 >= 2个基因，且必须包含 candidate_csv 中的基因。
    """
    
    # 1. 读取基础防御基因池 (Pool) - 用于构建已知骨架 (范围更广的集合，即所有的 defense 基因)
    print(f"读取防御基因池: {defense_pool_csv_path}")
    pool_df = pd.read_csv(defense_pool_csv_path)
    # 严格筛选：只有 binary_label 为 defense 的基因被视为有效的防御系统组件
    defense_pool_df = pool_df[pool_df['binary_label'] == 'defense'].copy()
    
    # 已知防御基因集合 (Known Defense)
    known_defense_set = set(defense_pool_df['gene_name'].astype(str))
    
    # 建立防御基因的元数据查找表
    defense_info_map = {}
    for _, row in defense_pool_df.iterrows():
        defense_info_map[row['gene_name']] = row.to_dict()

    print(f"参与骨架构建的 Defense 基因总数 (Pool): {len(known_defense_set)}")

    # 2. 读取候选基因 (Candidates) - 这是我们必须寻找的核心目标 (高置信度筛选出来的基因)
    print(f"读取候选靶标基因列表: {candidate_csv_path}")
    cand_df = pd.read_csv(candidate_csv_path)
    # 强制将 candidate 标识为 defense
    cand_df['binary_label'] = 'defense'
    candidate_gene_set = set(cand_df['gene_name'].astype(str))
    print(f"必须包含的罕见/高置信度候选基因数 (Candidate): {len(candidate_gene_set)}")
    
    # 将候选基因也加入元数据表，如果已存在则更新它（以 candidate 的信息为主）
    for _, row in cand_df.iterrows():
        defense_info_map[row['gene_name']] = row.to_dict()

    # 合并两个集合，定义所有允许参与聚类的节点
    valid_clustering_genes = known_defense_set.union(candidate_gene_set)


    # 获取涉及的文件列表
    file_prefixes = set()
    all_relevant_df = pd.concat([defense_pool_df, cand_df]) if not cand_df.empty else defense_pool_df
    
    for file_name in all_relevant_df['file_name'].unique():
        if isinstance(file_name, str):
            # 去掉 .csv
            prefix = file_name.split('.csv')[0]
            # 去掉预测脚本自动加的后缀，还原真实的基因组名前缀
            prefix = prefix.replace('_predictions', '')
            prefix = prefix.split('_protein_labels.pkl')[0]
            file_prefixes.add(prefix)
            
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    sequences_dir = os.path.join(output_dir, "defense_clusters_sequences")
    os.makedirs(sequences_dir, exist_ok=True)
    
    # 3. 扫描基因组获取坐标
    print("正在扫描基因组文件以获取基因坐标...")
    genome_data = find_all_gene_positions(genome_dir, file_prefixes)
    
    real_position_clusters = []
    
    print("开始仅基于 Defense Genes 进行聚类 (Single Linkage Clustering)...")
    
    for prefix, data in genome_data.items():
        all_genes_in_file = data['genes']
        file_path = data['file_path']
        file_basename = os.path.basename(file_path)
        
        # 【核心逻辑】：只提取有效的 Defense/Candidate Genes 作为节点
        valid_nodes = []
        for g in all_genes_in_file:
            gid = g['gene_id']
            if gid in valid_clustering_genes:
                valid_nodes.append(g)
        
        valid_nodes.sort(key=lambda x: x['start'])
        
        # 按 Contig 分组
        genes_by_contig = defaultdict(list)
        for gene in valid_nodes:
            genes_by_contig[gene['contig_id']].append(gene)
            
        for contig_id, node_genes in genes_by_contig.items():
            if not node_genes:
                continue
                
            current_cluster = [node_genes[0]]
            
            for i in range(1, len(node_genes)):
                prev_gene = current_cluster[-1]
                current_gene = node_genes[i]
                
                # 计算两个防御节点之间的物理距离
                distance = current_gene['start'] - prev_gene['end']
                
                if distance <= max_gap:
                    current_cluster.append(current_gene)
                else:
                    process_cluster(current_cluster, candidate_gene_set, defense_info_map, 
                                  real_position_clusters, prefix, file_basename, contig_id)
                    current_cluster = [current_gene]
            
            # 处理最后一个簇
            if current_cluster:
                process_cluster(current_cluster, candidate_gene_set, defense_info_map, 
                              real_position_clusters, prefix, file_basename, contig_id)

    print(f"找到 {len(real_position_clusters)} 个符合条件的防御系统簇")
    
    # 保存结果
    if real_position_clusters:
        save_cluster_csv(real_position_clusters, output_dir)
        extract_cluster_sequences(real_position_clusters, sequences_dir)
        create_summary_report(real_position_clusters, output_dir)
    else:
        print("未找到任何符合条件的簇。")
    
    return real_position_clusters

def process_cluster(cluster_genes, candidate_set, info_map, results_list, prefix, file_name, genome_id):
    """
    处理一个物理簇，判断是否保留并添加标签
    """
    if not cluster_genes:
        return

    # 【条件1】：总数必须大于等于2
    if len(cluster_genes) < 2:
        return

    # 统计簇内成分
    has_candidate = False
    candidates_in_cluster = []
    has_known_def = False
    cluster_rows = []
    
    for gene in cluster_genes:
        gid = gene['gene_id']
        is_cand = gid in candidate_set
        
        if is_cand:
            has_candidate = True
            candidates_in_cluster.append(gid)

        if gid in info_map:
            row = info_map[gid].copy()
            row['type'] = 'candidate' if is_cand else 'known_defense'
        else:
            row = {
                'gene_name': gid,
                'file_name': file_name,
                'multi_label': 'Target_Candidate' if is_cand else 'Known_Context',
                'binary_label': 'defense',
                'binary_confidence': 1.0,
                'type': 'candidate' if is_cand else 'known_defense'
            }
        cluster_rows.append(row)

        if not is_cand:
            # 只有 multi_label 不是 'other' 且有明确注释的防卫基因才被视为 KNOWN_DEF
            if str(row.get('multi_label', '')).lower() != 'other':
                has_known_def = True
    
    # 【条件2】：必须包含候选基因
    if has_candidate:
        if has_known_def:
            cluster_type = "HYBRID" 
        else:
            cluster_type = "PURE_CANDIDATE" 

        cluster_info = {
            'cluster_type': cluster_type, 
            'genome_id': genome_id,
            'file_prefix': prefix,
            'file_name': file_name,
            'genes': [g['gene_id'] for g in cluster_genes],
            'positions': [(g['start'], g['end']) for g in cluster_genes],
            'records': [{'id': g['record_id'], 'description': g['record_desc'], 'seq': g['record_seq']} for g in cluster_genes],
            'rows': cluster_rows,
            'gene_count': len(cluster_genes),
            'candidate_count': len(candidates_in_cluster)
        }
        results_list.append(cluster_info)

def save_cluster_csv(clusters, output_dir):
    real_pos_csv = os.path.join(output_dir, "real_position_clusters_with_context_TPMC-A.csv")
    data = []
    
    for cluster in clusters:
        start_pos = min([pos[0] for pos in cluster['positions']])
        end_pos = max([pos[1] for pos in cluster['positions']])
        cluster_id = f"{cluster['genome_id']}_{start_pos}-{end_pos}"
        
        for i, gene_name in enumerate(cluster['genes']):
            row = cluster['rows'][i]
            position = cluster['positions'][i]
            
            entry = {
                'cluster_id': cluster_id,
                'cluster_type': cluster['cluster_type'],
                'genome_id': cluster['genome_id'],
                'file_name': cluster['file_name'],
                'gene_name': gene_name,
                'start_position': position[0],
                'end_position': position[1],
                'gene_type': row['type'],
                'system_size': cluster['gene_count'],
                'candidate_count': cluster['candidate_count'],
                'multi_label': row.get('multi_label', ''),
                'binary_confidence': row.get('binary_confidence', 0)
            }
            data.append(entry)
            
    df = pd.DataFrame(data)
    df.to_csv(real_pos_csv, index=False)
    print(f"CSV已保存到 {real_pos_csv}")

def extract_cluster_sequences(clusters, output_dir):
    for cluster in clusters:
        start_pos = min([pos[0] for pos in cluster['positions']])
        end_pos = max([pos[1] for pos in cluster['positions']])
        cluster_id = f"{cluster['genome_id']}_{start_pos}-{end_pos}"
        
        output_file = os.path.join(output_dir, f"{cluster_id}.faa")
        
        with open(output_file, "w") as out_f:
            for i, record in enumerate(cluster['records']):
                gene_type = cluster['rows'][i]['type']
                out_f.write(f">{record['id']} type={gene_type} cluster_type={cluster['cluster_type']} {record['description']}\n{record['seq']}\n")

def create_summary_report(clusters, output_dir):
    summary_file = os.path.join(output_dir, "real_position_clusters_summary.txt")
    with open(summary_file, "w") as f:
        f.write(f"基于实际位置筛选，找到 {len(clusters)} 个包含候选基因的簇\n")
        f.write(f"筛选条件：1. 仅包含 defense 基因  2. 必须包含 candidate  3. 基因数 >= 2\n")
        
        hybrid_count = sum(1 for c in clusters if c['cluster_type'] == 'HYBRID')
        pure_count = sum(1 for c in clusters if c['cluster_type'] == 'PURE_CANDIDATE')
        f.write(f"统计：混合簇 (HYBRID): {hybrid_count} 个，纯新簇 (PURE_CANDIDATE): {pure_count} 个\n\n")
        
        for i, cluster in enumerate(clusters):
            start_pos = min([pos[0] for pos in cluster['positions']])
            end_pos = max([pos[1] for pos in cluster['positions']])
            cluster_id = f"{cluster['genome_id']}_{start_pos}-{end_pos}"
            
            f.write(f"簇 {i+1}: {cluster_id}  [{cluster['cluster_type']}]\n")
            f.write(f"基因组: {cluster['genome_id']}\n")
            f.write(f"总基因数: {cluster['gene_count']} (Candidates: {cluster['candidate_count']})\n")
            f.write(f"位置范围: {start_pos} - {end_pos}\n")
            f.write("基因列表:\n")
            
            for j, gene_name in enumerate(cluster['genes']):
                pos = cluster['positions'][j]
                row = cluster['rows'][j]
                
                label = row.get('multi_label', 'Unknown')
                
                if row['type'] == 'candidate':
                    type_mark = "[★ CANDIDATE]"
                elif str(label).lower() == 'other':
                    type_mark = "[? POTENTIAL]"
                else:
                    type_mark = "[  KNOWN_DEF]"
                    
                conf_str = f"conf: {float(row.get('binary_confidence', 0)):.4f}" 
                
                f.write(f"  {type_mark} {gene_name} | Pos: {pos[0]}-{pos[1]} | {label} | {conf_str}\n")
            
            f.write("\n")

    print(f"汇总报告已保存到 {summary_file}")

if __name__ == "__main__":
    # 【已互换】目标：寻找高置信度(>0.99)且类别为other的新颖防卫基因 (严格筛选出的靶标)
    candidate_csv = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/01_candidate_screening/TPMC-A_20260415/TPMC-A_comparison_results_high_conf.csv"
    
    # 【已互换】组件池：所有确认为防卫基因的集合，用于建立聚类的架构 (基础盘)
    defense_pool_csv = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/01_candidate_screening/TPMC-A_20260415/TPMC-A_all_defense.csv"
    
    genome_dir = "/data5/zhanghaohong/projects/TPMC_BGC/prodigal/TPMC-A"
    output_dir = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/02_physical_clustering/TPMC-A_20260420"
    
    # max_gap=1500 表示两个防御基因之间最大间隔。
    analyze_defense_genes_with_real_positions(candidate_csv, defense_pool_csv, genome_dir, output_dir, max_gap=1500)
    print("分析完成!")