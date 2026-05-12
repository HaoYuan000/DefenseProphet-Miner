#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文件名: system_architecture_analysis_TPMC-S.py
# 描述: 综合物理和同源聚类结果，分析系统组成模式，挖掘高置信度新系统候选并生成共现网络。

import pandas as pd
import os
import networkx as nx
from collections import Counter

def main():
    print("开始 TPMC-S 新防御系统深度分析与架构挖掘...")
    
    # === 1. 配置路径 ===
    base_dir = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang"
    
    # 输入文件 (依赖前面步骤的输出)
    physical_csv = os.path.join(base_dir, "02_physical_clustering/TPMC-S_20260420/real_position_clusters_with_context_TPMC-S.csv")
    homology_detail_csv = os.path.join(base_dir, "03_homology_clustering/TPMC-S_20260420/TPMC-S_homology_clusters_detail.csv")
    homology_stats_csv = os.path.join(base_dir, "03_homology_clustering/TPMC-S_20260420/TPMC-S_homology_cluster_statistics.csv")
    
    # 输出目录
    output_dir = os.path.join(base_dir, "04_system_analysis/TPMC-S_20260420")
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(physical_csv) or not os.path.exists(homology_detail_csv):
        print("错误: 找不到输入文件，请确保前面的物理聚类和同源聚类步骤已成功完成。")
        return

    # === 2. 加载数据 ===
    print("正在加载数据...")
    phy_df = pd.read_csv(physical_csv)
    hom_df = pd.read_csv(homology_detail_csv)
    hom_stats_df = pd.read_csv(homology_stats_csv)

    # 建立映射: gene_name -> homology_cluster_id
    gene_to_fam = dict(zip(hom_df['gene_name'], hom_df['homology_cluster_id']))
    
    # 将同源家族ID添加到物理聚类表中
    phy_df['homology_cluster_id'] = phy_df['gene_name'].map(gene_to_fam).fillna(0).astype(int)

    # === 3. 分析系统组成模式 (Architecture Analysis) ===
    print("分析防御系统组合架构...")
    
    # 按物理簇分组，收集每个簇包含的同源家族列表
    # 过滤掉 Homology ID 为 0 (未聚类) 的基因
    valid_phy = phy_df[phy_df['homology_cluster_id'] > 0]
    
    cluster_compositions = valid_phy.groupby('cluster_id')['homology_cluster_id'].apply(lambda x: sorted(list(x))).reset_index()
    
    # 统计每种组合出现的频率 (将列表转为元组以便哈希统计)
    composition_counts = Counter(tuple(x) for x in cluster_compositions['homology_cluster_id'])
    
    # 将结果转为 DataFrame
    comp_data = []
    for comp, count in composition_counts.items():
        # 判断这个组合里是否包含 Potential_Novel_Family
        is_novel = False
        novel_members = []
        for fam_id in set(comp):
            row = hom_stats_df[hom_stats_df['homology_cluster_id'] == fam_id]
            if not row.empty and row.iloc[0]['homology_type'] == 'Potential_Novel_Family':
                is_novel = True
                novel_members.append(fam_id)
        
        comp_data.append({
            'composition': str(comp),
            'family_count': len(comp),
            'system_frequency': count,
            'has_novel_family': is_novel,
            'novel_families': str(sorted(list(set(novel_members))))
        })
    
    comp_df = pd.DataFrame(comp_data).sort_values('system_frequency', ascending=False)
    comp_df.to_csv(os.path.join(output_dir, "TPMC-S_system_architecture_frequency.csv"), index=False)
    
    # === 4. 挖掘核心候选系统 (Core Novel Candidates) ===
    print("提炼高置信度新系统重现模式...")
    
    # 筛选条件：
    # 1. 必须包含 'Potential_Novel_Family' 类型的同源家族
    # 2. 该组合出现频率 >= 2 (在不同基因组中重复出现)
    target_systems = comp_df[(comp_df['has_novel_family'] == True) & (comp_df['system_frequency'] >= 2)]
    
    print(f"发现 {len(target_systems)} 种跨基因组重复出现且包含潜在新基因的系统架构")
    
    # 生成详细报告
    report_file = os.path.join(output_dir, "TPMC-S_novel_system_candidates_report.txt")
    with open(report_file, "w", encoding='utf-8') as f:
        f.write("=== TPMC-S 高置信度潜在新防御系统候选挖掘报告 ===\n\n")
        f.write(f"筛选标准：跨基因组出现 >= 2次，且含有全新同源家族组成的架构。\n")
        f.write("="*60 + "\n\n")
        
        for idx, row in target_systems.iterrows():
            comp_str = row['composition'] 
            freq = row['system_frequency']
            novel_fams = row['novel_families']
            
            f.write(f"【系统架构组合】 出现频次: {freq} 次\n")
            f.write(f"  家族组成ID: {comp_str}\n")
            f.write(f"  包含的新颖家族ID: {novel_fams}\n")
            
            target_comp = eval(comp_str)
            
            # 找到匹配该架构的物理簇ID
            matching_clusters = []
            for _, c_row in cluster_compositions.iterrows():
                if tuple(c_row['homology_cluster_id']) == target_comp:
                    matching_clusters.append(c_row['cluster_id'])
                    if len(matching_clusters) >= 3: break # 每个架构只展示前3个例子
            
            f.write("  代表性物理簇及基因明细:\n")
            for cid in matching_clusters:
                genes = phy_df[phy_df['cluster_id'] == cid]
                genome_id = genes.iloc[0]['genome_id']
                f.write(f"    ▶ 簇ID: {cid} (归属基因组: {genome_id})\n")
                for _, g in genes.iterrows():
                    fam_id = g['homology_cluster_id']
                    fam_info = hom_stats_df[hom_stats_df['homology_cluster_id'] == fam_id]
                    fam_type = fam_info.iloc[0]['homology_type'] if not fam_info.empty else "Unknown"
                    if fam_type == 'Potential_Novel_Family':
                        mark = "[★ NOVEL]"
                    elif fam_type == 'Expanded_Known_Family':
                        mark = "[+ EXPAND]"
                    else:
                        mark = "[  KNOWN]"
                        
                    f.write(f"      {mark} {g['gene_name']} (Family: {fam_id}, Type: {g['gene_type']})\n")
            f.write("\n" + "-"*60 + "\n\n")

    # === 5. 共现网络分析 (Co-occurrence Network) ===
    print("构建同源家族共现网络以导入 Cytoscape ...")
    G = nx.Graph()
    
    for comp_tuple, count in composition_counts.items():
        # 在同一个物理簇内的家族两两之间建立连接
        unique_fams = sorted(list(set(comp_tuple)))
        for i in range(len(unique_fams)):
            for j in range(i+1, len(unique_fams)):
                u, v = unique_fams[i], unique_fams[j]
                if G.has_edge(u, v):
                    G[u][v]['weight'] += count
                else:
                    G.add_edge(u, v, weight=count)

    # 导出网络边列表
    edge_data = []
    for u, v, data in G.edges(data=True):
        u_info = hom_stats_df[hom_stats_df['homology_cluster_id'] == u]
        v_info = hom_stats_df[hom_stats_df['homology_cluster_id'] == v]
        
        u_type = u_info.iloc[0]['homology_type'] if not u_info.empty else "Unknown"
        v_type = v_info.iloc[0]['homology_type'] if not v_info.empty else "Unknown"
        
        edge_data.append({
            'Source': u, 'Target': v, 'Weight': data['weight'],
            'Source_Type': u_type, 'Target_Type': v_type
        })
    
    if edge_data:
        pd.DataFrame(edge_data).to_csv(os.path.join(output_dir, "TPMC-S_family_cooccurrence_network.csv"), index=False)

    print(f"\n所有系统深度分析成功完成!")
    print(f"输出目录: {output_dir}")
    print(f"请重点查阅详细文字报告: {report_file}")

if __name__ == "__main__":
    main()
