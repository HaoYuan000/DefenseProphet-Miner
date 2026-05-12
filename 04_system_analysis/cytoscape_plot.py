import pandas as pd
import os

def prepare_cytoscape_files():
    # 路径配置
    base_dir = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang"
    raw_network_csv = os.path.join(base_dir, "04_system_analysis/TPMC-S_20260420/TPMC-S_family_cooccurrence_network.csv")
    stats_csv = os.path.join(base_dir, "03_homology_clustering/TPMC-S_20260420/TPMC-S_homology_cluster_statistics.csv")
    out_dir = os.path.join(base_dir, "05_visualizations")
    os.makedirs(out_dir, exist_ok=True)
    
    # 1. 过滤网络表 (去除杂乱低频连线)
    df_edges = pd.read_csv(raw_network_csv)
    # 设定阈值，只保留共现 >= 3 次的组合，瞬间清爽！
    filtered_edges = df_edges[df_edges['Weight'] >= 3].copy()
    out_edges = os.path.join(out_dir, "Cytoscape_Edges_Filtered.csv")
    filtered_edges.to_csv(out_edges, index=False)
    
    # 2. 获取当前网络中留下的所有节点 ID
    nodes_in_network = set(filtered_edges['Source']).union(set(filtered_edges['Target']))
    
    # 3. 制作节点注释表
    df_stats = pd.read_csv(stats_csv)
    node_data = []
    
    for node_id in nodes_in_network:
        row = df_stats[df_stats['homology_cluster_id'] == node_id]
        if not row.empty:
            node_type = row.iloc[0]['cluster_type']
            # 获取代表性注释，如果是全新家族就是 Novel，如果是已知就提取其防卫类型
            annotations = str(row.iloc[0].get('representative_annotations', 'Unknown'))
            if node_type == 'PURE_CANDIDATE':
                label = f"★ Novel_F{node_id}"
                cat = "Novel"
            else:
                # 简单截取注释前30个字符作为名字，防止太长
                short_anno = annotations.split(',')[0][:30] if annotations else f"Known_F{node_id}"
                label = f"F{node_id}: {short_anno}"
                cat = "Known"
        else:
            label = f"Fam_{node_id}"
            cat = "Unknown"
            
        node_data.append({'Node_ID': node_id, 'Label': label, 'Category': cat})
        
    df_nodes = pd.DataFrame(node_data)
    out_nodes = os.path.join(out_dir, "Cytoscape_Nodes_Annotated.csv")
    df_nodes.to_csv(out_nodes, index=False)
    
    print(f"✅ Cytoscape 文件已生成在 {out_dir}:")
    print(f"   1. 边文件: Cytoscape_Edges_Filtered.csv ({len(filtered_edges)} 条边)")
    print(f"   2. 节点文件: Cytoscape_Nodes_Annotated.csv ({len(nodes_in_network)} 个节点)")

if __name__ == "__main__":
    prepare_cytoscape_files()