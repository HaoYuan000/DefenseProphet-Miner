import pandas as pd
import os
import ast
from collections import Counter

# ================= 配置路径 =================
BASE_DIR = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang"

def process_known_functions(env_name):
    print(f"\n{'='*50}")
    print(f"🔍 正在解析系统已知组件功能: {env_name}")
    print(f"{'='*50}")
    
    # 1. 候选系统列表 (由 extract_top_candidates.py 生成)
    cand_summary_file = os.path.join(BASE_DIR, f"06_candidate_selection/{env_name}_20260420/{env_name}_candidate_systems_summary.csv")
    
    # 2. 同源聚类映射表 (包含 gene_name <-> homology_cluster_id)
    hom_detail_file = os.path.join(BASE_DIR, f"03_homology_clustering/{env_name}_20260420/{env_name}_homology_clusters_detail.csv")
    
    # 3. 原始的高置信度注释表 (包含 gene_name <-> multi_label)
    raw_label_file = os.path.join(BASE_DIR, f"01_candidate_screening/{env_name}_20260415/{env_name}_all_defense.csv")

    output_dir = os.path.join(BASE_DIR, f"06_candidate_selection/{env_name}_20260420")
    output_file = os.path.join(output_dir, f"{env_name}_candidate_known_functions.csv")

    # 检查文件完整性
    for fpath in [cand_summary_file, hom_detail_file, raw_label_file]:
        if not os.path.exists(fpath):
            print(f"❌ 错误: 找不到依赖文件 {fpath}，跳过该环境。")
            return

    df_cand = pd.read_csv(cand_summary_file)
    df_hom = pd.read_csv(hom_detail_file)
    df_raw = pd.read_csv(raw_label_file)

    # 预处理：建立查找字典 gene_name -> biological_label
    # 兼容处理列名 (可能是 multi_label 或 predicted_label)
    label_col = 'multi_label' if 'multi_label' in df_raw.columns else 'predicted_label'
    gene_to_label = dict(zip(df_raw['gene_name'], df_raw[label_col]))

    results = []
    print(f"正在翻译 {len(df_cand)} 个候选系统组合的已知组件...")

    for _, row in df_cand.iterrows():
        try:
            composition = ast.literal_eval(row['composition']) 
            novel_fams = ast.literal_eval(row['novel_families'])
        except:
            continue

        # 找出 Known IDs (即存在于架构中，但不属于新颖家族的ID)
        known_ids = [fid for fid in composition if fid not in novel_fams]
        
        if not known_ids:
            known_desc = "Pure_Novel_System"
        else:
            func_descriptions = []
            for k_id in known_ids:
                # 找到该同源簇下的所有基因
                genes_in_cluster = df_hom[df_hom['homology_cluster_id'] == k_id]['gene_name'].tolist()
                
                # 查这些基因的防御系统名称
                labels = []
                for g in genes_in_cluster:
                    lbl = gene_to_label.get(g)
                    # 忽略 'other', 'Other', 'unknown' 等非具体名称
                    if lbl and str(lbl).lower() != "other" and str(lbl).lower() != "unknown": 
                        labels.append(lbl)
                
                # 统计最常见的防御类型
                if labels:
                    most_common = Counter(labels).most_common(1)[0]
                    func_descriptions.append(f"Fam_{k_id}={most_common[0]}")
                else:
                    func_descriptions.append(f"Fam_{k_id}=Unknown_Defense")
            
            known_desc = "; ".join(func_descriptions)

        results.append({
            "Architecture": row['composition'],
            "Frequency": row['system_frequency'],
            "Novel_Subunits": str(novel_fams),
            "Anchored_Known_System": known_desc
        })

    # 输出并保存
    df_res = pd.DataFrame(results)
    df_res.to_csv(output_file, index=False)
    
    print("-" * 50)
    print("✅ 分析完成！预览 Top 5 结果：")
    print(df_res.head(5).to_string(index=False))
    print("-" * 50)
    print(f"结果已保存至: {output_file}")

if __name__ == "__main__":
    for env in ["TPMC-A", "TPMC-S"]:
        process_known_functions(env)