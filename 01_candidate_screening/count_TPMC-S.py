import os
import pandas as pd
from tqdm import tqdm
from collections import Counter
import concurrent.futures

def process_file(args):
    """处理单个CSV文件的函数"""
    file, results_folder = args
    csv_path = os.path.join(results_folder, file)
    
    local_all_defense = []      # 存放所有 defense 基因
    local_high_conf = []        # 存放高置信度特定基因
    local_total_proteins = 0
    local_label_counts = Counter()
    
    try:
        # 读取CSV文件，指定 header=0 过滤掉第一行的列名
        df = pd.read_csv(csv_path, header=0, on_bad_lines='skip')
        local_total_proteins = len(df)
        
        # 跳过列数严重不足的异常数据
        if len(df.columns) >= 3:
            # 统计二分类标签分布 (第3列，索引为2，注意去掉header后原来的行不再算作数据)
            cleaned_labels = df.iloc[:, 2].dropna().astype(str).str.strip()
            local_label_counts.update(cleaned_labels.tolist())
            
            # 如果列数满足详细筛选的要求，则进行检查
            if len(df.columns) >= 6:
                for index, row in df.iterrows():
                    gene_name = str(row.iloc[1]).strip()
                    binary_label = str(row.iloc[2]).strip()
                    multi_label = str(row.iloc[4]).strip()
                    
                    try:
                        binary_prob = float(row.iloc[3])
                    except ValueError:
                        continue
                        
                    # 只要是 defense 就提取并记录
                    if binary_label == 'defense':
                        result_dict = {
                            'gene_name': gene_name,
                            'multi_label': multi_label,
                            'binary_label': binary_label,
                            'binary_confidence': binary_prob,
                            'file_name': file
                        }
                        local_all_defense.append(result_dict)
                        
                        # 在 defense 的基础上，判断是否满足高置信度且多分类为other的强条件
                        if binary_prob > 0.99 and multi_label == 'other':
                            local_high_conf.append(result_dict)
                        
        return True, file, local_total_proteins, local_label_counts, local_all_defense, local_high_conf

    except Exception as e:
        return False, file, str(e), Counter(), [], []


if __name__ == '__main__':
    # 文件路径 - 修改为 TPMC-S 相关路径
    results_folder = '/home/xiongxinghao/data6/project/anti-virus/infer_model/results/TPMC_S_predictions_results'
    
    # 输出路径统一加上 bacteria 前缀以区分（如果是其他物种可对应修改）
    all_defense_output = '/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/01_candidate_screening/TPMC-S/TPMC-S_all_defense.csv'
    high_conf_output = '/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/01_candidate_screening/TPMC-S//TPMC-S_comparison_results_high_conf.csv'
    stats_output_path = '/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/01_candidate_screening/TPMC-S/TPMC-S_prediction_stats.txt'

    # 检查文件夹是否存在
    if not os.path.exists(results_folder):
        print(f"文件夹不存在：{results_folder}")
        exit(1)

    # --- 初始化变量 ---
    global_all_defense = []
    global_high_conf = []
    total_genomes = 0
    total_proteins = 0
    label_counts = Counter()

    # 获取所有需要处理的CSV文件
    csv_files = [f for f in os.listdir(results_folder) if f.endswith('.csv')]
    tasks = [(f, results_folder) for f in csv_files]

    # --- 多进程处理 ---
    max_workers = 32
    print(f"开始处理 TPMC-S 数据，使用核心数/进程数: {max_workers}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_file, task): task for task in tasks}
        
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(tasks), desc="多进程处理与筛选"):
            success, file, proteins_or_error, local_counts, results_all, results_high = future.result()
            
            if success:
                total_genomes += 1
                total_proteins += proteins_or_error
                label_counts.update(local_counts)
                global_all_defense.extend(results_all)
                global_high_conf.extend(results_high)
            else:
                print(f"\n处理文件出错：{file}, 错误信息：{proteins_or_error}")

    # --- 输出统计结果 ---
    print("\n" + "="*40)
    print(f"总共成功处理 {total_genomes} 个文件 (总计找到 {len(csv_files)} 个csv)")
    print(f"总共包含 {total_proteins} 个蛋白质预测 (已剔除表头行)")
    print("\n[二分类标签分布]:")

    stats_content = []
    stats_content.append("="*40)
    stats_content.append(f"总共成功处理 {total_genomes} 个文件 (总计找到 {len(csv_files)} 个csv)")
    stats_content.append(f"总共包含 {total_proteins} 个蛋白质预测 (已剔除表头行)")
    stats_content.append("\n[二分类标签分布]:")

    if total_proteins > 0:
        for label, count in label_counts.most_common():
            stat_line = f"  {label}: {count} ({count/total_proteins*100:.2f}%)"
            print(stat_line)
            stats_content.append(stat_line)
    else:
        print("没有找到蛋白质预测数据。")
        stats_content.append("没有找到蛋白质预测数据。")
        
    # 追加提取结果的数据到统计文件
    stats_content.append("\n" + "-"*40)
    stats_content.append(f"共提取出 {len(global_all_defense)} 个二分类为 defense 的基因。")
    stats_content.append(f"共筛选出 {len(global_high_conf)} 个高置信度(>0.99)且多分类为other的 defense 基因。")

    with open(stats_output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(stats_content) + '\n')
    print(f"\n统计结果已保存到：{stats_output_path}")

    # --- 保存提取结果 ---
    print("-" * 40)
    if global_all_defense:
        pd.DataFrame(global_all_defense).to_csv(all_defense_output, index=False)
        print(f"所有 [defense] 基因结果已保存到：{all_defense_output} (共 {len(global_all_defense)} 个)")
    else:
        print("\n没有找到任何二分类为 defense 的基因。")
        
    if global_high_conf:
        pd.DataFrame(global_high_conf).to_csv(high_conf_output, index=False)
        print(f"[高置信度+other] 基因结果已保存到：{high_conf_output} (共 {len(global_high_conf)} 个)")
    else:
        print("\n没有找到符合高置信度条件的基因。")
    print("="*40 + "\n")