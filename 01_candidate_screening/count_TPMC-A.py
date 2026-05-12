import os
import pandas as pd
from tqdm import tqdm
from collections import Counter
import concurrent.futures

def process_file(args):
    """处理单个CSV文件的函数"""
    file, results_folder = args
    csv_path = os.path.join(results_folder, file)
    
    local_comparison_results = []
    local_total_proteins = 0
    local_label_counts = Counter()
    
    try:
        # 读取CSV文件，没有表头
        df = pd.read_csv(csv_path, header=None, on_bad_lines='skip')
        local_total_proteins = len(df)
        
        # 跳过列数严重不足的异常数据
        if len(df.columns) >= 3:
            # 统计二分类标签分布 (第3列，索引为2)
            cleaned_labels = df[2].dropna().astype(str).str.strip()
            local_label_counts.update(cleaned_labels.tolist())
            
            # 如果列数满足详细筛选的要求，则进行检查
            if len(df.columns) >= 6:
                for index, row in df.iterrows():
                    gene_name = str(row[1]).strip()
                    binary_label = str(row[2]).strip()
                    multi_label = str(row[4]).strip()
                    
                    try:
                        binary_prob = float(row[3])
                    except ValueError:
                        continue
                        
                    # 条件：二分类是defense 且 置信度>0.99 且 多分类为 'other'
                    if binary_label == 'defense' and binary_prob > 0.99 and multi_label == 'other':
                        local_comparison_results.append({
                            'gene_name': gene_name,
                            'multi_label': multi_label,
                            'binary_label': binary_label,
                            'binary_confidence': binary_prob,
                            'file_name': file
                        })
                        
        return True, file, local_total_proteins, local_label_counts, local_comparison_results

    except Exception as e:
        return False, file, str(e), Counter(), []

if __name__ == '__main__':
    # 文件路径
    results_folder = '/home/xiongxinghao/data6/project/anti-virus/infer_model/results/TPMC_A_predictions_results'
    output_path = '/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/01_candidate_screening/archaea_comparison_results_high_conf.csv'
    stats_output_path = '/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang/01_candidate_screening/archaea_prediction_stats.txt'

    # 检查文件夹是否存在
    if not os.path.exists(results_folder):
        print(f"文件夹不存在：{results_folder}")
        exit(1)

    # --- 初始化变量 ---
    comparison_results = []
    total_genomes = 0
    total_proteins = 0
    label_counts = Counter()

    # 获取所有需要处理的CSV文件
    csv_files = [f for f in os.listdir(results_folder) if f.endswith('.csv')]
    tasks = [(f, results_folder) for f in csv_files]

    # --- 多进程处理 ---
    # 使用操作系统默认的 CPU 核心数。如果想指定核心数，例如8核，可以修改为 max_workers=8
    max_workers = 32
    print(f"开始处理，使用核心数/进程数: {max_workers}")

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # submit 并获取 future 结果
        futures = {executor.submit(process_file, task): task for task in tasks}
        
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(tasks), desc="多进程度处理与筛选"):
            success, file, proteins_or_error, local_counts, local_results = future.result()
            
            if success:
                total_genomes += 1
                total_proteins += proteins_or_error
                label_counts.update(local_counts)
                comparison_results.extend(local_results)
            else:
                print(f"\n处理文件出错：{file}, 错误信息：{proteins_or_error}")

    # --- 输出统计结果 ---
    print("\n" + "="*40)
    print(f"总共成功处理 {total_genomes} 个文件 (总计找到 {len(csv_files)} 个csv)")
    print(f"总共包含 {total_proteins} 个蛋白质预测")
    print("\n[二分类标签分布]:")

    stats_content = []
    stats_content.append("="*40)
    stats_content.append(f"总共成功处理 {total_genomes} 个文件 (总计找到 {len(csv_files)} 个csv)")
    stats_content.append(f"总共包含 {total_proteins} 个蛋白质预测")
    stats_content.append("\n[二分类标签分布]:")

    if total_proteins > 0:
        for label, count in label_counts.most_common():
            stat_line = f"  {label}: {count} ({count/total_proteins*100:.2f}%)"
            print(stat_line)
            stats_content.append(stat_line)
    else:
        print("没有找到蛋白质预测数据。")
        stats_content.append("没有找到蛋白质预测数据。")
        
    with open(stats_output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(stats_content) + '\n')
    print(f"\n统计结果已保存到：{stats_output_path}")

    # --- 保存筛选结果 ---
    print("-" * 40)
    if comparison_results:
        output_df = pd.DataFrame(comparison_results)
        output_df.to_csv(output_path, index=False)
        print(f"\n对比结果已保存到：{output_path}")
        print(f"共筛选出 {len(comparison_results)} 个高置信度(>0.99)二分类为defense且多分类为other的基因。")
    else:
        print("\n没有找到符合条件的高置信度基因。")
    print("="*40 + "\n")