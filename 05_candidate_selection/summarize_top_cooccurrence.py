import csv
import os
from collections import defaultdict

# ================= 配置路径 =================
BASE_DIR = "/home/xiongxinghao/data6/project/anti-virus/liyuanhao/qinzang"

def analyze_top_cooccurrence(env_name):
    input_file = os.path.join(BASE_DIR, f"06_candidate_selection/{env_name}_20260420/{env_name}_candidate_known_functions.csv")
    
    if not os.path.exists(input_file):
        print(f"❌ 错误: 找不到 {input_file}。请先运行 annotate_candidate_functions.py。")
        return

    freq_map = defaultdict(int)
    
    with open(input_file, mode='r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            known_func = row.get("Anchored_Known_System", "").strip()
            
            # 排除纯新系统（无已知防御基因）或者空值
            if known_func and known_func != "Pure_Novel_System":
                try:
                    freq = int(row.get("Frequency", 0))
                except ValueError:
                    freq = 0
                    
                # 累加出现频次
                freq_map[known_func] += freq
                
    # 按照频次从高到低排序
    sorted_combinations = sorted(freq_map.items(), key=lambda x: x[1], reverse=True)
    
    print(f"\n{'='*75}")
    print(f"🌟 {env_name} 中候选新系统最常 '抱大腿' 的已知防御机制 (Top 15)")
    print(f"{'='*75}")
    print(f"{'排名':<4} | {'锚定的已知防御系统 (家族=机制)':<48} | {'共现总频次'}")
    print("-" * 75)
    
    # 输出前 15 名
    for idx, (system_comb, freq) in enumerate(sorted_combinations[:15]):
        # 截断过长字符串以保持终端打印格式对齐
        sys_str = system_comb[:46] + ".." if len(system_comb) > 48 else system_comb
        print(f"Top{idx+1:<3}| {sys_str:<48} | {freq}")
        
    print("-" * 75)

if __name__ == "__main__":
    for env in ["TPMC-A", "TPMC-S"]:
        analyze_top_cooccurrence(env)