# DefenseProphet-Miner

本项目是用于**微生物新型防御系统的高通量挖掘与生态演化网络解析**的计算分析工作流代码库。

本研究聚焦于突破传统的序列同源比对局限，整合了深度学习模型 DefenseProphet 与多源宏基因组数据集（包括公共基准数据库 GTDB 与青藏高原极端生境数据集），致力于鉴定结构未知、演化机制全新的抗噬菌体免疫组件。

---

## ⚠️ 数据文件说明 (Data Availability)

**请注意：** 
由于 GitHub 对单个文件的上传大小有 100MB 的限制，为了保证代码库的正常托管，本仓库仅移除了以下三个体积过大的中间态特征矩阵文件：
- `bacteria_comparison_results_high_conf.csv`
- `TPMC-A_all_defense.csv`
- `multi_label_other_bacteria.csv`

除上述三个超大文件外，本项目的核心代码流、分析管线以及其余轻量级数据文件均完整保留。

---

## 核心分析管线 (Core Pipeline)

本代码库主要包含以下四个串联的推断与统计清洗步骤：

1. **候选基因的高置信度筛选 (Candidate Screening)**
   - 基于二分类打分（>0.99）与多分类防御亚型判定（Other），通过向量化矩阵操作快速清洗低置信度基线序列，界定出强潜力防御靶标。
2. **物理空间成簇分析 (Spatial Clustering)**
   - 依据“防御岛”理论，执行基于 1500 bp 一维物理间距的一维成簇算法。将邻近节点划分为纯新系统簇（PURE_CANDIDATE）与混合嵌合簇（HYBRID）。
3. **远源特征去冗余聚类 (Homology Clustering removal)**
   - 调度 CD-HIT 算法以极低的同源性阈值（40%）进行聚类分析，剔除重复片段并提取具有代表性的同源核心家族。
4. **复合网络拓扑解析 (Topological Network Analysis)**
   - 构建独创的拓扑表征规律（如 `K-U#` 结构表达式），转化复杂且异质的防御复合体构成，进而解构三种微生物免疫元件的重组演化模式（防线嵌合、内部扩增、毒性重塑）。

## 依赖环境 (Requirements)

本数据清洗分析流推荐在以下环境下运行：
- 操作系统：CentOS 7 / macOS
- 编程语言：Python 3.11.7
- 核心依赖库：`pandas`, `numpy`, `biopython`, `torch`, `networkx`
- 外部生信软件：`CD-HIT`
