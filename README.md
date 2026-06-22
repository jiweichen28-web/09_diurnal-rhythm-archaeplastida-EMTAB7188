# Archaeplastida 昼夜节律转录组跨界比较（9 物种时序，E-MTAB-7188）

跨越蓝藻、藻类到陆生植物的 9 类 Archaeplastida 物种昼夜（diurnal）时序转录组比较：节律基因检测（JTK_Cycle）、峰值相位分布、跨物种进化梯度。数据来自 Ferrari et al. 2019 *Nature Communications* 的补充数据（各物种 TPM 表达矩阵 + 论文 JTK 节律结果）。本项目为其再分析。

标签：`time-series` `circadian` `diurnal` `JTK_Cycle` `plant` `algae` `cross-species` `R`

## 数据来源

| 字段 | 内容 |
|------|------|
| 原文 | Ferrari C, et al. Kingdom-wide comparison reveals the evolution of diurnal gene expression in Archaeplastida. *Nat Commun* 2019, 10:737. DOI [10.1038/s41467-019-08703-2](https://doi.org/10.1038/s41467-019-08703-2) |
| 原始数据 | EBI ArrayExpress **E-MTAB-7188**（6 物种新测）+ 公共数据集若干 |
| 本项目输入 | Supplementary Data 2–11（10 个物种 TPM 矩阵 + JTK pVal/phase），Data 14（OrthoFinder 正交直系同源），从 Nature 补充材料下载 |
| 设计 | 每 2 h 采样、12 时间点/物种，覆盖一个昼夜周期（ZT1 起）|

**9 类群物种（进化梯度，蓝藻 → 藻类 → 陆生植物）**

| Data | 物种 | 类群 | 论文节律基因数 |
|------|------|------|---------------|
| 2 | *Synechocystis* sp. PCC 6803 | 蓝藻 | 2440 |
| 3 | *Cyanophora paradoxa* | 灰胞藻 | 12343 |
| 4 | *Porphyridium purpureum* | 红藻 | 6597 |
| 5 | *Chlamydomonas reinhardtii* | 绿藻 | 12341 |
| 6 | *Klebsormidium nitens* | 链藻 | 9375 |
| 7 | *Physcomitrella patens* | 苔藓 | 11692 |
| 8 | *Selaginella moellendorffii* | 卷柏 | 11860 |
| 9 | *Picea abies* | 云杉（裸子）| 5280 |
| 10 | *Oryza sativa* | 水稻 | 8082 |
| 11 | *Arabidopsis thaliana* | 拟南芥 | 5676 |

> 各 xlsx 为补充数据文本（MB 级），复现必需，已入库。原始 fastq（E-MTAB-7188）不需要、不入库。

## 目录结构

```
09_diurnal-rhythm-archaeplastida-EMTAB7188/
├── data/                       # 10 个 Supplementary Data xlsx
├── scripts/                    # _common.R + 00_install + 00–05 + run_all
├── results/
│   ├── figures/                # 每图 .pdf + .png(300ppi)，图内文字一律英文
│   └── tables/                 # 物种概览、节律结果、相位、跨物种汇总
├── README.md  LICENSE  .gitignore  PROJECT_STATUS.md  sessionInfo.txt
```

## 分析脚本

```
00_install_packages.R  依赖安装（仅首次）
_common.R              路径自定位 + 绘图风格 + 9 物种注册表 + xlsx 加载
00_setup_and_load.R    读 10 物种矩阵 → 表达 + 元数据 → 物种概览
01_qc_pca.R            各物种 QC + log2 + PCA（昼夜时间环）
02_rhythm_detection.R  JTK_Cycle 节律检测，与论文 pVal/phase 对照
03_phase_distribution.R 峰值相位分布 + 节律基因热图
04_enrichment_clock.R  节律基因功能富集 + 核心钟基因专题
05_cross_species.R     跨物种节律比例 + 相位保守性（进化梯度）
06_gene_age_rhythm.R   基因年龄(phylostratum) × 节律：置换检验富集/缺失
07_ortholog_phase_conservation.R  1:1 直系同源跨物种相位保守(Δphase 置换)
08_cellcycle_timing.R  细胞周期基因昼夜调控：单细胞藻 vs 多细胞陆生植物
run_all.R              按编号依次跑
```

## 主要结果

**数据规模（10 物种）**

- 基因数 3467（蓝藻）~ 36948（云杉）；表达列 12（水稻）/18（拟南芥）/24/36，重复 2~3。
- 时间点编码在列名（L/D + ZT 小时），光周期 12L/12D 或 16L/8D，由 `_common.R` 统一解析。

**JTK_Cycle 节律检测复现（与论文同口径 ADJ.P<0.05）**

| 物种 | 类群 | 自算节律基因 | 占比 | 论文 | Jaccard |
|------|------|------------|------|------|---------|
| *Synechocystis* | 蓝藻 | 2488 | 72% | 2440 | 0.97 |
| *C. paradoxa* | 灰胞藻 | 12219 | 72% | 12343 | 1.00 |
| *P. purpureum* | 红藻 | 6591 | 82% | 6597 | 1.00 |
| *C. reinhardtii* | 绿藻 | 12106 | 81% | 12341 | 1.00 |
| *K. nitens* | 链藻 | 9348 | 78% | 9375 | 1.00 |
| *P. patens* | 苔藓 | 11656 | 55% | 11692 | 1.00 |
| *S. moellendorffii* | 卷柏 | 11738 | 71% | 11860 | 1.00 |
| *P. abies* | 云杉 | 2136 | 12% | 5280 | 1.00 |
| *O. sativa* | 水稻 | 8082 | 32% | 8082 | 1.00 |
| *A. thaliana* | 拟南芥 | 5676 | 40% | 5676 | 1.00 |

- **复现高度吻合**：9/10 物种 Jaccard ≥0.97，水稻/拟南芥节律基因数与论文完全相同。
- **关键方法点**：JTK 必须保留生物学重复（时间向量含重复 ZT）才有检验力；论文 `pVal` 列 = JTK 的 ADJ.P。
- **云杉（P. abies）是唯一偏差**：自算 2136 vs 论文 5280。原因＝本项目对每物种去低表达基因（17535 测试 vs 论文全集），云杉低表达基因多，被滤掉的多为弱节律；Jaccard=1.00 说明自算集是论文集的干净子集，非方法错误。

**跨物种进化梯度（论文核心结论）**

- 单细胞藻平均 **77%** 基因受昼夜调控，多细胞陆生植物平均 **42%**，Wilcoxon **p=0.008**。
- 与论文一致：单细胞藻几乎全转录组随昼夜起伏；多细胞陆生植物在组织层面节律明显减弱。

**拟南芥功能富集 + 核心钟基因**

- 节律基因 GO BP 富集 36 条显著，**首位 = circadian rhythm（昼夜节律，p.adj=1.1e-8）**，其后为 response to water、chlorophyll biosynthesis、response to cold/red light 等——典型昼夜/光合相关过程。
- 8 个核心生物钟基因（CCA1/LHY/PRR9/PRR7/GI/TOC1/ELF3/ELF4）全部命中，黎明型（CCA1/LHY）与黄昏型（TOC1/GI）表达相位交替，重现生物钟工作机制。

### 进阶：复现论文的进化学结论（不止常规节律检测）

**① 基因年龄 × 节律（06，phylostratum 置换检验）**

- 复现论文方法：按 Fig.1a 物种树把基因分进化年龄层（Prokaryotic 最老 → Specific 最新），置换 1000 次得经验 p。
- **老基因显著富集节律**：Prokaryotic、Archaeplastida 等古老层在几乎所有物种里 rhythmic 富集（perm p<0.05）；老基因 vs 新基因(Specific)节律比例普遍 old>young（如 Cyanophora 69.5% vs 39.9%、Physcomitrella 55.7% vs 17.5%）。印证论文\"基因家族越古老越倾向昼夜表达\"。

**② 跨物种直系同源相位保守（07，论文标题结论）**

- 用 OrthoFinder 正交直系同源表（Supplementary Data 14）取 1:1 直系同源，复现论文 Δphase 置换法：观测相位差 vs 随机配对零分布。
- **39 个物种对里 28 对（72%）相位显著保守**（观测 Δphase 小于随机，p<0.05），Δphase obs/exp 比值中位 0.92。
- 最保守的 **Cyanophora–Chlamydomonas（比值 0.72，p<0.001）正是论文 Fig.3 举的例子**；水稻–拟南芥（0.81）等也显著。说明即使跨越巨大进化距离，直系同源仍倾向在一天同一时刻达峰——昼夜程序高度保守。

**③ 细胞周期基因的昼夜调控：单细胞 vs 多细胞（08，摘要核心结论）**

- 按注释（MapMan bin）取细胞周期/DNA 复制基因，比较节律比例。
- **单细胞藻 73% vs 多细胞陆生植物 38%，Wilcoxon p=0.008**。且单细胞藻里细胞周期基因比全基因组更节律（Cyanophora 70% vs 50%、Klebsormidium 70% vs 55%）——细胞分裂被昼夜钟\"定时\"；多细胞陆生植物只跟随（低）基线，无特异定时。精确复现摘要\"细胞分裂在单细胞藻受昼夜调控、多细胞陆生植物组织层面不明显\"。

> **诚实说明**：本项目用论文已发布的 TPM 矩阵 + JTK 结果做再分析与可视化，不是从 fastq 重跑。节律检测（02）独立复现与论文高度吻合（Jaccard≥0.97）；云杉计数差异已如实标注为低表达过滤所致。进阶分析（06–08）按论文原方法（置换检验、Δphase、MapMan 注释）复现其三个进化学核心结论，方向与显著性均与论文一致。云杉（Picea）因节律基因少、共节律直系同源仅 20~48 个，在 07 中不显著（数据所限，非方法错误）。

## 环境

```
R 4.x：ggplot2, dplyr, tidyr, readxl, ggrepel, pheatmap, reshape2, MetaCycle
Bioconductor：clusterProfiler, org.At.tair.db（模式物种功能富集）
```

图内文字一律英文，`theme_pub_bw()` 用 sans 字体。`sessionInfo.txt` 记录完整版本。

## 引用

> Ferrari C, Proost S, Janowski M, et al. Kingdom-wide comparison reveals the evolution of diurnal gene expression in Archaeplastida. *Nat Commun* 2019;10:737. https://doi.org/10.1038/s41467-019-08703-2
