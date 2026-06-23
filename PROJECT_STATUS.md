# PROJECT_STATUS — 09 Archaeplastida 昼夜节律转录组（E-MTAB-7188）

> 图内英文、代码注释中文、README/笔记中文。

## 一句话定位
9 类 Archaeplastida（蓝藻→藻类→陆生植物）昼夜时序转录组跨界比较：JTK_Cycle 节律检测 + 峰值相位分布 + 进化梯度。数据 = Ferrari 2019 Nat Commun 补充数据（10 物种 TPM 矩阵 + 论文 JTK 结果）。

## 流程位置：harness 8 阶段
- P0 INTAKE ✅　P1 SCAFFOLD ✅　P2 DATA ✅　P3 IMPLEMENT ✅　P4 RESULTS ✅
- P5 SHIP ⬜　P6 LEARN ⬜　P7 CLOSE ⬜（进行中）

## 主要结果（P4 已核）
- JTK 复现 9/10 物种 Jaccard≥0.97；水稻/拟南芥节律基因数与论文完全相同
- 单细胞藻 77% vs 多细胞陆生植物 42% 节律基因，Wilcoxon p=0.008（论文核心结论）
- 拟南芥 GO 首位 = circadian rhythm（p.adj 1.1e-8）；8 核心钟基因全命中，相位级联正确
- 云杉自算 2136 vs 论文 5280：低表达过滤所致，Jaccard=1.00（干净子集），已诚实标注
- **进阶（06-10，复现论文进化学结论）**：
  - 06 基因年龄×节律：老基因（Prokaryotic/Archaeplastida）显著富集节律，新基因(Specific)缺失，置换 p<0.05
  - 07 直系同源相位保守：39 物种对 28 对显著（Δphase<随机，p<0.05），Cpa-Cre=0.72 最佳（同论文 Fig3）
  - 08 细胞周期基因昼夜调控：单细胞 73% vs 多细胞 38%，Wilcoxon p=0.008（同摘要）
  - 09 WGCNA 模块昼夜性：节律模块比例单细胞>多细胞（Ppu/Cre 100% vs 云杉 27%），模块层面正交验证
  - 10 峰值时刻：9/10 物种黄昏/夜间峰；蓝藻 Syn(0.94)+红藻 Ppu(0.56) 双峰例外命中（同论文 Fig2）
- 交付：27 图（PDF+PNG）+ 29 表，run_all 零到出图可复现（00-10）

## 这是个什么项目（与前面项目的差异）
- 前 8 个项目多为两组对比（差异/诊断/预后）；本项目是**时序转录组**——核心是"基因表达随昼夜周期的节律"，方法是 JTK_Cycle 周期检测 + 相位分析，不是差异表达。
- **多物种比较**：9 类群覆盖蓝藻到陆生植物，论文核心结论是单细胞藻多数基因受昼夜调控、多细胞陆生植物组织层面不明显。
- 数据自带论文算好的 JTK pVal/phase → 自算结果可与论文对照（复现验证）。

## 数据计划（P2）
- 下载 Supplementary Data 2–11（MOESM5–14）共 10 个 xlsx → data/SupplementaryData2..11.xlsx
- 来源：Nature 补充材料 static-content.springer.com
- 用户已确认：全部下载
- 校验：每个 xlsx 基因数、表达列（L4/L8/L12/D4/D8/D12×重复）、pVal/phase 列存在

## 已建（P1 SCAFFOLD）
- scripts/：_common.R（9 物种注册表 + xlsx 加载 + JTK 对照常量）+ 00_install + 00–05 六 stub
- README.md 骨架（主要结果占位待 P4）
- .gitignore（白名单纳 xlsx，挡 fastq/大文件）
- LICENSE（MIT）
- 本 PROJECT_STATUS.md

## 下一步（P2 DATA）
下载 10 个 xlsx → 逐个校验结构 → 进 P3。

## 待确认映射（P2 已核验）
- 10 个 xlsx 全部下载并校验通过，pVal/phase 列齐全。
- **各物种结构差异大（P3 必须按列名解析时间点，不能假设固定网格）**：
  - 基因数 3467(Syn) ~ 36948(Pab)；表达列 12(Osa) / 18(Ath) / 24(Syn,Cre, 2重复) / 36(其余, 3重复)
  - 时间点编码在列名: 多数 L1..11/D1..11(2h间隔), Syn 是 L1.5..11.5, Osa/Ath 是 L4/8/12(4h间隔)
  - 光周期不同: 多数 12L/12D; Cpa/Ppa/Pab 是 16L/8D(L1..15 D1..7)
  - **重复后缀两种**: 多数 `.1/.2`, Selaginella(Data8) 是 `...A/B/C`
  - 列名解析器: `^[LD]\d+(\.5)?` 抓 L/D + ZT小时(可小数), 已验证 10 物种全部解析成功
  - ZT 换算: D 期 ZT = num + 光照时长(12 或 16)
