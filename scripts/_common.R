# _common.R —— 项目共享配置与函数 (R 版)
# 植物/藻类昼夜节律转录组跨界比较 (Ferrari et al. 2019, Nat Commun; E-MTAB-7188)
# 9 类 Archaeplastida 物种, 每物种昼夜时序表达矩阵 + 论文 JTK_Cycle 节律结果
# 固定绘图规范: 8pt 字体, theme_pub_bw, 每图 PDF + 300ppi PNG, 图内文字一律英文, 注释中文
# 路径全部相对项目根; 本文件自动定位项目根并 setwd, 从任意目录运行均可。

suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(readxl); library(pheatmap)})

# ── 自动定位项目根: 本脚本路径 -> 上溯一级(scripts -> 根) ──
.find_proj_root <- function() {
  a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 0) for (i in sys.nframe():1) { of <- sys.frame(i)$ofile; if (!is.null(of)) { f <- of; break } }
  if (length(f) == 0 || is.na(f[1]) || f[1] == "") return(NULL)
  normalizePath(file.path(dirname(normalizePath(f[1], winslash="/", mustWork=FALSE)), ".."), winslash="/", mustWork=FALSE)
}
.root <- .find_proj_root()
if (!is.null(.root) && dir.exists(file.path(.root, "data"))) setwd(.root)

# ── 输出目录 ──
FIG_DIR <- "results/figures"; TAB_DIR <- "results/tables"
for (d in c(FIG_DIR, TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive=TRUE)

# ── 数据集常量 ──
ACCESSION <- "E-MTAB-7188"
PAPER     <- "Ferrari et al. 2019, Nat Commun 10:737"

# 9 物种注册表: 短名 / 拉丁名 / 类群 / 对应 Supplementary Data 文件
# (Data 2-11 = 10 个矩阵; P. patens 同时有 12L/12D 与 16L/8D 两套, 论文按 9 类群讨论)
SPECIES <- data.frame(
  id    = c("Syn","Cpa","Ppu","Cre","Kni","Ppa","Smo","Pab","Osa","Ath"),
  latin = c("Synechocystis sp. PCC 6803","Cyanophora paradoxa","Porphyridium purpureum",
            "Chlamydomonas reinhardtii","Klebsormidium nitens","Physcomitrella patens",
            "Selaginella moellendorffii","Picea abies","Oryza sativa","Arabidopsis thaliana"),
  clade = c("Cyanobacteria","Glaucophyte","Rhodophyte","Chlorophyte(green alga)",
            "Charophyte","Moss","Lycophyte","Gymnosperm(spruce)","Monocot(rice)","Eudicot(Arabidopsis)"),
  data_file = c("SupplementaryData2.xlsx","SupplementaryData3.xlsx","SupplementaryData4.xlsx",
                "SupplementaryData5.xlsx","SupplementaryData6.xlsx","SupplementaryData7.xlsx",
                "SupplementaryData8.xlsx","SupplementaryData9.xlsx","SupplementaryData10.xlsx",
                "SupplementaryData11.xlsx"),
  stringsAsFactors = FALSE)

# 论文 JTK 报告的节律基因数(校验/对照用)
PAPER_RHYTHMIC <- c(Syn=2440, Cpa=12343, Ppu=6597, Cre=12341, Kni=9375,
                    Ppa=11692, Smo=11860, Pab=5280, Osa=8082, Ath=5676)

# 各物种光照时长(h): 多数 12L/12D; Cpa/Ppa/Pab 是 16L/8D
LIGHT_HOURS <- c(Syn=12, Cpa=16, Ppu=12, Cre=12, Kni=12, Ppa=16, Smo=12, Pab=16, Osa=12, Ath=12)

# ── 绘图风格 (8pt, 黑色描边, 无网格; 图内文字一律英文) ──
FONT_SIZE <- 8; PT <- FONT_SIZE / 2.845; LW <- 0.5 / 2.1333; PLOT_FAMILY <- "sans"
# 类群配色(蓝->绿->深绿, 进化梯度: 蓝藻->藻类->陆生植物)
CLADE_COL <- c(Syn="#2C3E80", Cpa="#3A5FCD", Ppu="#9B59B6", Cre="#1ABC9C", Kni="#16A085",
               Ppa="#27AE60", Smo="#2ECC71", Pab="#7F8C3A", Osa="#E67E22", Ath="#C0392B")
# 昼夜相位配色环(给相位图用)
PHASE_COL <- colorRampPalette(c("#FDE725","#21908C","#440154","#21908C","#FDE725"))(24)

theme_pub_bw <- function() {
  theme_bw(base_size=FONT_SIZE) +
    theme(text=element_text(size=FONT_SIZE, color="black", family=PLOT_FAMILY),
      plot.title=element_text(size=FONT_SIZE, hjust=0.5), axis.title=element_text(size=FONT_SIZE),
      axis.text=element_text(size=FONT_SIZE, color="black"), legend.title=element_text(size=FONT_SIZE),
      legend.text=element_text(size=FONT_SIZE), strip.text=element_text(size=FONT_SIZE),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
      panel.border=element_rect(color="black", fill=NA, linewidth=LW), axis.line=element_blank(),
      axis.ticks=element_line(color="black", linewidth=LW),
      panel.background=element_rect(fill="transparent", color=NA),
      plot.background=element_rect(fill="transparent", color=NA),
      legend.background=element_rect(fill="transparent", color=NA),
      legend.key=element_rect(fill="transparent", color=NA))
}

savefig <- function(p, name, w=5, h=4) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width=w, height=h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width=w, height=h, dpi=300)
  cat(sprintf("   [fig] %s.pdf / .png\n", name))
}
save_table <- function(df, name, row.names=FALSE) {
  write.csv(df, file.path(TAB_DIR, paste0(name, ".csv")), row.names=row.names, fileEncoding="UTF-8")
  cat(sprintf("   [tab] %s.csv  (%d x %d)\n", name, nrow(df), ncol(df)))
}

# ── 数据加载 ──
# 每个 Supplementary Data xlsx: 前3行说明, 第4行表头
#   元数据列: GeneName, Phylostratum, pVal(JTK p), phase(峰值相位1-24), Annotation
#   表达列: L4/L8/L12/D4/D8/D12 各若干重复(列名带 .1/.2 后缀)
.META_COLS <- c("GeneName","Phylostratum","pVal","phase","Annotation")

# 读单物种: 返回 list(meta=元数据df, expr=表达矩阵 基因x样本, design=样本设计df)
# design: 每个表达列 -> ld(L/D), zt(连续昼夜小时 0-24), tp(时间点标签)
load_species <- function(sp_id) {
  f <- file.path("data", SPECIES$data_file[SPECIES$id == sp_id])
  raw <- as.data.frame(read_excel(f, skip=3, col_names=TRUE))
  raw <- raw[!is.na(raw$GeneName) & raw$GeneName != "", ]
  expr_cols <- setdiff(colnames(raw), .META_COLS)
  expr <- as.matrix(sapply(raw[, expr_cols, drop=FALSE], function(z) suppressWarnings(as.numeric(z))))
  rownames(expr) <- raw$GeneName
  meta <- raw[, intersect(.META_COLS, colnames(raw)), drop=FALSE]
  meta$pVal  <- suppressWarnings(as.numeric(meta$pVal))    # "NR" -> NA
  meta$phase <- suppressWarnings(as.numeric(meta$phase))
  list(meta=meta, expr=expr, design=parse_design(expr_cols, sp_id))
}

# 解析表达列名 -> 样本设计. 列名形如 L<num>[.rep] 或 L<num>...<A/B/C>; num 可小数(1.5)
# zt = 连续昼夜小时: 光照期 L 直接用 num; 黑暗期 D 用 num + 光照时长
parse_design <- function(cols, sp_id) {
  lh <- LIGHT_HOURS[[sp_id]]
  m <- regmatches(cols, regexec("^([LD])([0-9]+(?:\\.5)?)", cols))
  ld  <- sapply(m, function(x) if (length(x)==3) x[2] else NA)
  num <- sapply(m, function(x) if (length(x)==3) as.numeric(x[3]) else NA)
  zt  <- ifelse(ld=="D", num + lh, num)
  data.frame(col=cols, ld=ld, num=num, zt=zt, tp=paste0(ld, num), stringsAsFactors=FALSE)
}

# BH-FDR
bh_fdr <- function(p) p.adjust(p, method="BH")
