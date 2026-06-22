# 04_enrichment_clock.R —— 节律基因 GO 富集(拟南芥) + 核心钟基因表达曲线
# 输入: results/tables/02_Ath_rhythm.csv + data/SupplementaryData11.xlsx
# 输出: results/tables/04_Ath_GO.csv, results/figures/04_Ath_GO.{pdf,png}
#       04_clock_genes.{pdf,png}(CCA1/LHY/TOC1/PRR7/PRR9/GI/ELF3/ELF4 昼夜表达)
# 注: 拟南芥 TAIR 位点有 GO 注释可做真富集; 核心钟基因展示黎明/黄昏交替表达(钟的工作原理)。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)

# 拟南芥核心生物钟基因(TAIR 位点 -> 常用名)
CLOCK <- c(AT2G46830="CCA1", AT1G01060="LHY", AT5G61380="TOC1", AT5G02810="PRR7",
           AT2G46790="PRR9", AT1G22770="GI", AT2G25930="ELF3", AT2G40080="ELF4")

main <- function() {
  cat("== 04 富集 + 钟基因 ==\n")
  d <- load_species("Ath")
  r <- read.csv(file.path(TAB_DIR, "02_Ath_rhythm.csv"), stringsAsFactors=FALSE)

  # ── 核心钟基因昼夜表达曲线(时间点均值) ──
  zt <- d$design$zt; uz <- sort(unique(zt))
  Xm <- sapply(uz, function(z) rowMeans(d$expr[, zt==z, drop=FALSE], na.rm=TRUE))
  rownames(Xm) <- toupper(rownames(d$expr)); colnames(Xm) <- uz
  ck <- intersect(names(CLOCK), rownames(Xm))
  long <- do.call(rbind, lapply(ck, function(g) {
    v <- Xm[g, ]; v <- (v - min(v)) / (max(v) - min(v) + 1e-9)   # 0-1 归一便于比较
    data.frame(gene=CLOCK[[g]], zt=as.numeric(colnames(Xm)), expr=v)
  }))
  long$gene <- factor(long$gene, levels=c("CCA1","LHY","PRR9","PRR7","GI","TOC1","ELF3","ELF4"))
  p_ck <- ggplot(long, aes(zt, expr, color=gene)) +
    annotate("rect", xmin=12, xmax=max(uz), ymin=-Inf, ymax=Inf, fill="grey85", alpha=0.5) +
    geom_line(linewidth=0.5) + geom_point(size=1) +
    scale_color_brewer(palette="Set1", name="Clock gene") +
    scale_x_continuous(breaks=uz) +
    labs(x="ZT (h), grey = dark", y="Scaled expression (0-1)",
         title="Arabidopsis core clock genes over the diurnal cycle") + theme_pub_bw()
  savefig(p_ck, "04_clock_genes", w=6, h=3.8)
  cat(sprintf("   核心钟基因 %d/%d 命中\n", length(ck), length(CLOCK)))

  # ── GO 富集(节律基因 vs 全部测试基因为背景) ──
  if (requireNamespace("clusterProfiler", quietly=TRUE) &&
      requireNamespace("org.At.tair.db", quietly=TRUE)) {
    suppressMessages({library(clusterProfiler); library(org.At.tair.db)})
    # TAIR 位点需大写(数据里是 At2g.. 混合大小写, org.At.tair.db 要 AT2G..)
    to_tair <- function(g) toupper(g[grepl("^AT[0-9CM]G[0-9]+$", g, ignore.case=TRUE)])
    rh_genes <- to_tair(r$GeneName[which(r$jtk_adjp < 0.05)])
    bg_genes <- to_tair(r$GeneName)
    cat(sprintf("   TAIR 位点: 节律 %d / 背景 %d\n", length(rh_genes), length(bg_genes)))
    ego <- tryCatch(enrichGO(gene=rh_genes, universe=bg_genes, OrgDb=org.At.tair.db,
             keyType="TAIR", ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1),
             error=function(e){cat("GO 富集出错:", conditionMessage(e), "\n"); NULL})
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      go_df <- as.data.frame(ego)
      save_table(go_df[, c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","Count")], "04_Ath_GO")
      top <- head(go_df[order(go_df$p.adjust), ], 15)
      top$Description <- factor(top$Description, levels=rev(top$Description))
      p_go <- ggplot(top, aes(-log10(p.adjust), Description, fill=Count)) +
        geom_col(width=0.7) + scale_fill_gradient(low="#9ECAE1", high="#08519C", name="Count") +
        labs(x="-log10 adjusted P", y=NULL,
             title="GO BP enrichment of Arabidopsis rhythmic genes") +
        theme_pub_bw() + theme(axis.text.y=element_text(size=6))
      savefig(p_go, "04_Ath_GO", w=7, h=4.5)
      cat(sprintf("   GO BP 富集 %d 条显著\n", nrow(go_df)))
    } else cat("   GO 富集无显著结果或出错, 跳过\n")
  } else cat("   org.At.tair.db 未装, 跳过 GO 富集\n")
  cat("== 04 完成 ==\n")
}

if (sys.nframe() == 0) main()
