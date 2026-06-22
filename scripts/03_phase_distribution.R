# 03_phase_distribution.R —— 峰值相位分布 + 节律基因热图(按相位排序)
# 输入: results/tables/02_<id>_rhythm.csv + data/*.xlsx(表达谱)
# 输出: results/figures/03_phase_all_species.{pdf,png}(各物种相位直方汇总)
#       03_<id>_rhythm_heatmap.{pdf,png}(模式物种节律基因相位热图)
#       results/tables/03_phase_summary.csv
# 注: 相位(phase 1-24)=峰值时刻; 多数物种黎明/黄昏富集, 是昼夜程序核心特征。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)

main <- function() {
  cat("== 03 相位分布 ==\n")
  phase_all <- list()
  for (sp in SPECIES$id) {
    r <- read.csv(file.path(TAB_DIR, sprintf("02_%s_rhythm.csv", sp)), stringsAsFactors=FALSE)
    rh <- r[which(r$jtk_adjp < 0.05), ]
    if (nrow(rh) == 0) next
    # 用论文 phase 列(已验证, 1-24); 缺失则退到自算 LAG
    ph <- rh$paper_phase; ph[is.na(ph)] <- rh$jtk_lag[is.na(ph)]
    phase_all[[sp]] <- data.frame(species=sp, phase=ph[!is.na(ph)])
    cat(sprintf("%-4s 节律基因 %d, 有相位 %d\n", sp, nrow(rh), sum(!is.na(ph))))
  }
  pa <- do.call(rbind, phase_all)
  pa$species <- factor(pa$species, levels=SPECIES$id)
  labs_sp <- setNames(SPECIES$latin, SPECIES$id)

  # 各物种相位直方(0-24)
  p <- ggplot(pa, aes(phase, fill=species)) +
    geom_histogram(binwidth=2, color="white", linewidth=0.1) +
    scale_fill_manual(values=CLADE_COL, guide="none") +
    facet_wrap(~species, scales="free_y", ncol=4, labeller=labeller(species=labs_sp)) +
    labs(x="Peak phase (ZT, h)", y="Rhythmic genes",
         title="Peak-phase distribution of rhythmic genes across species") +
    theme_pub_bw() + theme(strip.text=element_text(face="italic", size=6))
  savefig(p, "03_phase_all_species", w=9, h=6)

  # 相位汇总表: 各物种黎明峰(phase 0-6) / 黄昏峰(10-14) 占比
  summ <- do.call(rbind, lapply(split(pa, pa$species), function(s) {
    if(nrow(s)==0) return(NULL)
    data.frame(id=s$species[1], n=nrow(s),
      dawn_frac=round(mean(s$phase>=0 & s$phase<6)*100,1),
      dusk_frac=round(mean(s$phase>=10 & s$phase<16)*100,1))
  }))
  save_table(summ, "03_phase_summary")

  # 模式物种节律基因热图(拟南芥): 表达 z-score, 行按相位排序
  for (sp in c("Ath","Cre")) {
    d <- load_species(sp)
    r <- read.csv(file.path(TAB_DIR, sprintf("02_%s_rhythm.csv", sp)), stringsAsFactors=FALSE)
    rh <- r[which(r$jtk_adjp < 0.05 & !is.na(r$paper_phase)), ]
    rh <- rh[order(rh$paper_phase), ]
    # 时间点均值矩阵(按 zt 排序)
    zt <- d$design$zt; uz <- sort(unique(zt))
    Xm <- sapply(uz, function(z) rowMeans(d$expr[, zt==z, drop=FALSE], na.rm=TRUE))
    rownames(Xm) <- rownames(d$expr); colnames(Xm) <- paste0("ZT", uz)
    M <- Xm[rh$GeneName, , drop=FALSE]; M <- log2(M + 1)
    Mz <- t(scale(t(M))); Mz[is.na(Mz)] <- 0; Mz[Mz>2.5] <- 2.5; Mz[Mz< -2.5] <- -2.5
    for (ext in c("png","pdf")) {
      if (ext=="png") png(file.path(FIG_DIR, sprintf("03_%s_rhythm_heatmap.png", sp)), width=4, height=5, units="in", res=300)
      else pdf(file.path(FIG_DIR, sprintf("03_%s_rhythm_heatmap.pdf", sp)), width=4, height=5)
      pheatmap(Mz, cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=FALSE,
               color=colorRampPalette(c("#3B528B","white","#FDE725"))(50), fontsize=8,
               main=sprintf("%s rhythmic genes (sorted by phase)", SPECIES$latin[SPECIES$id==sp]))
      dev.off()
    }
    cat(sprintf("   [fig] 03_%s_rhythm_heatmap (%d 基因)\n", sp, nrow(rh)))
  }
  cat("== 03 完成 ==\n")
}

if (sys.nframe() == 0) main()
