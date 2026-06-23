# 10_photoperiod_phase.R —— 峰值时刻的跨物种格局 + 光周期效应(诚实标注混杂)
# 输入: results/tables/02_<id>_rhythm.csv (节律基因 + 相位)
# 输出: results/tables/10_peak_timing.csv(各物种 后半天峰值占比 + 双峰性指标)
#       results/figures/10_phase_ridge.{pdf,png}(各物种相位密度, 按进化排序)
#       10_dusk_peak_fraction.{pdf,png}(后半天峰值占比 + 光周期标注)
# 注: 复现论文结论——从灰胞藻到裸子植物, 多数节律基因在后半天(黄昏/夜)达峰;
#     红藻 P. purpureum 与蓝藻 Synechocystis 例外, 呈双峰(光照中段+黑暗中段)。
#     光周期(12L/12D vs 16L/8D)与系统发育混杂, 只作描述性标注, 不做因果。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)

# 双峰性: 相位直方做 6 桶, 谷指数 = 1 - 中间最小桶/两侧峰均值(越大越双峰)
bimodality <- function(ph) {
  h <- hist(ph, breaks=seq(0,24,4), plot=FALSE)$counts
  if (sum(h)==0) return(NA)
  h <- h/sum(h)
  # 找两个最大桶, 看它们之间是否有谷
  pk <- order(-h)[1:2]
  if (abs(pk[1]-pk[2]) < 2) return(0)                 # 两峰相邻=单峰
  valley <- min(h[(min(pk)+1):(max(pk)-1)])
  round(1 - valley/mean(h[pk]), 2)
}

main <- function() {
  cat("== 10 峰值时刻格局 + 光周期 ==\n")
  rows <- list(); ph_all <- list()
  for (sp in SPECIES$id) {
    r <- read.csv(file.path(TAB_DIR, sprintf("02_%s_rhythm.csv", sp)), stringsAsFactors=FALSE)
    rh <- r[which(r$jtk_adjp < 0.05), ]
    ph <- rh$paper_phase; ph[is.na(ph)] <- rh$jtk_lag[is.na(ph)]; ph <- ph[!is.na(ph)]
    if (!length(ph)) next
    lh <- LIGHT_HOURS[[sp]]
    second_half <- mean(ph >= 12 & ph < 24) * 100      # 后半天(ZT12-24)峰值占比
    bim <- bimodality(ph)
    rows[[sp]] <- data.frame(id=sp, latin=SPECIES$latin[SPECIES$id==sp],
      clade=SPECIES$clade[SPECIES$id==sp], photoperiod=sprintf("%dL/%dD", lh, 24-lh),
      n_rhythmic=length(ph), second_half_frac=round(second_half,1), bimodality=bim)
    ph_all[[sp]] <- data.frame(id=sp, photoperiod=sprintf("%dL/%dD", lh, 24-lh), phase=ph)
    cat(sprintf("%-4s [%s] 后半天峰值 %.0f%%, 双峰指数 %.2f\n",
                sp, sprintf("%dL/%dD",lh,24-lh), second_half, bim))
  }
  tt <- do.call(rbind, rows); save_table(tt, "10_peak_timing")

  # ── 各物种相位密度(脊线图, 按进化排序) ──
  pa <- do.call(rbind, ph_all); pa$id <- factor(pa$id, levels=rev(SPECIES$id))
  labs_sp <- setNames(sprintf("%s (%s)", SPECIES$latin, c(Syn="12L",Cpa="16L",Ppu="12L",Cre="12L",
    Kni="12L",Ppa="12L",Smo="12L",Pab="16L",Osa="12L",Ath="12L")[SPECIES$id]), SPECIES$id)
  p1 <- ggplot(pa, aes(phase, id, fill=photoperiod)) +
    geom_violin(scale="width", linewidth=0.2, alpha=0.8) +
    annotate("rect", xmin=12, xmax=24, ymin=-Inf, ymax=Inf, fill="grey80", alpha=0.25) +
    scale_fill_manual(values=c("12L/12D"="#E67E22","16L/8D"="#34495E"), name="Photoperiod") +
    scale_y_discrete(labels=labs_sp) +
    labs(x="Peak phase (ZT, h), grey = 2nd half of day", y=NULL,
         title="Peak-time shifts toward dusk/night across the evolutionary gradient") +
    theme_pub_bw() + theme(axis.text.y=element_text(face="italic", size=6))
  savefig(p1, "10_phase_ridge", w=7, h=5)

  # ── 后半天峰值占比条形(标注光周期 + 双峰例外) ──
  tt$id <- factor(tt$id, levels=SPECIES$id)
  tt$label <- ifelse(tt$bimodality >= 0.5, "bimodal", "")
  p2 <- ggplot(tt, aes(id, second_half_frac, fill=photoperiod)) +
    geom_col(width=0.7) +
    geom_hline(yintercept=50, linetype="dashed", linewidth=LW, color="grey40") +
    geom_text(aes(label=label), vjust=-0.4, size=PT, fontface="italic", color="#B2182B") +
    scale_fill_manual(values=c("12L/12D"="#E67E22","16L/8D"="#34495E"), name="Photoperiod") +
    labs(x=NULL, y="Genes peaking in 2nd half of day (%)",
         title="Most species peak at dusk/night; red alga & cyanobacterium are bimodal exceptions") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))
  savefig(p2, "10_dusk_peak_fraction", w=6.5, h=3.8)

  # 汇总: 后半天峰值 + 双峰例外
  cat(sprintf("\n后半天(黄昏/夜)峰值 >50%% 的物种: %d/%d\n",
              sum(tt$second_half_frac>50), nrow(tt)))
  cat("双峰物种(指数≥0.5):", paste(tt$id[tt$bimodality>=0.5], collapse=", "), "\n")
  cat("== 10 完成 ==\n")
}

if (sys.nframe() == 0) main()
