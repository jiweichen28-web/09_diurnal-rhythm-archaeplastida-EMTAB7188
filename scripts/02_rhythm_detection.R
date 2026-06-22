# 02_rhythm_detection.R —— JTK_CYCLE 节律检测(MetaCycle), 复现并对照论文 JTK 结果
# 输入: data/*.xlsx
# 输出: results/tables/02_<id>_rhythm.csv(自算 JTK ADJ.P/BH.Q/LAG + 论文 pVal/phase)
#       results/tables/02_rhythmic_counts.csv
#       results/figures/02_rhythmic_fraction.{pdf,png} 02_jtk_concordance.{pdf,png}
# 注: 关键——保留生物学重复(时间向量含重复 ZT), JTK 才有足够检验力;
#     论文 pVal 列 = JTK 的 ADJ.P, 自算 ADJ.P<0.05 与论文节律基因数完全吻合。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R"))
suppressMessages(library(MetaCycle)); set.seed(42)

# 对单物种跑 JTK: 保留重复, 时间向量含重复 ZT
run_jtk <- function(d, sp) {
  keep <- rowSums(is.na(d$expr)) == 0 & rowSums(d$expr > 1, na.rm=TRUE) >= 3
  X <- d$expr[keep, , drop=FALSE]; zt <- d$design$zt
  ord <- order(zt); Xo <- X[, ord, drop=FALSE]; zto <- zt[ord]
  tmp <- file.path(TAB_DIR, sprintf("_tmp_%s.csv", sp))
  write.csv(data.frame(GeneName=rownames(Xo), Xo, check.names=FALSE), tmp, row.names=FALSE)
  res <- meta2d(infile=tmp, filestyle="csv", outdir=TAB_DIR, timepoints=zto,
                cycMethod="JTK", minper=24, maxper=24, outputFile=FALSE, releaseNote=FALSE)$JTK
  unlink(tmp)
  # 合并论文 pVal/phase 对照
  pm <- d$meta[match(res$CycID, d$meta$GeneName), ]
  data.frame(GeneName=res$CycID, jtk_adjp=res$ADJ.P, jtk_q=res$BH.Q,
             jtk_lag=res$LAG, jtk_amp=res$AMP, jtk_per=res$PER,
             paper_pval=pm$pVal, paper_phase=pm$phase, stringsAsFactors=FALSE)
}

main <- function() {
  cat("== 02 JTK 节律检测 ==\n")
  counts <- list()
  for (sp in SPECIES$id) {
    d <- load_species(sp); r <- run_jtk(d, sp)
    save_table(r, sprintf("02_%s_rhythm", sp))
    self_n  <- sum(r$jtk_adjp < 0.05, na.rm=TRUE)    # 与论文同口径(JTK ADJ.P)
    paper_n <- PAPER_RHYTHMIC[[sp]]
    frac <- self_n / nrow(r) * 100
    # 一致性: 自算 vs 论文 节律基因集重合
    self_set  <- r$GeneName[which(r$jtk_adjp < 0.05)]
    paper_set <- r$GeneName[which(r$paper_pval < 0.05)]
    jac <- length(intersect(self_set,paper_set)) / length(union(self_set,paper_set))
    counts[[sp]] <- data.frame(id=sp, latin=SPECIES$latin[SPECIES$id==sp],
      clade=SPECIES$clade[SPECIES$id==sp], tested=nrow(r), rhythmic_self=self_n,
      rhythmic_frac=round(frac,1), rhythmic_paper=paper_n, jaccard_vs_paper=round(jac,3))
    cat(sprintf("%-4s 自算 %5d (%.0f%%) 论文 %5d  Jaccard %.2f\n", sp, self_n, frac, paper_n, jac))
  }
  cc <- do.call(rbind, counts); save_table(cc, "02_rhythmic_counts")
  cc$id <- factor(cc$id, levels=SPECIES$id)

  p1 <- ggplot(cc, aes(id, rhythmic_frac, fill=id)) +
    geom_col(width=0.7) + scale_fill_manual(values=CLADE_COL, guide="none") +
    geom_text(aes(label=sprintf("%.0f%%", rhythmic_frac)), vjust=-0.3, size=PT) +
    labs(x=NULL, y="Rhythmic genes (%)", title="Fraction of rhythmic genes (JTK ADJ.P<0.05)") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
  savefig(p1, "02_rhythmic_fraction", w=6, h=3.5)

  p2 <- ggplot(cc, aes(rhythmic_paper, rhythmic_self, color=id)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=LW, color="grey50") +
    geom_point(size=2.5) + ggrepel::geom_text_repel(aes(label=id), size=PT) +
    scale_color_manual(values=CLADE_COL, guide="none") +
    labs(x="Paper rhythmic genes", y="This study (JTK ADJ.P<0.05)",
         title="Rhythmic gene count: reproduction vs paper") + theme_pub_bw()
  savefig(p2, "02_jtk_concordance", w=4.5, h=4.5)
  cat("== 02 完成 ==\n")
}

if (sys.nframe() == 0) main()
