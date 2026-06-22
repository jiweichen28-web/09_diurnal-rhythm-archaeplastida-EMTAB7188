# 00_setup_and_load.R —— 读 10 物种 xlsx → 表达矩阵 + 元数据 → 物种概览
# 输入: data/SupplementaryData2..11.xlsx
# 输出: results/tables/00_species_overview.csv
# 注: 各物种结构差异大(时间点/重复/光周期不同), 由 _common 的 parse_design 统一解析。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R"))

main <- function() {
  cat("== 00 读取 10 物种 ==\n")
  ov <- do.call(rbind, lapply(SPECIES$id, function(sp) {
    d <- load_species(sp)
    n_tp  <- length(unique(d$design$tp))
    n_rep <- round(ncol(d$expr) / n_tp, 1)
    paper_rh <- PAPER_RHYTHMIC[[sp]]
    self_rh  <- sum(d$meta$pVal < 0.05, na.rm=TRUE)   # 论文给的 pVal<0.05 计数
    cat(sprintf("%-4s %-28s 基因%6d 样本%3d 时间点%3d 重复%.1f\n",
                sp, SPECIES$latin[SPECIES$id==sp], nrow(d$expr), ncol(d$expr), n_tp, n_rep))
    data.frame(id=sp, latin=SPECIES$latin[SPECIES$id==sp], clade=SPECIES$clade[SPECIES$id==sp],
               genes=nrow(d$expr), samples=ncol(d$expr), timepoints=n_tp, replicates=n_rep,
               light_hours=LIGHT_HOURS[[sp]],
               paper_rhythmic=paper_rh, paper_pval_lt05=self_rh)
  }))
  save_table(ov, "00_species_overview")
  cat(sprintf("\n合计 %d 物种, 基因数 %d~%d\n", nrow(ov), min(ov$genes), max(ov$genes)))
  cat("== 00 完成 ==\n")
}

if (sys.nframe() == 0) main()
