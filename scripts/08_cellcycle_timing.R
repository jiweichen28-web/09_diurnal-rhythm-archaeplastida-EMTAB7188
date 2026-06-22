# 08_cellcycle_timing.R —— 细胞周期基因的昼夜调控: 单细胞藻 vs 多细胞陆生植物
# 输入: data/*.xlsx (Annotation MapMan bin + pVal + phase)
# 输出: results/tables/08_cellcycle_rhythmic.csv(各物种 细胞周期基因 节律比例)
#       results/figures/08_cellcycle_rhythmic.{pdf,png}(单细胞 vs 多细胞 条形)
#       08_cellcycle_phase.{pdf,png}(节律细胞周期基因的峰值相位, 单细胞集中/多细胞分散)
# 注: 复现摘要结论——细胞分裂基因在单细胞藻受昼夜钟严格控制(夜间分裂), 多细胞陆生植物组织层面不明显。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)

UNICELLULAR <- c("Syn","Cpa","Ppu","Cre","Kni")
# MapMan/注释里细胞周期相关关键词
CC_PAT <- "cell\\.cycle|cell\\.division|cell cycle|cyclin|\\bCDK\\b|\\bCDC\\b|mitotic|DNA\\.synthesis"

main <- function() {
  cat("== 08 细胞周期基因昼夜调控 ==\n")
  rows <- list(); phase_rows <- list()
  for (sp in SPECIES$id) {
    d <- load_species(sp)
    ann <- d$meta$Annotation; ann[is.na(ann)] <- ""
    is_cc <- grepl(CC_PAT, ann, ignore.case=TRUE)
    is_rh <- !is.na(d$meta$pVal) & d$meta$pVal < 0.05
    n_cc <- sum(is_cc); n_rh <- sum(is_cc & is_rh)
    frac <- if (n_cc>0) n_rh/n_cc*100 else NA
    # 背景: 全基因组节律比例(对比细胞周期基因是否更/更不节律)
    bg_frac <- mean(is_rh) * 100
    rows[[sp]] <- data.frame(id=sp, latin=SPECIES$latin[SPECIES$id==sp],
      lifestyle=ifelse(sp %in% UNICELLULAR,"Unicellular","Multicellular"),
      n_cellcycle=n_cc, n_rhythmic=n_rh, cc_rhythmic_frac=round(frac,1),
      genome_rhythmic_frac=round(bg_frac,1))
    cat(sprintf("%-4s %-13s 细胞周期 %4d 基因, 节律 %.0f%% (全基因组 %.0f%%)\n",
                sp, ifelse(sp %in% UNICELLULAR,"单细胞","多细胞"), n_cc, frac, bg_frac))
    # 节律细胞周期基因的相位
    ph <- d$meta$phase[is_cc & is_rh & !is.na(d$meta$phase)]
    if (length(ph)) phase_rows[[sp]] <- data.frame(id=sp,
      lifestyle=ifelse(sp %in% UNICELLULAR,"Unicellular","Multicellular"), phase=ph)
  }
  rr <- do.call(rbind, rows); save_table(rr, "08_cellcycle_rhythmic")
  uni <- rr$cc_rhythmic_frac[rr$lifestyle=="Unicellular"]
  mul <- rr$cc_rhythmic_frac[rr$lifestyle=="Multicellular"]
  wt <- wilcox.test(uni, mul)
  cat(sprintf("\n细胞周期基因节律比例: 单细胞 %.0f%% vs 多细胞 %.0f%%  Wilcoxon p=%.3f\n",
              mean(uni,na.rm=TRUE), mean(mul,na.rm=TRUE), wt$p.value))

  # ── 条形图: 细胞周期 vs 全基因组 节律比例, 按生活方式分色 ──
  rr$id <- factor(rr$id, levels=SPECIES$id)
  pd <- rbind(data.frame(id=rr$id, frac=rr$cc_rhythmic_frac, set="Cell-cycle genes", lifestyle=rr$lifestyle),
              data.frame(id=rr$id, frac=rr$genome_rhythmic_frac, set="Whole genome", lifestyle=rr$lifestyle))
  p1 <- ggplot(pd, aes(id, frac, fill=set)) +
    geom_col(position="dodge", width=0.7) +
    scale_fill_manual(values=c("Cell-cycle genes"="#C0392B","Whole genome"="grey70"), name=NULL) +
    labs(x=NULL, y="Rhythmic genes (%)",
         title="Cell-cycle gene rhythmicity: tight in algae, weaker in land plants") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))
  savefig(p1, "08_cellcycle_rhythmic", w=7, h=3.8)

  # ── 相位分布: 节律细胞周期基因峰值时刻(单细胞应集中, 多细胞分散) ──
  pa <- do.call(rbind, phase_rows); pa$id <- factor(pa$id, levels=SPECIES$id)
  labs_sp <- setNames(SPECIES$latin, SPECIES$id)
  p2 <- ggplot(pa, aes(phase, fill=lifestyle)) +
    geom_histogram(binwidth=2, color="white", linewidth=0.1) +
    scale_fill_manual(values=c(Unicellular="#1ABC9C", Multicellular="#E67E22"), name=NULL) +
    facet_wrap(~id, scales="free_y", ncol=5, labeller=labeller(id=labs_sp)) +
    labs(x="Peak phase (ZT, h)", y="Cell-cycle genes",
         title="When cell-cycle genes peak: concentrated in algae, diffuse in land plants") +
    theme_pub_bw() + theme(strip.text=element_text(face="italic", size=5.5))
  savefig(p2, "08_cellcycle_phase", w=9, h=4)
  cat("== 08 完成 ==\n")
}

if (sys.nframe() == 0) main()
