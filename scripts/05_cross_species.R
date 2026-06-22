# 05_cross_species.R —— 跨物种比较: 节律基因比例(进化梯度) + 相位分布对比
# 输入: results/tables/02_rhythmic_counts.csv + 02_<id>_rhythm.csv
# 输出: results/figures/05_rhythmic_gradient.{pdf,png}(单细胞 vs 多细胞)
#       05_phase_density.{pdf,png}(各物种相位密度叠加)
#       results/tables/05_cross_species_summary.csv
# 注: 论文核心结论——单细胞藻多数基因受昼夜调控, 多细胞陆生植物组织层面节律弱。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)

# 单细胞 vs 多细胞分组(进化/形态)
UNICELLULAR <- c("Syn","Cpa","Ppu","Cre","Kni")        # 蓝藻 + 单细胞藻
MULTICELLULAR <- c("Ppa","Smo","Pab","Osa","Ath")      # 多细胞陆生植物

main <- function() {
  cat("== 05 跨物种比较 ==\n")
  cc <- read.csv(file.path(TAB_DIR, "02_rhythmic_counts.csv"), stringsAsFactors=FALSE)
  cc$lifestyle <- ifelse(cc$id %in% UNICELLULAR, "Unicellular", "Multicellular")
  cc$id <- factor(cc$id, levels=SPECIES$id)

  # ── 节律比例: 单细胞 vs 多细胞 ──
  p1 <- ggplot(cc, aes(id, rhythmic_frac, fill=lifestyle)) +
    geom_col(width=0.7) +
    scale_fill_manual(values=c(Unicellular="#1ABC9C", Multicellular="#E67E22"), name=NULL) +
    geom_text(aes(label=sprintf("%.0f%%", rhythmic_frac)), vjust=-0.3, size=PT) +
    labs(x=NULL, y="Rhythmic genes (%)",
         title="Diurnal control: unicellular algae vs multicellular land plants") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
  savefig(p1, "05_rhythmic_gradient", w=6.5, h=3.8)

  # 两组均值对比 + Wilcoxon
  uni <- cc$rhythmic_frac[cc$lifestyle=="Unicellular"]
  mul <- cc$rhythmic_frac[cc$lifestyle=="Multicellular"]
  wt <- wilcox.test(uni, mul)
  cat(sprintf("单细胞节律比例 %.0f%% vs 多细胞 %.0f%%  Wilcoxon p=%.3f\n",
              mean(uni), mean(mul), wt$p.value))

  # ── 各物种相位密度叠加(看峰值时刻分布是否保守) ──
  phase_all <- list()
  for (sp in SPECIES$id) {
    r <- read.csv(file.path(TAB_DIR, sprintf("02_%s_rhythm.csv", sp)), stringsAsFactors=FALSE)
    rh <- r[which(r$jtk_adjp < 0.05), ]
    ph <- rh$paper_phase; ph[is.na(ph)] <- rh$jtk_lag[is.na(ph)]
    ph <- ph[!is.na(ph)]
    if (length(ph)) phase_all[[sp]] <- data.frame(species=sp,
      lifestyle=ifelse(sp %in% UNICELLULAR,"Unicellular","Multicellular"), phase=ph)
  }
  pa <- do.call(rbind, phase_all)
  pa$species <- factor(pa$species, levels=SPECIES$id)
  p2 <- ggplot(pa, aes(phase, color=species)) +
    geom_density(linewidth=0.5) + facet_wrap(~lifestyle, ncol=1) +
    scale_color_manual(values=CLADE_COL, name="Species") +
    labs(x="Peak phase (ZT, h)", y="Density",
         title="Peak-phase density: dawn/dusk enrichment across species") +
    theme_pub_bw()
  savefig(p2, "05_phase_density", w=6, h=4.5)

  # ── 汇总表 ──
  summ <- cc[, c("id","latin","clade","lifestyle","tested","rhythmic_self","rhythmic_frac",
                 "rhythmic_paper","jaccard_vs_paper")]
  summ <- summ[order(match(summ$id, SPECIES$id)), ]
  save_table(summ, "05_cross_species_summary")
  cat("== 05 完成 ==\n")
}

if (sys.nframe() == 0) main()
