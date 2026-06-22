# 06_gene_age_rhythm.R —— 基因年龄(phylostratum) × 节律: 置换检验富集/缺失
# 输入: data/*.xlsx (Phylostratum + pVal 列)
# 输出: results/tables/06_phylo_enrichment.csv(各物种×各年龄层 观测/期望/经验p)
#       results/figures/06_phylo_rhythmic_frac.{pdf,png}(各年龄层节律比例, 按进化排序)
#       06_phylo_enrichment_heat.{pdf,png}(物种×年龄层 富集/缺失热图)
# 注: 复现论文方法——shuffle 基因-年龄层分配 1000 次得经验 p; 检验"老基因是否更可能有节律"。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)
NPERM <- 1000

# 单物种: 各年龄层节律基因观测数 + 置换零分布 -> 经验 p(富集/缺失)
# 注: pVal 列只对节律基因有值(非节律 = "NR" -> NA); 故 is_rh = pVal<0.05 且非 NA,
#     非节律基因要保留(它们是背景), 只按 Phylostratum 是否缺失过滤。
phylo_perm <- function(meta) {
  meta <- meta[!is.na(meta$Phylostratum) & meta$Phylostratum != "", ]
  is_rh <- !is.na(meta$pVal) & meta$pVal < 0.05
  strata <- intersect(PHYLO_ORDER, unique(meta$Phylostratum))
  obs <- sapply(strata, function(s) sum(is_rh & meta$Phylostratum==s))
  size <- sapply(strata, function(s) sum(meta$Phylostratum==s))
  # 置换: 打乱 is_rh 与年龄层的对应
  perm <- matrix(0, NPERM, length(strata), dimnames=list(NULL, strata))
  for (i in seq_len(NPERM)) {
    sh <- sample(is_rh)
    perm[i, ] <- sapply(strata, function(s) sum(sh & meta$Phylostratum==s))
  }
  p_enrich <- sapply(strata, function(s) (sum(perm[,s] >= obs[s]) + 1) / (NPERM + 1))
  p_deplete<- sapply(strata, function(s) (sum(perm[,s] <= obs[s]) + 1) / (NPERM + 1))
  exp_mean <- colMeans(perm)
  data.frame(stratum=strata, size=size, obs_rhythmic=obs, exp_rhythmic=round(exp_mean,1),
             obs_frac=round(obs/size*100,1), log2_obs_exp=round(log2((obs+1)/(exp_mean+1)),3),
             p_enrich=p_enrich, p_deplete=p_deplete, stringsAsFactors=FALSE)
}

main <- function() {
  cat("== 06 基因年龄 × 节律 (置换检验) ==\n")
  all_res <- list()
  for (sp in SPECIES$id) {
    d <- load_species(sp)
    r <- phylo_perm(d$meta); r$species <- sp
    all_res[[sp]] <- r
    sig <- r$stratum[r$p_enrich < 0.05]
    cat(sprintf("%-4s 富集年龄层(p<0.05): %s\n", sp, if(length(sig)) paste(sig,collapse=", ") else "无"))
  }
  res <- do.call(rbind, all_res); save_table(res, "06_phylo_enrichment")

  # 各物种: 节律比例 vs 年龄层(老->新)
  res$stratum <- factor(res$stratum, levels=PHYLO_ORDER)
  res$species <- factor(res$species, levels=SPECIES$id)
  labs_sp <- setNames(SPECIES$latin, SPECIES$id)
  p1 <- ggplot(res, aes(stratum, obs_frac, group=species, color=species)) +
    geom_line(linewidth=0.4) + geom_point(size=0.8) +
    scale_color_manual(values=CLADE_COL, name="Species", labels=labs_sp) +
    labs(x="Phylostratum (old → young)", y="Rhythmic genes (%)",
         title="Gene age vs rhythmicity: older strata are more rhythmic") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),
                           legend.text=element_text(face="italic", size=6))
  savefig(p1, "06_phylo_rhythmic_frac", w=7, h=4)

  # 富集/缺失热图: log2(obs/exp), 标显著
  res$signif <- ifelse(res$p_enrich<0.05, "+", ifelse(res$p_deplete<0.05, "-", ""))
  p2 <- ggplot(res, aes(stratum, species, fill=log2_obs_exp)) +
    geom_tile(color="white", linewidth=0.3) +
    geom_text(aes(label=signif), size=PT+1) +
    scale_fill_gradient2(low="#3B528B", mid="white", high="#B2182B", midpoint=0,
                         name="log2(obs/exp)") +
    scale_y_discrete(labels=labs_sp, limits=rev(SPECIES$id)) +
    labs(x="Phylostratum (old → young)", y=NULL,
         title="Rhythmic enrichment(+)/depletion(-) by gene age (perm p<0.05)") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),
                           axis.text.y=element_text(face="italic", size=6))
  savefig(p2, "06_phylo_enrichment_heat", w=7, h=4)

  # 汇总: 老基因(前4层) vs 新基因(Specific) 节律比例
  old_str <- PHYLO_ORDER[1:4]
  agg <- do.call(rbind, lapply(SPECIES$id, function(sp){
    s <- res[res$species==sp, ]
    of <- weighted.mean(s$obs_frac[s$stratum %in% old_str], s$size[s$stratum %in% old_str], na.rm=TRUE)
    yf <- s$obs_frac[s$stratum=="Specific"]
    data.frame(id=sp, old_frac=round(of,1), young_frac=ifelse(length(yf), yf, NA))
  }))
  cat("\n老基因(前4古老层) vs 新基因(Specific) 节律比例:\n")
  print(agg, row.names=FALSE)
  cat("== 06 完成 ==\n")
}

if (sys.nframe() == 0) main()
