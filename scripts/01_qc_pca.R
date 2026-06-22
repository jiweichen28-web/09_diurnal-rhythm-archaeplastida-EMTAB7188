# 01_qc_pca.R —— 各物种 log2 + PCA(样本按昼夜时间上色, 看"时间环")
# 输入: data/*.xlsx
# 输出: results/figures/01_PCA_all_species.{pdf,png}(9宫格), 01_<id>_PCA 单图可选
#       results/tables/01_pca_variance.csv
# 注: 时序转录组 PCA 常呈环形(样本沿昼夜周期首尾相接), 是数据质量的直观信号。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)

main <- function() {
  cat("== 01 QC + PCA ==\n")
  pca_list <- list(); var_rows <- list()
  for (sp in SPECIES$id) {
    d <- load_species(sp)
    keep <- rowSums(is.na(d$expr)) == 0 & rowSums(d$expr > 1, na.rm=TRUE) >= 3
    X <- log2(d$expr[keep, , drop=FALSE] + 1)       # log2(TPM+1), 去 NA 与低表达
    pr <- prcomp(t(X), center=TRUE, scale.=FALSE)
    ve <- pr$sdev^2 / sum(pr$sdev^2)
    df <- data.frame(PC1=pr$x[,1], PC2=pr$x[,2], zt=d$design$zt, ld=d$design$ld,
                     species=sp, latin=SPECIES$latin[SPECIES$id==sp])
    pca_list[[sp]] <- df
    var_rows[[sp]] <- data.frame(id=sp, PC1_var=round(ve[1]*100,1), PC2_var=round(ve[2]*100,1),
                                 genes_used=nrow(X))
    cat(sprintf("%-4s PC1 %4.1f%% PC2 %4.1f%% (基因 %d)\n", sp, ve[1]*100, ve[2]*100, nrow(X)))
  }
  save_table(do.call(rbind, var_rows), "01_pca_variance")

  all_df <- do.call(rbind, pca_list)
  all_df$species <- factor(all_df$species, levels=SPECIES$id)
  # 面板标题用拉丁名缩写
  labs_sp <- setNames(SPECIES$latin, SPECIES$id)
  p <- ggplot(all_df, aes(PC1, PC2, color=zt)) +
    geom_point(aes(shape=ld), size=1.3) +
    scale_color_gradientn(colors=c("#FDE725","#5DC863","#21908C","#3B528B","#440154"),
                          name="ZT (h)") +
    scale_shape_manual(values=c(L=16, D=17), name="Light/Dark") +
    facet_wrap(~species, scales="free", ncol=4, labeller=labeller(species=labs_sp)) +
    labs(title="Diurnal time-series PCA across 9 Archaeplastida species") +
    theme_pub_bw() + theme(strip.text=element_text(face="italic", size=6))
  savefig(p, "01_PCA_all_species", w=9, h=6)
  cat("== 01 完成 ==\n")
}

if (sys.nframe() == 0) main()
