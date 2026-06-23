# 09_wgcna_module_diurnal.R —— 共表达模块的昼夜性: 构网→模块→模块特征基因 JTK→峰值相位
# 输入: data/*.xlsx
# 输出: results/tables/09_<id>_modules.csv(模块大小+特征基因节律 p+峰值相位)
#       09_module_summary.csv(各物种 模块数/节律模块数/相位覆盖)
#       results/figures/09_<id>_ME_phase.{pdf,png}(模块特征基因昼夜曲线, 按相位排序)
#       09_module_diurnal_frac.{pdf,png}(各物种节律模块比例)
# 注: 共表达模块=表达模式相似的基因群; 若模块特征基因(ME)自身有节律, 说明该模块整体随昼夜起伏。
#     聚焦样本足够(≥24)的物种; ME 用 JTK 检验节律, 模块相位=ME 峰值时刻。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)
suppressMessages({library(WGCNA); library(MetaCycle)}); options(stringsAsFactors=FALSE)

# 对单物种: 构共表达网络 -> 模块 -> ME -> JTK 检验 ME 节律
run_wgcna <- function(sp) {
  d <- load_species(sp)
  keep <- rowSums(is.na(d$expr))==0 & rowSums(d$expr>1, na.rm=TRUE) >= ncol(d$expr)/2
  X <- log2(d$expr[keep,,drop=FALSE] + 1)
  # 取变异最大的 3000 基因(WGCNA 标准做法, 控时间)
  v <- apply(X, 1, var); X <- X[order(-v)[1:min(3000,nrow(X))], ]
  datExpr <- t(X)                                    # 样本 x 基因
  # 软阈值: 自动选(目标 R^2>0.8); 选不到时按 WGCNA 经验值(样本数决定)回退
  sft <- pickSoftThreshold(datExpr, powerVector=c(1:10,12,14,16,18,20), verbose=0)
  pw <- sft$powerEstimate
  if (is.na(pw) || pw < 4) {
    n <- nrow(datExpr)                                # 经验回退(WGCNA FAQ, unsigned)
    pw <- if (n < 20) 9 else if (n < 30) 8 else 7
  }
  net <- blockwiseModules(datExpr, power=pw, TOMType="unsigned", minModuleSize=30,
           mergeCutHeight=0.25, numericLabels=TRUE, verbose=0, maxBlockSize=3500)
  MEs <- net$MEs                                      # 模块特征基因(样本 x 模块)
  zt <- d$design$zt
  # 每个模块 ME: 按 zt 排序跑 JTK(保留重复)
  ord <- order(zt); MEo <- MEs[ord,,drop=FALSE]; zto <- zt[ord]
  res <- lapply(colnames(MEs), function(m) {
    if (m=="ME0") return(NULL)                        # ME0=未分配灰色模块, 跳过
    tmp <- file.path(TAB_DIR, sprintf("_tmp_wg_%s.csv", sp))
    write.csv(data.frame(id=m, t(MEo[,m,drop=FALSE]), check.names=FALSE), tmp, row.names=FALSE)
    jk <- suppressWarnings(meta2d(infile=tmp, filestyle="csv", outdir=TAB_DIR, timepoints=zto,
            cycMethod="JTK", minper=24, maxper=24, outputFile=FALSE, releaseNote=FALSE)$JTK)
    unlink(tmp)
    data.frame(module=m, size=sum(net$colors==as.integer(sub("ME","",m))),
               ME_adjp=jk$ADJ.P, ME_phase=jk$LAG)
  })
  list(tab=do.call(rbind, res), sp=sp)
}

main <- function() {
  cat("== 09 WGCNA 模块昼夜性 ==\n")
  use_sp <- SPECIES$id[SPECIES$id %in% c("Cpa","Ppu","Cre","Kni","Ppa","Smo","Pab","Syn")]
  summ <- list()
  for (sp in use_sp) {
    r <- tryCatch(run_wgcna(sp), error=function(e){cat(sp,"出错:",conditionMessage(e),"\n"); NULL})
    if (is.null(r) || is.null(r$tab)) next
    tab <- r$tab; tab$rhythmic <- tab$ME_adjp < 0.05
    save_table(tab, sprintf("09_%s_modules", sp))
    n_mod <- nrow(tab); n_rh <- sum(tab$rhythmic, na.rm=TRUE)
    summ[[sp]] <- data.frame(id=sp, latin=SPECIES$latin[SPECIES$id==sp],
      lifestyle=ifelse(sp %in% c("Syn","Cpa","Ppu","Cre","Kni"),"Unicellular","Multicellular"),
      n_modules=n_mod, n_rhythmic=n_rh, frac_rhythmic=round(n_rh/n_mod*100,1))
    cat(sprintf("%-4s 模块 %d, 节律模块 %d (%.0f%%)\n", sp, n_mod, n_rh, n_rh/n_mod*100))
    # 节律模块按峰值相位排序, 画 模块×峰值相位 条形
    rh_tab <- tab[which(tab$rhythmic), ]; rh_tab <- rh_tab[order(rh_tab$ME_phase), ]
    if (nrow(rh_tab) >= 1) {
      pdat <- rh_tab; pdat$module <- factor(pdat$module, levels=pdat$module)
      p <- ggplot(pdat, aes(ME_phase, module, fill=ME_phase)) +
        geom_col(width=0.7) +
        scale_fill_gradientn(colors=c("#FDE725","#21908C","#440154","#21908C","#FDE725"),
                             limits=c(0,24), name="Phase") +
        labs(x="Module peak phase (ZT, h)", y="Rhythmic module",
             title=sprintf("%s: co-expression modules by peak time", SPECIES$latin[SPECIES$id==sp])) +
        theme_pub_bw() + theme(plot.title=element_text(face="italic"))
      savefig(p, sprintf("09_%s_ME_phase", sp), w=5, h=3.5)
    }
  }
  ss <- do.call(rbind, summ); save_table(ss, "09_module_summary")

  ss$id <- factor(ss$id, levels=SPECIES$id)
  p <- ggplot(ss, aes(id, frac_rhythmic, fill=lifestyle)) +
    geom_col(width=0.7) +
    scale_fill_manual(values=c(Unicellular="#1ABC9C", Multicellular="#E67E22"), name=NULL) +
    geom_text(aes(label=sprintf("%d/%d", n_rhythmic, n_modules)), vjust=-0.3, size=PT) +
    labs(x=NULL, y="Rhythmic modules (%)",
         title="Fraction of co-expression modules that are diurnally rhythmic") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))
  savefig(p, "09_module_diurnal_frac", w=6, h=3.8)
  cat("== 09 完成 ==\n")
}

if (sys.nframe() == 0) main()
