# 07_ortholog_phase_conservation.R —— 跨物种 1:1 直系同源相位保守(论文标题结论)
# 输入: data/SupplementaryData14_orthogroups.xlsx + 各物种 phase
# 输出: results/tables/07_ortholog_phase_pairs.csv(物种对×共节律直系同源数×相位相关)
#       results/figures/07_phase_heatmap_examples.{pdf,png}(代表物种对 相位2D密度)
#       07_phase_conservation_matrix.{pdf,png}(全物种对 相位相关热图)
# 注: 1:1 直系同源(两物种各恰好1基因) + 两边都节律; 相位为环形(1-24h), 用环形相关。
#     论文核心: 即使跨10亿年进化, 直系同源仍倾向在一天同一时刻达峰(相位保守)。

source(file.path(dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))])), "_common.R")); set.seed(42)
suppressMessages(library(readxl))
NPERM <- 1000

# 环形相位差(周期24h): min(|a-b|, 24-|a-b|)
circ_diff <- function(a, b) { d <- abs(a - b); pmin(d, 24 - d) }

# 复现论文方法: 观测 Δphase 均值 vs 置换零分布(打乱配对) -> 经验 p
# 显著相似 = 观测 Δphase 比随机配对更小
dphase_perm <- function(pa, pb) {
  obs <- mean(circ_diff(pa, pb))
  perm <- replicate(NPERM, mean(circ_diff(pa, sample(pb))))
  p <- (sum(perm <= obs) + 1) / (NPERM + 1)
  list(obs=obs, exp=mean(perm), p=p)
}

# 各物种: GeneName -> phase(仅节律基因, pVal<0.05)
phase_lookup <- function(sp) {
  d <- load_species(sp)
  rh <- !is.na(d$meta$pVal) & d$meta$pVal < 0.05 & !is.na(d$meta$phase)
  setNames(d$meta$phase[rh], d$meta$GeneName[rh])
}

# orthogroup 表 -> 物种列内基因数
n_in <- function(cell) if (is.na(cell)||cell=="") 0L else length(strsplit(cell, ",\\s*")[[1]])
first_gene <- function(cell) trimws(strsplit(cell, ",\\s*")[[1]][1])

main <- function() {
  cat("== 07 直系同源相位保守 ==\n")
  og <- as.data.frame(read_excel(file.path("data","SupplementaryData14_orthogroups.xlsx"), skip=3))
  # 物种列名 -> id (按 SPECIES$latin 匹配 orthogroup 表头)
  sp_col <- setNames(SPECIES$id, SPECIES$latin)
  colnames(og)[match(names(sp_col), colnames(og))] <- sp_col[names(sp_col)[names(sp_col) %in% colnames(og)]]
  # 兜底: 逐列按拉丁名前缀匹配
  for (sp in SPECIES$id) {
    lat <- SPECIES$latin[SPECIES$id==sp]
    hit <- which(colnames(og)==lat | startsWith(colnames(og), substr(lat,1,15)))
    if (length(hit)) colnames(og)[hit[1]] <- sp
  }
  ph <- lapply(SPECIES$id, phase_lookup); names(ph) <- SPECIES$id

  pairs <- t(combn(SPECIES$id, 2))
  res <- list()
  for (k in seq_len(nrow(pairs))) {
    a <- pairs[k,1]; b <- pairs[k,2]
    if (!all(c(a,b) %in% colnames(og))) next
    na <- sapply(og[[a]], n_in); nb <- sapply(og[[b]], n_in)
    o11 <- which(na==1 & nb==1)                          # 1:1 直系同源
    ga <- sapply(og[[a]][o11], first_gene); gb <- sapply(og[[b]][o11], first_gene)
    pa <- ph[[a]][ga]; pb <- ph[[b]][gb]                  # 两边相位(节律才有)
    ok <- !is.na(pa) & !is.na(pb)
    n_both <- sum(ok)
    if (n_both >= 30) {
      dp <- dphase_perm(pa[ok], pb[ok])
      res[[k]] <- data.frame(sp1=a, sp2=b, n_1to1=length(o11), n_both_rhythmic=n_both,
                 dphase_obs=round(dp$obs,2), dphase_exp=round(dp$exp,2),
                 ratio=round(dp$obs/dp$exp,3), emp_p=dp$p)
    }
  }
  rr <- do.call(rbind, res); save_table(rr, "07_ortholog_phase_pairs")
  n_sig <- sum(rr$emp_p < 0.05)
  cat(sprintf("物种对(共节律≥30) %d; 显著相似(Δphase 小于随机, p<0.05) %d 对; 比值中位 %.2f\n",
              nrow(rr), n_sig, median(rr$ratio)))

  # ── 相位保守热图矩阵(全物种对): 用 Δphase 比值(obs/exp, <1=保守) ──
  m <- matrix(NA, length(SPECIES$id), length(SPECIES$id), dimnames=list(SPECIES$id, SPECIES$id))
  for (k in seq_len(nrow(rr))) { m[rr$sp1[k],rr$sp2[k]] <- rr$ratio[k]; m[rr$sp2[k],rr$sp1[k]] <- rr$ratio[k] }
  diag(m) <- NA
  mdf <- expand.grid(sp1=SPECIES$id, sp2=SPECIES$id)
  mdf$ratio <- mapply(function(i,j) m[i,j], mdf$sp1, mdf$sp2)
  mdf$sp1 <- factor(mdf$sp1, levels=SPECIES$id); mdf$sp2 <- factor(mdf$sp2, levels=rev(SPECIES$id))
  p1 <- ggplot(mdf, aes(sp1, sp2, fill=ratio)) +
    geom_tile(color="white", linewidth=0.3) +
    geom_text(aes(label=ifelse(is.na(ratio),"",sprintf("%.2f",ratio))), size=PT-0.3) +
    scale_fill_gradient2(low="#B2182B", mid="white", high="#3B528B", midpoint=1,
                         na.value="grey90", name="Δphase\nobs/exp") +
    labs(x=NULL, y=NULL,
         title="Ortholog phase conservation: Δphase obs/exp (<1 = conserved)") +
    theme_pub_bw() + theme(axis.text.x=element_text(angle=45,hjust=1))
  savefig(p1, "07_phase_conservation_matrix", w=6.5, h=5.5)

  # ── 代表物种对相位 2D 热图(论文 Fig3 对角线): 选 4 对 ──
  ex_pairs <- list(c("Cpa","Cre"), c("Cre","Kni"), c("Osa","Ath"), c("Cpa","Ath"))
  ex_df <- list()
  for (pr in ex_pairs) {
    a<-pr[1]; b<-pr[2]
    na <- sapply(og[[a]], n_in); nb <- sapply(og[[b]], n_in); o11 <- which(na==1&nb==1)
    ga <- sapply(og[[a]][o11], first_gene); gb <- sapply(og[[b]][o11], first_gene)
    pa <- ph[[a]][ga]; pb <- ph[[b]][gb]; ok <- !is.na(pa)&!is.na(pb)
    if (sum(ok)>=20) ex_df[[paste(a,b,sep="-")]] <-
      data.frame(pair=sprintf("%s vs %s", a, b), pa=pa[ok], pb=pb[ok])
  }
  edf <- do.call(rbind, ex_df)
  p2 <- ggplot(edf, aes(pa, pb)) +
    geom_bin2d(bins=12) + scale_fill_gradient(low="#deebf7", high="#08306b", name="Orthologs") +
    geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=LW, color="red") +
    facet_wrap(~pair, ncol=2) +
    labs(x="Peak phase species 1 (ZT, h)", y="Peak phase species 2 (ZT, h)",
         title="1:1 ortholog peak phases align on the diagonal (phase conserved)") +
    theme_pub_bw()
  savefig(p2, "07_phase_heatmap_examples", w=6.5, h=6)
  cat("== 07 完成 ==\n")
}

if (sys.nframe() == 0) main()
