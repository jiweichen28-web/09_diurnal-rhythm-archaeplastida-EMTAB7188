# 00_install_packages.R —— 依赖安装 (仅首次)
# 运行: Rscript scripts/00_install_packages.R

cran <- c("ggplot2","dplyr","tidyr","readxl","ggrepel","RColorBrewer","pheatmap","reshape2")
# 节律检测: MetaCycle 含 JTK_CYCLE 实现(复现论文 JTK_Cycle 结果)
cran <- c(cran, "MetaCycle")
bioc <- c("clusterProfiler","org.At.tair.db")  # 拟南芥/水稻功能富集(模式物种有注释)

for (p in cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
if (length(bioc)) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
  for (p in bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
}
cat("依赖检查完成。\n")
