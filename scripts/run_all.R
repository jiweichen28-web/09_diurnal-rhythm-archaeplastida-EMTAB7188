# run_all.R —— 按编号依次跑 00→10 (聚合器)
# 注: 每步作为独立 Rscript 子进程运行(与单独跑一致), 避免 source 跳过 main 的陷阱。
# 运行: Rscript scripts/run_all.R   (02 JTK + 09 WGCNA 较慢, 全程约 10 分钟)

.DIR  <- dirname(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))]))
.RBIN <- file.path(R.home("bin"), "Rscript")
.STEPS <- c("00_setup_and_load.R","01_qc_pca.R","02_rhythm_detection.R",
            "03_phase_distribution.R","04_enrichment_clock.R","05_cross_species.R",
            "06_gene_age_rhythm.R","07_ortholog_phase_conservation.R","08_cellcycle_timing.R",
            "09_wgcna_module_diurnal.R","10_photoperiod_phase.R")

for (.s in .STEPS) {
  cat(sprintf("\n===== 运行 %s =====\n", .s))
  code <- system2(.RBIN, shQuote(file.path(.DIR, .s)), stdout="", stderr="")
  if (code != 0) stop(sprintf("步骤 %s 失败 (exit %d)", .s, code))
}
cat("\n===== run_all 全部完成 =====\n")
