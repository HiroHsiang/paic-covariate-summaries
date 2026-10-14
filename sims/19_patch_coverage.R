# 阶段 3 · 覆盖率子集的补跑（决策 22）：R5 下 M3b / M5 / M5s 的 bootstrap 区间用修订后的重建法重算，
# 网格与种子与 sims/11 完全相同（20261200 + k），结果写到 output/coverage_patch/，由汇总脚本替换 output/coverage 中对应的行
# （method ∈ {M3b, M5, M5s} 且 R == "R5"；R1 下不涉及重建，不变）。约 1 小时（12 核）。可续跑。
# 用法：Rscript sims/19_patch_coverage.R [nsim] [ncores]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 300
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
B <- 100
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/coverage_patch"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "B", B, "\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))
grid <- expand.grid(em_mult = c(1, 2), shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid$id <- sprintf("%s_%s_em%d_binary_anch", grid$D, grid$shape, grid$em_mult)
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = "binary", anchored = TRUE, N_AC = 300, N_BC = 300, em_mult = g$em_mult)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261200 + k, ncores = ncores, R_set = "R5",
                 methods = c("M3b", "M5", "M5s"), B = B, do_boot = TRUE, progress_file = log_file)
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.0f s/rep\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)])), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
