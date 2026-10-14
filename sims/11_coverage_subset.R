# 阶段 3 · 覆盖率子集（protocol v1.0 决策 16，2026-10-09 修订）：
# G-computation（M3a/M3b）与 G-MAIC（M5）的 bootstrap 置信区间，全部 24 个锚定二分类 DGM
# （D1–D4 × 3 种形状 × EM 强度 2 档），R1 与 R5，300 次重复，B = 100。约 3 小时（12 核）。可续跑。
# MAIC 族的"重抽目标值"bootstrap 不在主研究中运行（试点已回答，见 output/pilot/coverage_B.csv）。
# 用法：Rscript sims/11_coverage_subset.R [nsim] [ncores]   （第 3 个参数忽略；B 固定为 100）

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 300
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
B <- 100                                                                    # 估计 SE 用 100 次足够（Efron），成本减半
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/coverage"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "B", B, "\n", file = log_file, append = TRUE)
grid <- expand.grid(em_mult = c(1, 2), shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid$id <- sprintf("%s_%s_em%d_binary_anch", grid$D, grid$shape, grid$em_mult)
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = "binary", anchored = TRUE, N_AC = 300, N_BC = 300, em_mult = g$em_mult)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261200 + k, ncores = ncores, R_set = c("R1", "R5"),
                 methods = c("M3a", "M3b", "M5", "M5s"), B = B, do_boot = TRUE, progress_file = log_file)   # M5s：soft G-MAIC
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.0f s/rep\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)])), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
