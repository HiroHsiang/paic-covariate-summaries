# 敏感性分析 B 的扩展（决策 26）：N_AC = N_BC = 1000 扩到全部 24 个锚定二分类强档 DGM（D1–D4 × 三种形状 × 锚定 / 无锚定），
# 方法含 M0、M1q、M7k1、M5s、M5sk1。sims/24 已跑的 8 个（D1/D3 × U 形/阈值）用相同的文件名跳过。种子 20261500 + 网格行号。
# 输出 output/sens_n1000/。可续跑。用法：Rscript sims/28_sens_n1000_full.R [nsim] [ncores]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/sens_n1000"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start(full)", format(Sys.time()), "nsim", nsim, "N 1000\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))

grid <- expand.grid(anchored = c(TRUE, FALSE), shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid$id <- sprintf("%s_%s_em2_binary_%s_N1000", grid$D, grid$shape, ifelse(grid$anchored, "anch", "unanch"))
grid$k <- 100 + seq_len(nrow(grid))                                         # 与 sims/24 的种子（20261500 + 1..8）不重叠

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }                  # sims/24 已跑的 8 个直接沿用
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = "binary", anchored = g$anchored, N_AC = 1000, N_BC = 1000, em_mult = 2)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261500 + g$k, ncores = ncores,
                 R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"),
                 methods = c("M0", "M1", "M1m", "M1c", "M1q", "M7", "M7k1", "M2", "M3a", "M3b", "M5", "M5s", "M5sk1"), do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id; res$N <- 1000
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.1f s/rep, errors %d\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)]), sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
