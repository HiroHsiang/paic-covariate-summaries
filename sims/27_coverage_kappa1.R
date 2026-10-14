# 覆盖率子集的补充（决策 26）：κ = 1 的 soft G-MAIC（M5sk1）的 bootstrap 区间，与 sims/11 相同的 24 个锚定二分类 DGM、R1/R5、300 次、B = 100、
# 种子 20261200 + k。M7k1 的稳健 SE 覆盖率在 sims/26 里已经有。输出 output/coverage_kappa1/。可续跑。
# 用法：Rscript sims/27_coverage_kappa1.R [nsim] [ncores]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 300
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
B <- 100
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/coverage_kappa1"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "B", B, "\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))
grid <- expand.grid(em_mult = c(1, 2), shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid$id <- sprintf("%s_%s_em%d_binary_anch", grid$D, grid$shape, grid$em_mult)
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = "binary", anchored = TRUE, N_AC = 300, N_BC = 300, em_mult = g$em_mult)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261200 + k, ncores = ncores, R_set = c("R1", "R5"),
                 methods = "M5sk1", B = B, do_boot = TRUE, progress_file = log_file)
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.0f s/rep\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)])), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
