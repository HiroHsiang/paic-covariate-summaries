# 敏感性分析 B（决策 24，2026-10-11）：样本量 N_AC = N_BC = 1000（主分析 300）。偏倚不随 N 缩小、SE 随 N 缩小，
# 看 RMSE 的排名是否翻向"多纳入信息"的方法。8 个二分类 DGM：D1 / D3 × U 形 / 阈值 × 锚定 / 无锚定，强档；全部报告方式；
# 方法含 M0、M1q、M5s。真值与主运行相同（超总体，output/main/truths.rds）。种子 20261500 + k。输出 output/sens_n1000/。可续跑。
# 用法：Rscript sims/24_sens_n1000.R [nsim] [ncores]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/sens_n1000"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "N 1000\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))

grid <- expand.grid(anchored = c(TRUE, FALSE), shape = c("quad", "threshold"), D = c("D1", "D3"), stringsAsFactors = FALSE)
grid$id <- sprintf("%s_%s_em2_binary_%s_N1000", grid$D, grid$shape, ifelse(grid$anchored, "anch", "unanch"))

for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = "binary", anchored = g$anchored, N_AC = 1000, N_BC = 1000, em_mult = 2)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261500 + k, ncores = ncores,
                 R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"),
                 methods = c("M0", "M1", "M1m", "M1c", "M1q", "M7", "M2", "M3a", "M3b", "M5", "M5s"), do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id; res$N <- 1000
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.1f s/rep, errors %d\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)]), sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
