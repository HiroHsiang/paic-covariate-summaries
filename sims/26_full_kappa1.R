# 补充运行（决策 26，2026-10-12）：κ = 1 的 sMAIC（M7k1）与 κ = 1 权重的 soft G-MAIC（M5sk1）跑全情景。
# 依据：κ 敏感性（sims/23）显示 κ 1.96 → 1 把 U 形 EM 下的偏倚减半、ESS 只降 9%、右偏下几乎不把 SD 过头带回来；N = 1000 时 κ 1.96 的残留偏倚开始拖累 RMSE。
# 全部 96 个 DGM、7 种报告方式，网格顺序与种子与 sims/10 相同 → 可与主运行逐次重复配对。输出 output/full_kappa1/。可续跑。
# 用法：Rscript sims/26_full_kappa1.R [nsim] [ncores] [endpoint 过滤: binary|tte|all]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
ep_filter <- if (length(args) >= 3) args[3] else "all"
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/full_kappa1"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))

grid <- expand.grid(endpoint = c("binary", "tte"), anchored = c(TRUE, FALSE), em_mult = c(1, 2),
                    shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid <- grid[order(grid$endpoint != "binary", grid$D, grid$shape, grid$em_mult, !grid$anchored), ]
grid$k <- seq_len(nrow(grid))
grid$id <- sprintf("%s_%s_em%d_%s_%s", grid$D, grid$shape, grid$em_mult, grid$endpoint, ifelse(grid$anchored, "anch", "unanch"))
if (ep_filter != "all") grid <- grid[grid$endpoint == ep_filter, ]

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = g$endpoint, anchored = g$anchored, N_AC = 300, N_BC = 300, em_mult = g$em_mult)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261100 + g$k, ncores = ncores,
                 R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"), methods = c("M7k1", "M5sk1"), do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.1f s/rep, errors %d\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)]), sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
