# 阶段 3 · 补充运行：M0（未调整 Bucher / 朴素比较）与 M5s（soft G-MAIC：sMAIC 权重 + 回归标准化），
# 在主运行的全部 DGM 上，使用与 sims/10 相同的网格顺序和种子 → 每次重复的数据与主运行完全一致，可配对比较。
# 用法：Rscript sims/15_supp_m0_softgmaic.R [nsim] [ncores] [endpoint 过滤: binary|tte|all]
# 产出：output/supp/<dgm_id>.rds（列结构与主运行相同）。可续跑。

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
ep_filter <- if (length(args) >= 3) args[3] else "all"
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/supp"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))                                           # 主进程先建好偏正态分位数表，随 .sn_cache 导出给工作进程

# 与 sims/10 完全相同的网格与顺序（种子 = 20261100 + 行号 k）
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
                 R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"), methods = c("M0", "M5s"), do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.1f s/rep, errors %d\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)]), sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
