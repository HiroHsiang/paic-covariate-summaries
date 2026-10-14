# 阶段 3 · 补跑（决策 22，2026-10-11）：重建边际分布的方法（M3b、M5）在用到分位数 / 范围的报告方式（R2、R2r、R4、R5）下，
# 用修订后的重建法（把均值、SD 与分位数一起拟合，见 R/07）重新计算；网格顺序与种子与 sims/10 完全相同 → 每次重复的数据一致，
# 结果写到 output/main_patch/<dgm_id>.rds，由 sims/14 在汇总时替换主运行中对应的行（method ∈ {M3b, M5} 且 R ∈ {R2, R2r, R4, R5}）。
# 用法：Rscript sims/17_patch_reconstruct.R [nsim] [ncores] [endpoint 过滤: binary|tte|all]。可续跑。

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
ep_filter <- if (length(args) >= 3) args[3] else "all"
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/main_patch"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
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
                 R_set = c("R2", "R2r", "R4", "R5"), methods = c("M3b", "M5"), do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.1f s/rep, errors %d\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)]), sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
