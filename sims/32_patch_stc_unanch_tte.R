# 阶段 3 · 补跑（2026-10-12，用户要求"补"）：无锚定生存的代入式 STC（M2）原先未实现（返回 NA），现按 R/07 的新实现
# （A 臂 Weibull AFT 在 BC 均值处的 RMST(24) 减去发表的 B 臂 RMST）在 24 个无锚定生存 DGM 上重算，全部 7 种报告方式。
# 网格顺序与种子与 sims/10 完全相同 → 每次重复的数据一致；结果写到 output/main_patch_stc/<dgm_id>.rds，
# 由 sims/14 在汇总时替换主运行中对应的行（method == "M2"，endpoint == "tte"，anchored == FALSE）。
# 用法：Rscript sims/32_patch_stc_unanch_tte.R [nsim] [ncores]。可续跑。

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/main_patch_stc"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))                                           # 与其他驱动一致：先建偏正态分位数表再导出

# 与 sims/10 完全相同的网格与顺序（种子 = 20261100 + 行号 k），只取无锚定生存
grid <- expand.grid(endpoint = c("binary", "tte"), anchored = c(TRUE, FALSE), em_mult = c(1, 2),
                    shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid <- grid[order(grid$endpoint != "binary", grid$D, grid$shape, grid$em_mult, !grid$anchored), ]
grid$k <- seq_len(nrow(grid))
grid$id <- sprintf("%s_%s_em%d_%s_%s", grid$D, grid$shape, grid$em_mult, grid$endpoint, ifelse(grid$anchored, "anch", "unanch"))
grid <- grid[grid$endpoint == "tte" & !grid$anchored, ]

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = g$endpoint, anchored = g$anchored, N_AC = 300, N_BC = 300, em_mult = g$em_mult)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261100 + g$k, ncores = ncores,
                 R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"), methods = "M2", do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.1f s/rep, errors %d\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)]), sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
