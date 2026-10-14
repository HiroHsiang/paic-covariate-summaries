# 敏感性分析 · 重叠差（决策 23，2026-10-11；用户定 0.45 SD：0.6 时 ESS 只剩 20% 太极端）：BC 均值平移 0.45 SD（54.5）而非 0.3 SD（53），尺度比仍 1.2；
# 二分类，锚定 + 无锚定，D1–D4 × 三种形状，强档 EM（b_EM × 2，不重新校准：同样的效应修饰、两个试验离得更远），24 个 DGM × nsim。
# N = 300 时实际 ESS：只均值约 72%，均值+SD 35–49%，全精确 27–36%，sMAIC 40–56%（主网格：84–87% / 52–57% / 40–44% / 56–66%）。输出 output/sens_shift/；可续跑。
# 用法：Rscript sims/20_sens_shift.R [nsim] [ncores]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
delta_sd <- 0.45
pr <- readRDS("output/poplimit/dgm_params.rds")
pr$specs <- lapply(names(pr$specs), function(D) dgm_spec(D, delta_sd = delta_sd)); names(pr$specs) <- c("D1", "D2", "D3", "D4")
out_dir <- "output/sens_shift"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "delta_sd", delta_sd, "\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))
saveRDS(pr, file.path(out_dir, "dgm_params_shift.rds"))

grid <- expand.grid(anchored = c(TRUE, FALSE), shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid$em_mult <- 2; grid$endpoint <- "binary"
grid$id <- sprintf("%s_%s_em2_binary_%s_shift045", grid$D, grid$shape, ifelse(grid$anchored, "anch", "unanch"))

# 超总体真值（平移后的 BC 人群）
truth_file <- file.path(out_dir, "truths.rds")
truths <- if (file.exists(truth_file)) readRDS(truth_file) else NULL
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; key <- sprintf("%s_%s_em2_binary", g$D, g$shape)
  if (!is.null(truths) && key %in% truths$key) next
  tr <- dgm_truth(pr$specs[[g$D]], g$shape, pr$b_em[[g$shape]] * 2, "binary", NULL, n = 5e5)
  truths <- rbind(truths, data.frame(key = key, D = g$D, shape = g$shape, em_mult = 2, endpoint = "binary",
                                     estimand = names(tr$truth), truth = tr$truth, mcse = tr$mcse, row.names = NULL))
  saveRDS(truths, truth_file)
}

for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = "binary", anchored = g$anchored, N_AC = 300, N_BC = 300, em_mult = 2)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261400 + k, ncores = ncores,
                 R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"),
                 methods = c("M0", "M1", "M1m", "M1c", "M1q", "M7", "M2", "M3a", "M3b", "M5", "M5s"), do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id; res$delta_sd <- delta_sd
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min, %.1f s/rep, errors %d\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              mean(res$secs[!duplicated(res$rep)]), sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
