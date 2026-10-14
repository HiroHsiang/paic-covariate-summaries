# 补算（决策 39，2026-10-13）：无锚定生存比较中加权方法的 RMST 差的标准误改用 IJ（权重固定的稳健标准误），以便报告覆盖率。
# 原因：主运行 weighted_arm() 的 se(rmean) 把权重当频数权重，不随 ESS 下降而增大（R/06 rmst_ij() 的说明）。
# 做法：24 个无锚定生存 DGM，网格顺序与种子与 sims/10 相同（20261100 + k）→ 每次重复的数据与主运行完全相同；
#       对 7 种报告方式 × 加权方法（M1、M1m、M1c、M1q、M7、M7k1、M7k05）重拟权重，算 A 臂加权 KM 的 RMST(24) 及 IJ 标准误，
#       再与发表的 B 臂 RMST 的标准误合并。估计值应与主运行逐次相同（汇总时核对）。输出 output/patch_rmst_se/<id>.rds。可续跑。
# 用法：Rscript sims/41_patch_rmst_se.R [nsim] [ncores]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 4
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/patch_rmst_se"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "cores", ncores, "\n", file = log_file, append = TRUE)
invisible(sn_quantile_grid(SN_U))

grid <- expand.grid(endpoint = c("binary", "tte"), anchored = c(TRUE, FALSE), em_mult = c(1, 2),
                    shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid <- grid[order(grid$endpoint != "binary", grid$D, grid$shape, grid$em_mult, !grid$anchored), ]
grid$k <- seq_len(nrow(grid))
grid$id <- sprintf("%s_%s_em%d_%s_%s", grid$D, grid$shape, grid$em_mult, grid$endpoint, ifelse(grid$anchored, "anch", "unanch"))
grid <- grid[grid$endpoint == "tte" & !grid$anchored, ]

R_SET <- c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5")
METHODS <- c("M1", "M1m", "M1c", "M1q", "M7", "M7k1", "M7k05")

#' 一次重复：各报告方式 × 加权方法的 A 臂 RMST 与 IJ 标准误（规则同 analyse_replicate：报了均值时不跑 M1c；M1q 只在 R3/R3m/R4/R5）
rmst_patch_replicate <- function(rep) {
  pub <- rep$pub_bc$effect[["est"]]; se_pub <- rep$pub_bc$effect[["se"]]; d <- rep$ipd_ac
  rows <- list()
  for (R in R_SET) {
    vis <- visible_table(rep$pub_bc$table, R, 1); has_mean <- !is.null(vis$x1$mean)
    for (m in METHODS) {
      if (m == "M1c" && has_mean) next
      if (m == "M1q" && !R %in% c("R3", "R3m", "R4", "R5")) next
      fit <- tryCatch(fit_maic_weights(d[, COVS], vis, MAIC_RECIPES[[m]]), error = function(e) NULL)
      if (is.null(fit) || !fit$feasible) { rows[[length(rows) + 1]] <- data.frame(method = m, R = R, est_ab = NA_real_, se_ij = NA_real_, se_km = NA_real_); next }
      a <- d$trt == "A" & fit$w > 0
      ij <- rmst_ij(d$time[a], d$status[a], fit$w[a])
      km <- weighted_arm(d, fit$w, "tte")                                   # 主运行用的标准误（对照用）
      rows[[length(rows) + 1]] <- data.frame(method = m, R = R, est_ab = ij[["est"]] - pub,
                                             se_ij = sqrt(ij[["se"]]^2 + se_pub^2), se_km = sqrt(km[["se"]]^2 + se_pub^2))
    }
  }
  do.call(rbind, rows)
}

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  spec <- pr$specs[[g$D]]; b_em <- pr$b_em[[g$shape]] * g$em_mult; pars <- pars_for(spec)
  worker <- function(r) {
    rep <- gen_replicate(spec, g$shape, b_em, "tte", pr$cens_rate, FALSE, 300, 300, pars = pars)   # 与 run_dgm 相同：先生成数据（同一随机数流）
    x <- rmst_patch_replicate(rep); x$rep <- r
    if (r %% 100 == 0) cat(g$D, g$shape, g$em_mult, "rep", r, format(Sys.time(), "%H:%M:%S"), "\n", file = log_file, append = TRUE)
    x
  }
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores); doParallel::registerDoParallel(cl)
    objs <- setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed")
    objs <- objs[vapply(objs, function(x) is.function(get(x, globalenv())) || object.size(get(x, globalenv())) < 5e6, logical(1))]
    parallel::clusterExport(cl, c(objs, "spec", "b_em", "pars", "g"), envir = environment())
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(survival); library(sn); library(estmeansd) }) })
    res <- doRNG::`%dorng%`(foreach::foreach(r = seq_len(nsim), .combine = rbind, .options.RNG = 20261100 + g$k), worker(r))
    parallel::stopCluster(cl)
  } else {
    set.seed(20261100 + g$k); res <- do.call(rbind, lapply(seq_len(nsim), worker))
  }
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %.1f min\n", g$id, as.numeric(difftime(Sys.time(), t0, units = "mins"))), file = log_file, append = TRUE)
}
cat("all done", format(Sys.time()), "\n", file = log_file, append = TRUE)
