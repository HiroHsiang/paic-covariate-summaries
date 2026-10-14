# 阶段 3 · 主运行（protocol v0.3 第 3.7 节）：96 个 DGM = D(4) × G(3) × EM 强度(2) × 结局(2) × 锚定(2)，
# 每个 DGM nsim = 1000，全部 7 种报告方式 × 8 种方法，不跑 bootstrap（覆盖率子集见 sims/11）。
# 可续跑：已有 output/main/<dgm>.rds 的 DGM 跳过。二分类先跑（快），生存后跑。
# 用法：Rscript sims/10_main_run.R [nsim] [ncores] [endpoint 过滤: binary|tte|all]

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(foreach); library(doParallel); library(doRNG) })
args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else 12
ep_filter <- if (length(args) >= 3) args[3] else "all"
pr <- readRDS("output/poplimit/dgm_params.rds")
out_dir <- "output/main"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress.log"); cat("start", format(Sys.time()), "nsim", nsim, "cores", ncores, "\n", file = log_file, append = TRUE)
t_all <- Sys.time()

# ---- DGM 网格（二分类在前） ----------------------------------------------------------
grid <- expand.grid(endpoint = c("binary", "tte"), anchored = c(TRUE, FALSE), em_mult = c(1, 2),
                    shape = c("linear", "quad", "threshold"), D = c("D1", "D2", "D3", "D4"), stringsAsFactors = FALSE)
grid <- grid[order(grid$endpoint != "binary", grid$D, grid$shape, grid$em_mult, !grid$anchored), ]
if (ep_filter != "all") grid <- grid[grid$endpoint == ep_filter, ]
grid$id <- sprintf("%s_%s_em%d_%s_%s", grid$D, grid$shape, grid$em_mult, grid$endpoint, ifelse(grid$anchored, "anch", "unanch"))
cat(nrow(grid), "DGMs in this run\n")

# ---- 超总体真值（缓存到 truths.rds） -------------------------------------------------
truth_file <- file.path(out_dir, "truths.rds")
truths <- if (file.exists(truth_file)) readRDS(truth_file) else NULL
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; key <- sprintf("%s_%s_em%d_%s", g$D, g$shape, g$em_mult, g$endpoint)   # 真值与锚定与否无关
  if (!is.null(truths) && key %in% truths$key) next
  tr <- dgm_truth(pr$specs[[g$D]], g$shape, pr$b_em[[g$shape]] * g$em_mult, g$endpoint, pr$cens_rate,
                  n = if (g$endpoint == "tte") 2e5 else 5e5)
  truths <- rbind(truths, data.frame(key = key, D = g$D, shape = g$shape, em_mult = g$em_mult, endpoint = g$endpoint,
                                     estimand = names(tr$truth), truth = tr$truth, mcse = tr$mcse, row.names = NULL))
  saveRDS(truths, truth_file)
}
cat("truths ready:", length(unique(truths$key)), "keys\n", file = log_file, append = TRUE)

# ---- 逐 DGM 运行（可续跑） ------------------------------------------------------------
for (k in seq_len(nrow(grid))) {
  g <- grid[k, ]; f_out <- file.path(out_dir, paste0(g$id, ".rds"))
  if (file.exists(f_out)) { cat("skip", g$id, "\n"); next }
  t0 <- Sys.time()
  dgm <- list(D = g$D, shape = g$shape, endpoint = g$endpoint, anchored = g$anchored, N_AC = 300, N_BC = 300, em_mult = g$em_mult)
  res <- run_dgm(dgm, pr, nsim = nsim, seed = 20261100 + k, ncores = ncores,
                 R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"),
                 methods = c("M1", "M1m", "M1c", "M7", "M2", "M3a", "M3b", "M5"), do_boot = FALSE, progress_file = log_file)
  res$dgm_id <- g$id
  saveRDS(res, f_out)
  cat(sprintf("done %s: %d rows, %.1f min, %.1f s/rep, errors %d\n", g$id, nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "mins")), mean(res$secs[!duplicated(res$rep)]),
              sum(grepl("^error", res$infeasible_type))), file = log_file, append = TRUE)
}
write_env_manifest(file.path(out_dir, "environment.md"))
writeLines(c("# Main run manifest", paste0("Run: ", format(t_all), " to ", format(Sys.time())),
             paste0("nsim = ", nsim, "; cores = ", ncores, "; DGMs = ", nrow(grid), " (", ep_filter, ")"),
             paste0("b_em (moderate): ", paste(names(pr$b_em), round(pr$b_em, 4), collapse = ", "), "; strong = 2x; cens_rate = ", signif(pr$cens_rate, 4)),
             "Seeds: 20261100 + DGM index (grid order: binary first; D, shape, em_mult, anchored)",
             paste0("Code: R/00-08 as of ", format(max(file.info(list.files("R", full.names = TRUE))$mtime)))),
           file.path(out_dir, "manifest.md"))
cat("all done", format(Sys.time()), "total", round(difftime(Sys.time(), t_all, units = "hours"), 2), "h\n", file = log_file, append = TRUE)
