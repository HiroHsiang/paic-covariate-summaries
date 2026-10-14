# 阶段 3 · 正式汇总（protocol v1.2 第 6 节；2026-10-12 重写）。
# 把所有运行并成一套性能表：主运行（output/main，M3b/M5 在 R2/R2r/R4/R5 下用 output/main_patch 替换）、补充运行（M0、M5s）、
# M1q（extra_m1q）、κ = 1 版本（full_kappa1：M7k1、M5sk1）、ML-NMR（mlnmr：M9 @R1 正态边际 / @R5 重建边际）、
# G 类方法的 bootstrap 区间（coverage + coverage_patch + coverage_kappa1，300 次）、三个敏感性分析（sens_shift、sens_n1000、sens_kappa）。
# 性能指标与 MCSE 按 Morris 等 2019 表 6 的公式手算（偏倚、经验 SE、RMSE、覆盖率、ESS、不可行率），并与 rsimsum 在一个 DGM 上核对；
# 另给相对"只均值"（M1m@R1）与"均值+SD"（M1@R1）的逐次重复配对差及其 MCSE。
# 产出：output/summary/perf_long.csv、perf_paired.csv、perf_gfamily_coverage.csv、perf_mlnmr.csv、perf_sens_*.csv、summary_tables.md、README.md
# 用法：Rscript sims/14_main_summary.R   （约 3–5 分钟）

truths <- readRDS("output/main/truths.rds")
dir.create("output/summary", showWarnings = FALSE)
t0 <- Sys.time()
prim <- function(ep, an) if (ep == "binary") "logOR" else if (an) "logHR_diff" else "RMST_diff"
secn <- function(ep, an) if (ep == "binary") "RD" else if (an) "RMST_diff" else NA
keep_cols <- c("method", "R", "rep", "est_ab", "est_sec", "ess", "feasible", "fallback", "infeasible_type", "se_fixedw", "se_sandwich", "se_boot_joint",
               "truth_sample", "truth_sample_sec", "D", "shape", "endpoint", "anchored", "em_mult", "dgm_id")
rd <- function(f) { d <- readRDS(f); d[, intersect(keep_cols, names(d))] }

# ---- 1. 逐 DGM 组装逐次重复的长表 -------------------------------------------------------------
ids <- sub("[.]rds$", "", list.files("output/main", pattern = "^D.*[.]rds$"))
assemble <- function(id) {
  d <- rd(file.path("output/main", paste0(id, ".rds")))
  pf <- file.path("output/main_patch", paste0(id, ".rds"))
  if (file.exists(pf)) { p <- rd(pf); sel <- d$method %in% c("M3b", "M5") & d$R %in% c("R2", "R2r", "R4", "R5"); d <- rbind(d[!sel, ], p[, names(d)]) }
  pf <- file.path("output/main_patch_stc", paste0(id, ".rds"))                 # 无锚定生存的 STC（sims/32，2026-10-12 补跑）：替换主运行的 NA 行
  if (file.exists(pf)) { p <- rd(pf); sel <- d$method == "M2"; d <- rbind(d[!sel, ], p[, names(d)]) }
  for (dir in c("output/supp", "output/extra_m1q", "output/full_kappa1", "output/full_kappa05")) { f <- file.path(dir, paste0(id, ".rds")); if (file.exists(f)) { x <- rd(f); d <- rbind(d, x[, names(d)]) } }
  pf <- file.path("output/patch_rmst_se", paste0(id, ".rds"))                  # 决策 39：无锚定 RMST 差改用 IJ 稳健标准误（sims/41）
  if (file.exists(pf)) { p <- readRDS(pf)
    m <- match(paste(d$method, d$R, d$rep), paste(p$method, p$R, p$rep)); hit <- !is.na(m)
    if (any(hit) && max(abs(d$est_ab[hit] - p$est_ab[m[hit]]), na.rm = TRUE) > 1e-8) stop("RMST 补算的估计值与主运行不一致：", id)
    d$se_fixedw[hit] <- p$se_ij[m[hit]] }
  d$dgm_id <- id; d
}
measures <- function(err, se = NULL, lo = NULL, hi = NULL, truth = NULL) {       # Morris 等 2019 表 6
  ok <- is.finite(err); e <- err[ok]; n <- length(e)
  if (n < 2) return(c(n = n, bias = NA, bias_mcse = NA, empSE = NA, empSE_mcse = NA, rmse = NA, rmse_mcse = NA, cover = NA, cover_mcse = NA, se_ratio = NA))
  b <- mean(e); s <- sd(e); mse <- mean(e^2); mse_mcse <- sqrt(sum((e^2 - mse)^2) / (n * (n - 1)))
  cov <- cov_mcse <- se_ratio <- NA
  if (!is.null(se)) { s2 <- se[ok]; ok2 <- is.finite(s2); if (sum(ok2) > 10) { cv <- abs(e[ok2]) <= qnorm(0.975) * s2[ok2]; cov <- mean(cv); cov_mcse <- sqrt(cov * (1 - cov) / sum(ok2)); se_ratio <- mean(s2[ok2]) / sd(e[ok2]) } }
  if (!is.null(lo)) { l <- lo[ok]; h <- hi[ok]; ok2 <- is.finite(l) & is.finite(h); cv <- l[ok2] <= truth & truth <= h[ok2]; cov <- mean(cv); cov_mcse <- sqrt(cov * (1 - cov) / sum(ok2)) }
  c(n = n, bias = b, bias_mcse = s / sqrt(n), empSE = s, empSE_mcse = s / sqrt(2 * (n - 1)), rmse = sqrt(mse), rmse_mcse = mse_mcse / (2 * sqrt(mse)),
    cover = cov, cover_mcse = cov_mcse, se_ratio = se_ratio)
}
se_for <- function(d) ifelse(d$method %in% c("M0", "M1", "M1m", "M1c", "M1q", "M7", "M7k1", "M7k05", "M2"), d$se_fixedw, NA_real_)   # G 类方法的区间见覆盖率子集

perf <- list(); paired <- list()
for (id in ids) {
  d <- assemble(id); ep <- d$endpoint[1]; an <- d$anchored[1]
  key <- sprintf("%s_%s_em%d_%s", d$D[1], d$shape[1], d$em_mult[1], ep)
  tr <- truths$truth[truths$key == key & truths$estimand == prim(ep, an)]; tr_mcse <- truths$mcse[truths$key == key & truths$estimand == prim(ep, an)]
  d$err <- d$est_ab - tr; d$se_use <- se_for(d); d$err_sample <- d$est_ab - d$truth_sample
  g <- split(d, list(d$R, d$method), drop = TRUE)
  perf[[id]] <- do.call(rbind, lapply(g, function(x) {
    m <- measures(x$err, x$se_use); ms <- measures(x$err_sample, x$se_use)
    data.frame(dgm_id = id, D = x$D[1], shape = x$shape[1], em_mult = x$em_mult[1], endpoint = ep, anchored = an, R = x$R[1], method = x$method[1],
               truth = tr, truth_mcse = tr_mcse, t(m), bias_sample_truth = ms[["bias"]], cover_sample_truth = ms[["cover"]],
               ess_mean = mean(x$ess[is.finite(x$err)], na.rm = TRUE), ess_p5 = unname(quantile(x$ess[is.finite(x$err)], 0.05, na.rm = TRUE)),
               infeasible = mean(!is.finite(x$est_ab)), fallback = mean(x$fallback[is.finite(x$err)] > 0, na.rm = TRUE),
               se_sandwich_ratio = if (all(is.na(x$se_sandwich))) NA else mean(x$se_sandwich, na.rm = TRUE) / sd(x$err[is.finite(x$err)]), row.names = NULL) }))
  # 配对差：同一次重复里 方法 − 参照（M1m@R1 与 M1@R1）
  ref1 <- d[d$method == "M1m" & d$R == "R1", c("rep", "err")]; names(ref1)[2] <- "err_ref1"
  ref2 <- d[d$method == "M1" & d$R == "R1", c("rep", "err")]; names(ref2)[2] <- "err_ref2"
  dd <- merge(merge(d, ref1, by = "rep"), ref2, by = "rep")
  paired[[id]] <- do.call(rbind, lapply(split(dd, list(dd$R, dd$method), drop = TRUE), function(x) {
    f <- function(v) { ok <- is.finite(v); c(mean(v[ok]), sd(v[ok]) / sqrt(sum(ok))) }
    a <- f(x$err - x$err_ref1); b <- f(x$err - x$err_ref2); a2 <- f(x$err^2 - x$err_ref1^2); b2 <- f(x$err^2 - x$err_ref2^2)
    data.frame(dgm_id = id, D = x$D[1], shape = x$shape[1], em_mult = x$em_mult[1], endpoint = ep, anchored = an, R = x$R[1], method = x$method[1],
               dbias_vs_means = a[1], dbias_vs_means_mcse = a[2], dmse_vs_means = a2[1], dmse_vs_means_mcse = a2[2],
               dbias_vs_meanSD = b[1], dbias_vs_meanSD_mcse = b[2], dmse_vs_meanSD = b2[1], dmse_vs_meanSD_mcse = b2[2], row.names = NULL) }))
}
perf <- do.call(rbind, perf); paired <- do.call(rbind, paired)

# ---- 2. ML-NMR（24 个锚定二分类 DGM；后验可信区间覆盖率） ------------------------------------------
if (dir.exists("output/mlnmr")) {   # 公开仓库不含 ML-NMR（sims/12）；目录不存在时跳过
ml <- do.call(rbind, lapply(list.files("output/mlnmr", pattern = "^D.*[.]rds$", full.names = TRUE), function(f) { d <- readRDS(f)
  key <- sprintf("%s_%s_em%d_binary", d$D[1], d$shape[1], d$em_mult[1]); tr <- truths$truth[truths$key == key & truths$estimand == "logOR"]
  d$err <- d$est_ab - tr
  do.call(rbind, lapply(split(d, d$R), function(x) { m <- measures(x$err, lo = x$cri_lo, hi = x$cri_hi, truth = tr)
    data.frame(dgm_id = d$dgm_id[1], D = x$D[1], shape = x$shape[1], em_mult = x$em_mult[1], endpoint = "binary", anchored = TRUE, R = x$R[1], method = "M9",
               marginals = if (x$R[1] == "R1") "normal" else "reconstructed", truth = tr, t(m), post_sd_ratio = mean(x$sd_post, na.rm = TRUE) / sd(x$err, na.rm = TRUE),
               rhat_bad = mean(x$infeasible_type == "rhat>1.05"), n_int = 1000, row.names = NULL) })) }))
write.csv(ml, "output/summary/perf_mlnmr.csv", row.names = FALSE)
}

# ---- 3. G 类方法的 bootstrap 区间（覆盖率子集：锚定二分类，R1/R5，300 次，B = 100） -------------------
cov_rows <- list()
for (f in list.files("output/coverage", pattern = "^D.*[.]rds$", full.names = TRUE)) {
  d <- readRDS(f); id <- sub("[.]rds$", "", basename(f))
  pf <- file.path("output/coverage_patch", basename(f)); if (file.exists(pf)) { p <- readRDS(pf); sel <- d$method %in% c("M3b", "M5", "M5s") & d$R == "R5"; d <- rbind(d[!sel, ], p[, names(d)]) }
  kf <- file.path("output/coverage_kappa1", basename(f)); if (file.exists(kf)) { k <- readRDS(kf); d <- rbind(d, k[, names(d)]) }
  kf5 <- file.path("output/coverage_kappa05", basename(f)); if (file.exists(kf5)) { k5 <- readRDS(kf5); d <- rbind(d, k5[, names(d)]) }   # 决策 35：κ = 0.5 的 soft G-MAIC
  key <- sprintf("%s_%s_em%d_binary", d$D[1], d$shape[1], d$em_mult[1]); tr <- truths$truth[truths$key == key & truths$estimand == "logOR"]; d$err <- d$est_ab - tr
  cov_rows[[id]] <- do.call(rbind, lapply(split(d, list(d$R, d$method), drop = TRUE), function(x) { m <- measures(x$err, x$se_boot_joint)
    data.frame(dgm_id = id, D = x$D[1], shape = x$shape[1], em_mult = x$em_mult[1], R = x$R[1], method = x$method[1], n = m[["n"]], bias = m[["bias"]], empSE = m[["empSE"]],
               cover_boot = m[["cover"]], cover_boot_mcse = m[["cover_mcse"]], se_ratio_boot = m[["se_ratio"]], boot_fail = mean(x$boot_fail, na.rm = TRUE), row.names = NULL) }))
}
gcov <- do.call(rbind, cov_rows); write.csv(gcov, "output/summary/perf_gfamily_coverage.csv", row.names = FALSE)
# 并入主表：G 类方法在锚定二分类 R1/R5 的覆盖率来自子集（n = 300）
perf <- merge(perf, gcov[, c("dgm_id", "R", "method", "cover_boot", "cover_boot_mcse", "se_ratio_boot")], by = c("dgm_id", "R", "method"), all.x = TRUE)
gfam <- perf$method %in% c("M3a", "M3b", "M5", "M5s", "M5sk1", "M5sk05")
perf$cover[gfam] <- perf$cover_boot[gfam]; perf$cover_mcse[gfam] <- perf$cover_boot_mcse[gfam]; perf$se_ratio[gfam] <- perf$se_ratio_boot[gfam]
perf$cover_source <- ifelse(perf$method %in% c("M0", "M1", "M1m", "M1c", "M1q", "M7", "M7k1", "M7k05"), "robust SE (fixed weights), 1000 reps",
                     ifelse(perf$method == "M2", "model-based SE, 1000 reps", ifelse(is.na(perf$cover), "", "IPD bootstrap B=100 (reported summaries and simulated population fixed), coverage subset 300 reps")))
perf <- perf[order(perf$endpoint != "binary", !perf$anchored, perf$em_mult, perf$D, perf$shape, perf$R, perf$method), ]
write.csv(perf, "output/summary/perf_long.csv", row.names = FALSE); write.csv(paired, "output/summary/perf_paired.csv", row.names = FALSE)

# ---- 4. 与 rsimsum 核对一个 DGM --------------------------------------------------------------
chk <- NULL
if (requireNamespace("rsimsum", quietly = TRUE)) {
  d <- assemble("D1_quad_em2_binary_anch"); tr <- truths$truth[truths$key == "D1_quad_em2_binary" & truths$estimand == "logOR"]
  x <- d[d$method == "M1" & d$R == "R5" & is.finite(d$est_ab), ]; x$se <- x$se_fixedw; x$truth <- tr
  s <- rsimsum::simsum(data = x, estvarname = "est_ab", true = "truth", se = "se"); ss <- summary(s)$summ
  mine <- perf[perf$dgm_id == "D1_quad_em2_binary_anch" & perf$method == "M1" & perf$R == "R5", ]
  chk <- sprintf("rsimsum check (D1_quad_em2_binary_anch, M1@R5): bias %.4f (rsimsum %.4f), empSE %.4f (%.4f), coverage %.3f (%.3f)",
                 mine$bias, ss$est[ss$stat == "bias"], mine$empSE, ss$est[ss$stat == "empse"], mine$cover, ss$est[ss$stat == "cover"])
  cat(chk, "\n")
}

# ---- 5. 敏感性分析 --------------------------------------------------------------------------
sens <- function(dir, truth_tab, pattern = "^D.*[.]rds$", extra = NULL) {
  do.call(rbind, lapply(list.files(dir, pattern = pattern, full.names = TRUE), function(f) { d <- readRDS(f); d <- d[, intersect(c(keep_cols, extra), names(d))]
    key <- sprintf("%s_%s_em%d_%s", d$D[1], d$shape[1], d$em_mult[1], d$endpoint[1]); tr <- truth_tab$truth[truth_tab$key == key & truth_tab$estimand == prim(d$endpoint[1], d$anchored[1])]
    d$err <- d$est_ab - tr; d$se_use <- se_for(d)
    do.call(rbind, lapply(split(d, list(d$R, d$method), drop = TRUE), function(x) { m <- measures(x$err, x$se_use)
      data.frame(dgm_id = d$dgm_id[1], D = x$D[1], shape = x$shape[1], em_mult = x$em_mult[1], endpoint = x$endpoint[1], anchored = x$anchored[1], R = x$R[1], method = x$method[1],
                 truth = tr, t(m), ess_mean = mean(x$ess[is.finite(x$err)], na.rm = TRUE), infeasible = mean(!is.finite(x$est_ab)), row.names = NULL) })) }))
}
sh <- sens("output/sens_shift", readRDS("output/sens_shift/truths.rds")); sh$delta_sd <- 0.45; write.csv(sh, "output/summary/perf_sens_shift045.csv", row.names = FALSE)
nn <- sens("output/sens_n1000", truths); nn$N <- 1000; write.csv(nn, "output/summary/perf_sens_n1000.csv", row.names = FALSE)
if (dir.exists("output/sens_kappa")) { kp <- sens("output/sens_kappa", truths); write.csv(kp, "output/summary/perf_sens_kappa.csv", row.names = FALSE) }

# ---- 6. 人读的表（Markdown）：记分卡 + 每个 结局 × 锚定 × 强度 的 R1/R5 表 ------------------------
meth <- c(M0 = "Unadjusted (Bucher)", M1m = "MAIC means only", M1c = "MAIC converted mean+SD", M1 = "MAIC all summaries exact (+fallback)", M1q = "MAIC mean + quartiles + n(%) (SD unused)",
          M7 = "sMAIC kappa 1.96", M7k1 = "sMAIC kappa 1", M7k05 = "sMAIC kappa 0.5", M2 = "Plug-in STC", M3a = "G-computation normal marginals", M3b = "G-computation reconstructed marginals",
          M5 = "G-MAIC", M5s = "sG-MAIC kappa 1.96", M5sk1 = "sG-MAIC kappa 1", M5sk05 = "sG-MAIC kappa 0.5", M9 = "ML-NMR")
D_lab <- c(D1 = "D1 both normal", D2 = "D2 both right-skew", D3 = "D3 both left-skew", D4 = "D4 BC truncated"); s_lab <- c(linear = "linear", quad = "U-shape", threshold = "threshold")
fmt <- function(x, d = 3) ifelse(is.na(x), "—", formatC(x, format = "f", digits = d))
md <- c("# Summary tables (protocol v1.2)", paste0("Generated ", Sys.time(), ". Truth = superpopulation marginal effect; 1000 replicates per DGM unless stated. MCSE in perf_long.csv."), if (!is.null(chk)) chk else "", "")
score <- function(ep, an, Rsel) {
  x <- perf[perf$endpoint == ep & perf$anchored == an & perf$R == Rsel & perf$method %in% names(meth), ]
  ag <- do.call(rbind, lapply(split(x, x$method), function(y) data.frame(method = y$method[1], mean_abs_bias = mean(abs(y$bias)), max_abs_bias = max(abs(y$bias)), mean_empSE = mean(y$empSE),
                                                                         mean_rmse = mean(y$rmse), min_cover = suppressWarnings(min(y$cover, na.rm = TRUE)), mean_ess = mean(y$ess_mean, na.rm = TRUE), n_cells = nrow(y))))
  ag <- ag[order(match(ag$method, names(meth))), ]
  c(sprintf("## Scorecard: %s, %s, availability %s (%d cells = D x shape x EM)", ep, if (an) "anchored" else "unanchored", Rsel, max(ag$n_cells)), "",
    "| Method | mean abs bias | max abs bias | mean emp. SE | mean RMSE | min coverage | mean ESS |", "|---|---|---|---|---|---|---|",
    sprintf("| %s | %s | %s | %s | %s | %s | %s |", meth[ag$method], fmt(ag$mean_abs_bias), fmt(ag$max_abs_bias), fmt(ag$mean_empSE), fmt(ag$mean_rmse), fmt(ifelse(is.finite(ag$min_cover), ag$min_cover, NA), 2), fmt(ag$mean_ess, 0)), "")
}
for (ep in c("binary", "tte")) for (an in c(TRUE, FALSE)) for (Rsel in c("R1", "R5")) md <- c(md, score(ep, an, Rsel))
cell <- function(ep, an, em) {
  x <- perf[perf$endpoint == ep & perf$anchored == an & perf$em_mult == em & perf$R %in% c("R1", "R5") & perf$method %in% names(meth), ]
  out <- c(sprintf("## %s, %s, %s EM: R1 (mean (SD)) vs R5 (all summaries)", ep, if (an) "anchored" else "unanchored", if (em == 1) "moderate" else "strong"), "",
           "| Scenario | Method | R1 bias | R1 emp. SE | R1 RMSE | R1 cov | R5 bias | R5 emp. SE | R5 RMSE | R5 cov |", "|---|---|---|---|---|---|---|---|---|---|")
  for (D in names(D_lab)) for (s in names(s_lab)) for (m in names(meth)) {
    a <- x[x$D == D & x$shape == s & x$method == m & x$R == "R1", ]; b <- x[x$D == D & x$shape == s & x$method == m & x$R == "R5", ]; if (!nrow(a) && !nrow(b)) next
    g <- function(z, v) if (nrow(z)) z[[v]][1] else NA
    out <- c(out, sprintf("| %s, %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |", D_lab[[D]], s_lab[[s]], meth[[m]], fmt(g(a, "bias")), fmt(g(a, "empSE")), fmt(g(a, "rmse")), fmt(g(a, "cover"), 2), fmt(g(b, "bias")), fmt(g(b, "empSE")), fmt(g(b, "rmse")), fmt(g(b, "cover"), 2)))
  }
  c(out, "")
}
for (ep in c("binary", "tte")) for (an in c(TRUE, FALSE)) for (em in c(2, 1)) md <- c(md, cell(ep, an, em))
writeLines(md, "output/summary/summary_tables.md")
writeLines(c("# output/summary", "", "- perf_long.csv: one row per DGM x availability x method; bias/empSE/RMSE/coverage with MCSE; ESS; infeasible and fallback rates; coverage source per method.",
             "- perf_paired.csv: paired differences (same replicate) of error and squared error versus MAIC means-only (M1m@R1) and versus MAIC mean+SD (M1@R1), with MCSE.",
             "- perf_gfamily_coverage.csv: bootstrap-interval coverage of G-computation / G-MAIC / soft G-MAIC (anchored binary, R1/R5, 300 reps, B = 100; R5 rows from the reconstruction patch; kappa-1 soft G-MAIC from coverage_kappa1).",
             "- perf_mlnmr.csv: ML-NMR (n_int = 1000) on the 24 anchored binary DGMs; coverage of 95% credible intervals.",
             "- perf_sens_shift045.csv / perf_sens_n1000.csv / perf_sens_kappa.csv: sensitivity analyses (BC shift 0.45 SD; N = 1000; sMAIC kappa 1 / 0.5).",
             "- summary_tables.md: scorecards and R1-vs-R5 tables by endpoint x anchoring x EM strength.",
             "Sources: output/main (+ main_patch for M3b/M5 at R2/R2r/R4/R5; + main_patch_stc for plug-in STC in unanchored survival), supp (M0, M5s), extra_m1q (M1q), full_kappa1 (M7k1, M5sk1), mlnmr, coverage(+patch, kappa1), sens_*."), "output/summary/README.md")
cat("done:", nrow(perf), "perf rows;", nrow(paired), "paired rows;", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
