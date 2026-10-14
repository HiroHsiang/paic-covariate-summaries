# 分析调度（protocol v0.2 第 3.6–3.7 节）：对一次重复，在每种报告方式（可得性 R）下跑各方法，
# 并把所有结果拼成长表；再提供按 DGM 并行跑 nsim 次重复的驱动函数。
# 依赖：R/01–07（由调用脚本 source）。可打包代码：不用 setwd()，不写绝对路径。

#' 一次重复 × 一组报告方式 × 一组方法 → 长表（每行一个方法 × R）
#' @param methods 子集 of c("M1","M1m","M1c","M7","M2","M3a","M3b","M5")；M1c 只在无均值的报告方式下有意义
#' @param do_boot 是否跑 bootstrap 方差（慢；方差研究的子集才开）
analyse_replicate <- function(rep, endpoint, R_set = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"),
                              methods = c("M1", "M1m", "M1c", "M7", "M2", "M3a", "M3b", "M5"),
                              B = 200, do_boot = FALSE, kappa = 1.96, round_dp = 1, n_sim = 5000, anchored = TRUE) {
  rows <- lapply(R_set, function(R) {
    vis <- visible_table(rep$pub_bc$table, R, round_dp)
    has_mean <- !is.null(vis$x1$mean)
    out <- list()
    for (m in intersect(methods, c("M0", "M1", "M1m", "M1c", "M1q", "M7", "M7k1", "M7k05"))) {
      if (m == "M1c" && has_mean) next                                      # 报了均值时 M1c 与 M1 相同
      if (m == "M1q" && !R %in% c("R3", "R3m", "R4", "R5")) next            # 没有 SD 可放弃时（R1 = 只均值，R2/R2r 同 M1），M1q 与已有方法相同
      out[[m]] <- tryCatch(run_maic_method(rep, vis, m, endpoint, B, kappa, round_dp, do_boot, anchored),
                           error = function(e) err_row(m, conditionMessage(e)))
    }
    if ("M2" %in% methods) out$M2 <- tryCatch(run_stc_method(rep, vis, endpoint, anchored), error = function(e) err_row("M2", conditionMessage(e)))
    gm <- intersect(methods, c("M3a", "M3b", "M5", "M5s", "M5sk1", "M5sk05"))
    Qs <- NULL                                                              # 重建的边际分布：同一张表只拟合一次，M3b / M5 / M5s 共用
    if (length(setdiff(gm, "M3a")) && (!is.null(vis$x1$q) || !is.null(vis$x1$range)))
      Qs <- tryCatch({ q <- lapply(COVS, function(v) reconstruct_marginal(vis[[v]])); names(q) <- COVS; q }, error = function(e) NULL)
    for (m in gm) {
      if (m == "M3b" && is.null(vis$x1$q) && is.null(vis$x1$range)) next    # 没有分位数/范围时 M3b 与 M3a 相同
      out[[m]] <- tryCatch(run_gcomp_method(rep, vis, m, endpoint, if (do_boot) B else 0, n_sim, kappa, round_dp, anchored, Qs),
                           error = function(e) err_row(m, conditionMessage(e)))
    }
    d <- do.call(rbind, out); d$R <- R; d
  })
  d <- do.call(rbind, rows)
  # 样本真值：主估计目标（锚定 binary → logOR，tte → logHR_diff；无锚定 tte → RMST_diff）与次要估计目标
  ts <- rep$truth_sample
  prim <- if (endpoint == "tte" && !anchored) "RMST_diff" else names(ts)[1]
  sec  <- switch(endpoint, binary = "RD", tte = if (anchored) "RMST_diff" else NA, NA)
  d$truth_sample <- ts[[prim]]; d$truth_sample_sec <- if (!is.na(sec)) ts[[sec]] else NA_real_
  d$pub_est <- rep$pub_bc$effect[["est"]]; d$pub_se <- rep$pub_bc$effect[["se"]]
  d
}

#' 出错时的占位行（保证长表列齐，错误信息进 infeasible_type）
err_row <- function(method, msg)
  data.frame(method = method, est_ac = NA_real_, est_ab = NA_real_, est_sec = NA_real_, se_sec = NA_real_, ess = NA_real_,
             feasible = FALSE, fallback = NA_integer_, infeasible_type = paste0("error: ", substr(msg, 1, 60)),
             se_fixedw = NA_real_, se_sandwich = NA_real_, se_boot_indep = NA_real_, se_boot_joint = NA_real_, boot_fail = NA_real_)

#' 按 DGM 并行跑 nsim 次重复（doParallel + doRNG，可重复）
#' @param dgm list(D, shape, endpoint, anchored, N_AC, N_BC, em_mult)；em_mult = 1（中等）或 2（强）；@param params dgm_params.rds 的内容
#' @return 长表，附 DGM 标识、重复编号、每次重复的耗时
run_dgm <- function(dgm, params, nsim, seed, ncores = 1, R_set, methods, B = 200, do_boot = FALSE,
                    kappa = 1.96, round_dp = 1, n_sim = 5000, progress_file = NULL) {
  if (is.null(dgm$em_mult)) dgm$em_mult <- 1
  spec <- params$specs[[dgm$D]]; b_em <- params$b_em[[dgm$shape]] * dgm$em_mult; cens <- params$cens_rate
  pars <- pars_for(spec)                                                    # 该 DGM 的中心化常数（AC 形状随 D 变）
  worker <- function(i) {
    t0 <- Sys.time()
    rep <- gen_replicate(spec, dgm$shape, b_em, dgm$endpoint, cens, dgm$anchored, dgm$N_AC, dgm$N_BC, pars = pars)
    d <- analyse_replicate(rep, dgm$endpoint, R_set, methods, B, do_boot, kappa, round_dp, n_sim, dgm$anchored)
    d$rep <- i; d$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(progress_file) && i %% 25 == 0) cat(dgm$D, dgm$shape, "rep", i, format(Sys.time(), "%H:%M:%S"), "\n", file = progress_file, append = TRUE)
    d
  }
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores); doParallel::registerDoParallel(cl); on.exit(parallel::stopCluster(cl))
    # 导出全局环境里的全部函数与（不太大的）对象：R/ 各文件定义的常量（SCALE、COVS、PARS、SN_U、.sn_um、.sn_cache …）工作进程都要用
    objs <- setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed")
    objs <- objs[vapply(objs, function(x) is.function(get(x, globalenv())) || object.size(get(x, globalenv())) < 5e6, logical(1))]
    parallel::clusterExport(cl, objs)
    parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(survival); library(sn); library(estmeansd); library(numDeriv) }) })
    res <- doRNG::`%dorng%`(foreach::foreach(i = seq_len(nsim), .combine = rbind, .options.RNG = seed), worker(i))
  } else {
    set.seed(seed); res <- do.call(rbind, lapply(seq_len(nsim), worker))
  }
  res$D <- dgm$D; res$shape <- dgm$shape; res$endpoint <- dgm$endpoint; res$anchored <- dgm$anchored; res$em_mult <- dgm$em_mult
  res
}
