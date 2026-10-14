# 总体极限（N → ∞）下的偏倚计算（protocol v1.2 第 3.7 节，阶段 1）。
# 思路：用固定种子的大样本近似两个试验的总体（三维协变量，与主研究同一 DGM）；
# 在 AC 总体样本上按"基线表"（BC 边际分布的总体值，可选舍入）做指数倾斜，得到加权后的协变量分布；
# 用加权分布算出的边际效应与 BC 总体下的边际效应之差，就是样本无限大时仍然残留的偏倚。
# 蒙特卡洛误差只来自 n_pop（默认 5e4）的抽样（logOR 约 0.005），同一格内各策略共用同一样本，可配对比较。
# 依赖：R/01–06（由调用脚本 source）：基线表与约束的翻译与样本层面完全相同（visible_table、build_constraints_multi）。

#' 阶段 1 评价的提取策略（约束配方）。use = 用到的汇总量；soft = 其中软匹配的；convert = 没报均值时把中位数换算成均值/SD；
#' avail = 该策略需要的报告方式（决定基线表上"看得到"哪些数）
POPLIMIT_METHODS <- list(
  none          = list(use = character(0),                                   avail = "R1"),
  mean          = list(use = "mean",                                         avail = "R1"),
  mean_sd       = list(use = c("mean", "sd"),                                avail = "R1"),
  q_only        = list(use = "quantiles",                                    avail = "R2"),
  conv_iqr_mean = list(use = "mean", convert = "iqr",                        avail = "R2"),   # 中位数/IQR 换算后只用均值
  conv_iqr_msd  = list(use = c("mean", "sd"), convert = "iqr",               avail = "R2"),   # 换算后用均值 + SD
  conv_rng_msd  = list(use = c("mean", "sd"), convert = "range",             avail = "R2r"),  # 中位数/范围换算后用均值 + SD
  mean_prop     = list(use = c("mean", "prop"),                              avail = "R3"),
  mean_q        = list(use = c("mean", "quantiles"),                         avail = "R4"),
  mean_sd_q     = list(use = c("mean", "sd", "quantiles"),                   avail = "R4"),
  all_exact     = list(use = c("mean", "sd", "quantiles", "prop", "range"),  avail = "R5"),
  smaic         = list(use = c("mean", "sd", "quantiles", "prop", "range"),
                       soft = c("sd", "quantiles"),                          avail = "R5"),   # 均值、比例精确；SD、分位数软；范围截取
  oracle_g      = list(use = "oracle",                                       avail = "oracle")  # 直接平衡 E_BC[Σ g(X_j)]（span oracle）
)

#' BC 总体的"基线表"：三个协变量同一边际，所以三张子表相同。结构与 baseline_table() 一致（未舍入；舍入在 visible_table）
#' @param N_BC 用于 sMAIC 容差与范围的总体类比（样本 min/max ≈ 1/(2N) 与 1 − 1/(2N) 分位数）
bc_table_pop <- function(spec, N_BC = 300, cut = PARS$cut, mis_cut = 60) {
  d <- marg_dist(spec$bc$x); mom <- marg_moments(spec$bc$x, cut)
  ct <- list(mean = mom$mean, sd = mom$sd, q = d$q(c(0.25, 0.5, 0.75)), range = d$q(c(1 / (2 * N_BC), 1 - 1 / (2 * N_BC))),
             n = N_BC, p_cut = 1 - d$p(cut), cut = cut, p_mis = 1 - d$p(mis_cut), mis_cut = mis_cut)
  list(x1 = ct, x2 = ct, x3 = ct, n = N_BC)
}

#' 一个格子：给定 DGM 规格、EM 形状、b_EM，返回各提取策略在各尺度上的极限偏倚
#' @param scales 子集 of c("meandiff", "RD", "logOR", "logHR", "RMST")
#' @param prop_at "true"：n (%) 按真实切点报告；"mis"：按错位切点 mis_cut 报告（阈值形状的切点错位扫描）
#' @param n_pop 两个总体样本的大小；seeds 固定，使所有格子共用同一组样本
poplimit_cell <- function(spec, shape, b_em, cens_rate, scales = c("meandiff", "RD", "logOR", "logHR", "RMST"),
                          methods = names(POPLIMIT_METHODS), N_BC = 300, mis_cut = 60, prop_at = c("true", "mis"),
                          kappa = 1.96, round_dp = 1, n_pop = 5e4, seeds = c(ac = 20261008, bc = 20261009),
                          tt = seq(0.1, PARS$admin_cens, by = 0.1), tau = 24) {
  prop_at <- match.arg(prop_at); pars <- pars_for(spec)
  X_ac <- population_draw(spec$ac, n_pop, seeds[["ac"]]); X_bc <- population_draw(spec$bc, n_pop, seeds[["bc"]])
  tb <- bc_table_pop(spec, N_BC, pars$cut, mis_cut)
  need_tte <- any(c("logHR", "RMST") %in% scales)
  if (need_tte) G <- G_cens(tt, cens_rate, pars$admin_cens)

  # 各尺度的"效应泛函"：输入（总体样本，权重）→ 边际效应；生存两臂只算一次
  effects <- function(X, w) {
    out <- c()
    if ("meandiff" %in% scales) out["meandiff"] <- marginal_meandiff(X, w, "A", "C", shape, b_em, pars)
    if (any(c("RD", "logOR") %in% scales)) { mb <- marginal_binary(X, w, "A", "C", shape, b_em, pars); out["RD"] <- mb[["RD"]]; out["logOR"] <- mb[["logOR"]] }
    if (need_tte) {
      sA <- marginal_surv(X, w, "A", shape, b_em, tt, pars); sC <- marginal_surv(X, w, "C", shape, b_em, tt, pars)
      out["logHR"] <- cox_limit(sA$S, sA$h, sC$S, sC$h, G); out["RMST"] <- rmst(sA$S, tt, tau) - rmst(sC$S, tt, tau)
    }
    out[scales]
  }
  truth <- effects(X_bc, rep(1, n_pop))
  # span oracle 的目标：E_BC[Σ_j g(X_j)] = 3 × E_BC[g(X)]（由 BC 边际数值积分，精确）
  xg <- seq(-20, 140, by = 0.01); fb <- marg_dist(spec$bc$x)$d(xg); fb <- fb / sum(fb)
  target_g <- 3 * sum(g_em(xg, shape, pars$cut) * fb)

  out <- lapply(methods, function(m) {
    rec <- POPLIMIT_METHODS[[m]]
    if (length(rec$use) == 0) {                                       # 不调整
      w <- rep(1, n_pop); feas <- TRUE; k <- 0L
    } else if (identical(rec$use, "oracle")) {                        # span oracle：平衡 E_BC[Σ g]
      Hg <- g_em(X_ac$x1, shape, pars$cut) + g_em(X_ac$x2, shape, pars$cut) + g_em(X_ac$x3, shape, pars$cut)
      tw <- tilt_weights(cbind(Hg), target_g, 0); w <- tw$w; feas <- tw$feasible; k <- 1L
    } else {
      vis <- visible_table(tb, rec$avail, round_dp, prop_at = prop_at)
      cs <- build_constraints_multi(X_ac, vis, rec, kappa, round_dp)
      base <- as.numeric(cs$mask)                                      # 范围：截掉报告范围之外的 AC 人群
      if (is.null(cs$H)) { w <- base; feas <- TRUE; k <- 0L }
      else { tw <- tilt_weights(cs$H, cs$targets, cs$tol, base = base / sum(base)); w <- tw$w; feas <- tw$feasible; k <- ncol(cs$H) }
    }
    est <- effects(X_ac, w)
    data.frame(method = m, scale = scales, estimate = est, truth = truth, bias = est - truth,
               ess_frac = sum(w)^2 / sum(w^2) / n_pop, feasible = feas, n_constraints = k, row.names = NULL)
  })
  do.call(rbind, out)
}
