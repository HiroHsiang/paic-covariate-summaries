# 回归类估计量（protocol v0.2 第 5 节）：M2 代入均值的 STC、M3a/M3b G-computation、M5 G-MAIC，
# 以及"由基线表重建边际分布"（第 5.5 节）。方差一律用 bootstrap（重抽 AC，发表的 SE 另加）。
# 依赖：R/01–06（由调用脚本 source）。可打包代码：不用 setwd()，不写绝对路径。

# ---- 1. 由可见汇总量重建一个连续协变量的边际分布（返回分位数函数 Q(u)） -------------
#' 族：偏正态 SN(xi, omega, alpha)；报告了范围时在 [min, max] 内截断。
#' 拟合规则（2026-10-11 修订，决策 22：把所有可见的汇总量都用上）：
#'   只有均值/SD（R1、R3）：正态；只有三个分位数（R2）：由三个分位数闭式确定（精确）；
#'   只有中位数 + 范围（R2r）：把 min / max 当作 1/(2n)、1 − 1/(2n) 分位数，闭式确定；
#'   均值/SD 与分位数同时可见（R4、R5）：以各汇总量的 SE 为权重做加权最小二乘，三个参数拟合全部 4–5 个汇总量（含截断），
#'   初值取分位数闭式解。原先只用分位数、忽略已报告的 SD：混合分布下把 SD 低估 15%，截断分布下把均值拉低 0.6（主运行中期核查）。
.sn_cache <- new.env()
SN_ALPHA <- seq(-8, 8, by = 0.05)                                           # alpha 网格（超出时取端点）
SN_U     <- c(seq(0.0005, 0.02, length.out = 40), seq(0.025, 0.975, length.out = 381), seq(0.98, 0.9995, length.out = 40))  # 概率网格，两端加密

#' 标准偏正态在给定概率处的分位数随 alpha 变化的表，按概率向量缓存（每个会话算一次）
sn_quantile_grid <- function(probs, alpha_grid = SN_ALPHA) {
  key <- paste(signif(probs, 8), collapse = "_")
  if (!is.null(.sn_cache[[key]])) return(.sn_cache[[key]])
  Z <- t(vapply(alpha_grid, function(a) sn::qsn(probs, xi = 0, omega = 1, alpha = a, solver = "RFB"), numeric(length(probs))))
  .sn_cache[[key]] <- list(alpha = alpha_grid, Z = Z)
  .sn_cache[[key]]
}
#' 标准偏正态在 SN_U 网格上的分位数（alpha 在表内线性插值；表每个会话算一次，约 35 秒，存在 .sn_cache$Zfull）
sn_qstd_grid <- function(alpha) {
  Z <- .sn_cache$Zfull
  if (is.null(Z)) { Z <- sn_quantile_grid(SN_U)$Z; .sn_cache$Zfull <- Z }
  a <- min(max(alpha, SN_ALPHA[1]), SN_ALPHA[length(SN_ALPHA)])
  i <- findInterval(a, SN_ALPHA, all.inside = TRUE); w <- (a - SN_ALPHA[i]) / (SN_ALPHA[i + 1] - SN_ALPHA[i])
  (1 - w) * Z[i, ] + w * Z[i + 1, ]
}
#' 偏正态（可截断）的分位数函数与汇总量
sn_Q <- function(xi, om, al, lo = -Inf, hi = Inf) {
  zq <- sn_qstd_grid(al)
  plo <- if (is.finite(lo)) sn::psn((lo - xi) / om, 0, 1, al) else 0
  phi <- if (is.finite(hi)) sn::psn((hi - xi) / om, 0, 1, al) else 1
  if (phi - plo < 1e-6) { plo <- 0; phi <- 1 }                               # 范围与拟合矛盾时放弃截断
  function(u) xi + om * approx(SN_U, zq, xout = plo + u * (phi - plo), rule = 2)$y
}
.sn_um <- (seq_len(500) - 0.5) / 500                                        # 算均值 / SD 用的等概率中点网格（500 点，SD 误差 < 0.5%）
sn_summaries <- function(Q) { x <- Q(c(.sn_um, 0.25, 0.5, 0.75)); m <- mean(x[1:500]); c(mean = m, sd = sqrt(mean((x[1:500] - m)^2)), q1 = x[501], q2 = x[502], q3 = x[503]) }

#' 由三个分位数闭式确定 (xi, omega, alpha)：(q3 − q2)/(q2 − q1) 只依赖 alpha → 查表反解；再由尺度与位置解 omega、xi
sn_fit_quantiles <- function(probs, qs) {
  g <- sn_quantile_grid(probs)
  r_grid <- (g$Z[, 3] - g$Z[, 2]) / (g$Z[, 2] - g$Z[, 1]); r_obs <- (qs[3] - qs[2]) / (qs[2] - qs[1])
  r_obs <- min(max(r_obs, min(r_grid)), max(r_grid))                       # 超出表范围时取端点（alpha = ±8）
  al <- approx(r_grid, g$alpha, xout = r_obs, ties = "ordered")$y
  z  <- vapply(1:3, function(j) approx(g$alpha, g$Z[, j], xout = al)$y, numeric(1))
  om <- (qs[3] - qs[1]) / (z[3] - z[1]); xi <- qs[2] - om * z[2]
  c(xi = xi, om = om, al = al)
}

reconstruct_marginal <- function(vc) {
  n <- vc$n
  has_q3 <- !is.null(vc$q) && !any(is.na(vc$q)); has_med <- !is.null(vc$q) && !is.na(vc$q[2])
  has_ms <- !is.null(vc$mean) && !is.null(vc$sd); has_rng <- !is.null(vc$range)
  lo <- if (has_rng) vc$range[1] else -Inf; hi <- if (has_rng) vc$range[2] else Inf
  if (has_ms && !has_med) {                                                 # R1 / R3：正态（截断用闭式 CDF）
    if (!has_rng) return(function(u) qnorm(u, vc$mean, vc$sd))
    plo <- pnorm(lo, vc$mean, vc$sd); phi <- pnorm(hi, vc$mean, vc$sd)
    return(function(u) qnorm(plo + u * (phi - plo), vc$mean, vc$sd))
  }
  if (!has_q3 && !(has_med && has_rng)) stop("可见汇总量不足以重建边际分布")
  th <- if (has_q3) sn_fit_quantiles(c(0.25, 0.5, 0.75), vc$q)
        else sn_fit_quantiles(c(1 / (2 * n), 0.5, 1 - 1 / (2 * n)), c(lo, vc$q[2], hi))   # R2r：范围当作极端分位数
  if (has_ms && has_q3) {                                                   # R4 / R5：加权最小二乘，均值、SD、三个分位数一起拟合
    s_ref <- vc$sd
    tg <- c(mean = vc$mean, sd = vc$sd, q1 = vc$q[1], q2 = vc$q[2], q3 = vc$q[3])
    se <- c(s_ref / sqrt(n), s_ref * sqrt(2 / (n - 1)), 1.36 * s_ref / sqrt(n), 1.25 * s_ref / sqrt(n), 1.36 * s_ref / sqrt(n))
    obj <- function(p) { s <- sn_summaries(sn_Q(p[1], exp(p[2]), p[3], lo, hi)); sum(((s - tg) / se)^2) }
    op <- tryCatch(optim(c(th[["xi"]], log(th[["om"]]), th[["al"]]), obj, method = "Nelder-Mead", control = list(maxit = 200, reltol = 1e-6)),
                   error = function(e) NULL)
    if (!is.null(op) && all(is.finite(op$par)) && op$value <= obj(c(th[["xi"]], log(th[["om"]]), th[["al"]])))
      th <- c(xi = op$par[1], om = exp(op$par[2]), al = op$par[3])
  }
  sn_Q(th[["xi"]], th[["om"]], th[["al"]], lo, hi)
}

#' 模拟 BC 协变量（n_sim 行）：边际由可见汇总量重建（或正态），相关结构借用 AC IPD
#' @param marginals "normal"（M3a：均值/SD，必要时由中位数换算）或 "reconstruct"（M3b/M5：全部可见汇总量）
#' @param Qs 预先重建好的分位数函数 list(x1, x2, x3)（同一张表在 M3b / M5 / M5s 之间共用，省去重复拟合）
simulate_bc_covariates <- function(vis, ipd_ac, n_sim = 5000, marginals = c("normal", "reconstruct"), Qs = NULL) {
  marginals <- match.arg(marginals)
  Xac <- ipd_ac[, COVS]
  R <- cor(Xac); R[!is.finite(R)] <- 0; diag(R) <- 1                        # AC 的相关矩阵（借给 BC）
  Z <- matrix(rnorm(n_sim * 3), n_sim, 3) %*% chol(R); U <- pnorm(Z)
  qf <- function(vc) {
    if (marginals == "normal") {
      if (is.null(vc$mean)) {                                               # 只有中位数：换算
        cv <- if (!is.null(vc$q) && !any(is.na(vc$q))) estmeansd::qe.mean.sd(q1.val = vc$q[1], med.val = vc$q[2], q3.val = vc$q[3], n = vc$n)
              else estmeansd::qe.mean.sd(min.val = vc$range[1], med.val = vc$q[2], max.val = vc$range[2], n = vc$n)
        vc$mean <- cv$est.mean; vc$sd <- cv$est.sd }
      function(u) qnorm(u, vc$mean, vc$sd)
    } else reconstruct_marginal(vc)
  }
  out <- lapply(seq_along(COVS), function(j) { v <- COVS[j]
    Q <- if (marginals == "reconstruct" && !is.null(Qs[[v]])) Qs[[v]] else qf(vis[[v]]); Q(U[, j]) }); names(out) <- COVS
  as.data.frame(out)
}

# ---- 2. 结局模型（分析者的模型：协变量线性，治疗 × 效应修饰因子交互） ----------------
#' 锚定：y ~ trt × (x1 + x2 + x3)；无锚定（只有 A 臂）：y ~ x1 + x2 + x3
outcome_formula <- function(endpoint, anchored = TRUE) {
  rhs <- if (anchored) "trt * (x1 + x2 + x3)" else "x1 + x2 + x3"
  switch(endpoint, binary = as.formula(paste("y ~", rhs)), continuous = as.formula(paste("y ~", rhs)),
         tte = as.formula(paste("survival::Surv(time, status) ~", rhs)))
}
fit_outcome_model <- function(ipd, endpoint, w = NULL, anchored = TRUE) {
  if (anchored) ipd$trt <- factor(ipd$trt, levels = c("C", "A"))            # 系数 = A vs C
  ipd$.w <- if (is.null(w)) rep(1, nrow(ipd)) else w                        # 权重放进数据框，供公式环境取用
  ipd <- ipd[ipd$.w > 0, ]                                                  # 被范围截掉（权重 0）的个体不进模型（survreg 不接受零权重）
  fml <- outcome_formula(endpoint, anchored)
  switch(endpoint,
    binary     = glm(fml, family = if (all(ipd$.w == 1)) binomial else quasibinomial, data = ipd, weights = .w),   # 加权拟合只用于 G-MAIC（quasibinomial 只是免去非整数警告，点估计相同）；方差走 bootstrap
    continuous = lm(fml, data = ipd, weights = .w),
    tte        = survival::survreg(fml, data = ipd, weights = .w, dist = "weibull"))
}

# ---- 3. M2：代入均值的 STC（TSD 18 中心化；条件估计目标） --------------------------
#' BC 各协变量的"代入值"：报了均值用均值，否则用中位数近似
bc_centres <- function(vis) vapply(COVS, function(v) if (!is.null(vis[[v]]$mean)) vis[[v]]$mean else reconstruct_marginal(vis[[v]])(0.5), numeric(1))

#' 把协变量在 BC 均值处中心化后拟合，trt 主效应即 BC 均值处的条件效应
stc_plugin <- function(ipd_ac, vis, endpoint) {
  m <- bc_centres(vis)
  d <- ipd_ac; for (v in COVS) d[[v]] <- d[[v]] - m[[v]]
  f <- fit_outcome_model(d, endpoint)
  if (endpoint == "tte") c(est = -unname(coef(f)["trtA"]) / f$scale, se = NA_real_)        # Weibull AFT → PH 的 log HR
  else c(est = unname(coef(f)["trtA"]), se = unname(sqrt(vcov(f)["trtA", "trtA"])))
}

# ---- 4. G-computation：结局模型 + 在模拟的 BC 协变量上标准化 --------------------------
#' 用拟合好的模型在 Xsim 上预测 A、C 两臂的结局，取平均后作边际对比；同时给出 A 臂单独的边际结局（无锚定用）
#' 无锚定时模型不含 trt，只预测 A 臂：返回 armA（二分类 = 风险；生存 = RMST(τ)；连续 = 均值）
standardise <- function(fit, Xsim, endpoint, cens_rate_ac = NULL, tau = 24, pars = PARS, anchored = TRUE) {
  if (!anchored) {
    return(switch(endpoint,
      binary     = c(est = NA_real_, RD = NA_real_, armA = mean(predict(fit, Xsim, type = "response"))),
      continuous = c(est = NA_real_, armA = mean(predict(fit, Xsim))),
      tte = { lpA <- predict(fit, Xsim, type = "lp"); n <- nrow(Xsim)
        Tsim <- exp(lpA + fit$scale * log(rexp(n))); Cens <- pmin(rexp(n, cens_rate_ac), pars$admin_cens)
        d <- data.frame(time = pmin(Tsim, Cens), status = as.numeric(Tsim <= Cens))
        km <- survival::survfit(survival::Surv(time, status) ~ 1, data = d); tb <- summary(km, rmean = tau)$table
        c(est = NA_real_, RMST = NA_real_, armA = unname(tb["rmean"])) }))
  }
  XA <- cbind(trt = factor("A", levels = c("C", "A")), Xsim); XC <- cbind(trt = factor("C", levels = c("C", "A")), Xsim)
  switch(endpoint,
    binary = { pA <- mean(predict(fit, XA, type = "response")); pC <- mean(predict(fit, XC, type = "response"))
               c(est = qlogis(pA) - qlogis(pC), RD = pA - pC, armA = pA) },
    continuous = { mA <- mean(predict(fit, XA)); c(est = mA - mean(predict(fit, XC)), armA = mA) },
    tte = {                                                                 # 由 Weibull AFT 模拟合成队列的事件时间，再拟合 Cox
      lpA <- predict(fit, XA, type = "lp"); lpC <- predict(fit, XC, type = "lp"); n <- nrow(Xsim)
      Tsim <- exp(c(lpA, lpC) + fit$scale * log(rexp(2 * n)))
      Cens <- pmin(rexp(2 * n, cens_rate_ac), pars$admin_cens)
      d <- data.frame(time = pmin(Tsim, Cens), status = as.numeric(Tsim <= Cens), trt = rep(c(1, 0), each = n))
      cx <- survival::coxph(survival::Surv(time, status) ~ trt, data = d)
      km <- survival::survfit(survival::Surv(time, status) ~ trt, data = d); tb <- summary(km, rmean = tau)$table
      c(est = unname(coef(cx)), RMST = unname(tb[2, "rmean"] - tb[1, "rmean"]), armA = unname(tb[2, "rmean"])) })
}

#' AC 试验的删失率估计（指数假设）：供生存结局的 G-computation 模拟用
cens_rate_from <- function(ipd_ac) { d <- ipd_ac; sum(d$status == 0 & d$time < PARS$admin_cens) / sum(d$time) }

#' G-computation（M3a 正态边际 / M3b 重建边际）或 G-MAIC（M5：M1 权重加权拟合 + 重建边际）
#' @return 一行 data.frame：est_ac、est_ab、se（bootstrap + 发表 SE）
#' M5 = G-MAIC（Gao 2025：全部精确的 MAIC 权重 + 回归标准化）；M5s = soft G-MAIC（权重换成 sMAIC 的软约束权重，
#' 即 ITC 场景下的 approximate residual balancing 思路，Athey/Imbens/Wager 2018）
run_gcomp_method <- function(rep, vis, method = c("M3a", "M3b", "M5", "M5s", "M5sk1", "M5sk05"), endpoint, B = 200, n_sim = 5000, kappa = 1.96, round_dp = 1,
                             anchored = TRUE, Qs = NULL) {
  method <- match.arg(method)
  marg <- if (method == "M3a") "normal" else "reconstruct"
  cr <- if (endpoint == "tte") cens_rate_from(rep$ipd_ac) else NULL
  Xsim <- simulate_bc_covariates(vis, rep$ipd_ac, n_sim, marg, Qs)         # 每次重复固定一套合成协变量
  w_recipe <- switch(method, M5 = MAIC_RECIPES$M1, M5s = MAIC_RECIPES$M7, M5sk1 = MAIC_RECIPES$M7k1, M5sk05 = MAIC_RECIPES$M7k05, NULL)   # M5sk1 / M5sk05：soft G-MAIC，κ = 1 / 0.5 的权重（决策 26、35）
  w_of <- function(ipd) if (!is.null(w_recipe)) {
    f <- fit_maic_weights(ipd[, COVS], vis, w_recipe, kappa, round_dp)
    if (f$feasible) f$w else NULL } else rep(1, nrow(ipd))
  pub <- rep$pub_bc$effect[["est"]]; se_pub <- rep$pub_bc$effect[["se"]]; ps <- rep$pub_bc$effect_sec
  # 把标准化结果换算成 A vs B：锚定用两臂对比减发表的 B vs C；无锚定用 A 臂边际结局对发表的 B 臂
  to_ab <- function(st) {
    if (anchored) c(ab = st[["est"]] - pub, sec = switch(endpoint, binary = st[["RD"]] - ps[["est"]], tte = st[["RMST"]] - ps[["est"]], NA_real_))
    else switch(endpoint,
      binary = c(ab = qlogis(st[["armA"]]) - qlogis(pub), sec = st[["armA"]] - pub),
      c(ab = st[["armA"]] - pub, sec = NA_real_)) }
  one <- function(ipd) { w <- w_of(ipd); if (is.null(w)) return(c(ab = NA_real_, sec = NA_real_))
    to_ab(standardise(fit_outcome_model(ipd, endpoint, w, anchored), Xsim, endpoint, cr, anchored = anchored)) }
  est <- one(rep$ipd_ac)
  bs <- if (B > 0) { idx_by_arm <- split(seq_len(nrow(rep$ipd_ac)), rep$ipd_ac$trt)
    vapply(seq_len(B), function(b) { idx <- unlist(lapply(idx_by_arm, function(i) sample(i, replace = TRUE)))
      tryCatch(one(rep$ipd_ac[idx, ])[["ab"]], error = function(e) NA_real_) }, numeric(1)) } else NA_real_
  # 无锚定二分类的发表 SE 要换到 logit 尺度后再并入；锚定直接并入
  se_pub_ab <- if (anchored) se_pub else if (endpoint == "binary") se_pub / (pub * (1 - pub)) else se_pub
  data.frame(method = method, est_ac = NA_real_, est_ab = est[["ab"]], est_sec = est[["sec"]], se_sec = NA_real_, ess = NA_real_,
             feasible = is.finite(est[["ab"]]), fallback = NA_integer_, infeasible_type = "none",
             se_fixedw = NA_real_, se_sandwich = NA_real_, se_boot_indep = NA_real_,
             se_boot_joint = if (B > 0) sqrt(sd(bs, na.rm = TRUE)^2 + se_pub_ab^2) else NA_real_,
             boot_fail = if (B > 0) mean(is.na(bs)) else NA_real_)
}

#' M2 的一行（条件估计目标；SE 为模型法）。无锚定：在 BC 均值处预测 A 臂结局，对发表的 B 臂
#' 无锚定生存（2026-10-12 补）：A 臂 Weibull AFT 在 BC 均值处的生存曲线积分到 tau 得 RMST，减去发表的 B 臂 RMST(tau)
run_stc_method <- function(rep, vis, endpoint, anchored = TRUE, tau = 24) {
  pub <- rep$pub_bc$effect[["est"]]; se_pub <- rep$pub_bc$effect[["se"]]
  if (anchored) {
    s <- stc_plugin(rep$ipd_ac, vis, endpoint); est_ab <- s[["est"]] - pub; se <- sqrt(s[["se"]]^2 + se_pub^2)
  } else {
    newA <- as.data.frame(as.list(bc_centres(vis)))                        # 在 BC 的均值（或中位数）处预测
    f <- fit_outcome_model(rep$ipd_ac, endpoint, anchored = FALSE)
    if (endpoint == "tte") {
      lp <- unname(predict(f, newA, type = "lp"))                           # Weibull AFT：S(t) = exp(-(t / exp(lp))^(1 / scale))
      est_ab <- integrate(function(t) exp(-(t / exp(lp))^(1 / f$scale)), 0, tau)$value - pub
    } else if (endpoint == "binary") { pA <- predict(f, newA, type = "response"); est_ab <- qlogis(pA) - qlogis(pub) }
    else est_ab <- unname(predict(f, newA)) - pub
    se <- NA_real_
  }
  data.frame(method = "M2", est_ac = NA_real_, est_ab = est_ab, est_sec = NA_real_, se_sec = NA_real_, ess = NA_real_,
             feasible = is.finite(est_ab), fallback = NA_integer_, infeasible_type = "none",
             se_fixedw = se, se_sandwich = NA_real_, se_boot_indep = NA_real_, se_boot_joint = NA_real_, boot_fail = NA_real_)
}
