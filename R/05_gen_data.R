# 数据生成（protocol v1.2 第 3.5–3.6 节）：一次重复 = AC 试验 IPD + BC 试验 IPD（内部保留）+ BC 的"发表物"。
# 发表物 = B vs C 的效应估计（带 SE）+ 基线表（三个连续协变量各自的均值、SD、四分位数、范围、切点以上比例）。
# 可得性（R1–R5）和舍入在分析时由 visible_table() 施加，这样同一份数据能配对比较不同的报告方式。
# 依赖：R/01_dgm_spec.R、R/02_outcome_model.R。可打包代码：不用 setwd()，不写绝对路径。

COVS <- c("x1", "x2", "x3")                                            # 协变量列名（三个都是连续变量）

#' Bender 法生成 Weibull–Cox 事件时间：T = scale0 · (−log U / exp(η))^(1/shape)
r_weibull_cox <- function(eta, pars = PARS) {
  U <- runif(length(eta))
  weib_scale0(pars) * (-log(U) / exp(eta))^(1 / pars$weib_shape)
}

#' 生成一个试验的 IPD（等比例随机分配）
#' @param spec_trial dgm_spec()$ac 或 $bc；@param arms 如 c("A", "C")；单臂试验传 "A"
#' @param pars 须是 pars_for(spec) 的结果（含该 DGM 的中心化常数）
#' @return data.frame(trt, x1–x3, 结局列：binary → y；tte → time, status；continuous → y)
gen_trial <- function(spec_trial, n, arms, shape, b_em, endpoint, cens_rate = NULL, pars = PARS) {
  X   <- population_draw(spec_trial, n)
  trt <- sample(rep(arms, length.out = n))                       # 等比例随机分配
  eta <- lin_pred(X, trt, shape, b_em, endpoint, pars)
  out <- switch(endpoint,
    binary     = data.frame(y = rbinom(n, 1, plogis(eta))),
    continuous = data.frame(y = eta + rnorm(n)),
    tte        = {                                                # 随机删失 + 行政删失，两试验同参数
      Tlat <- r_weibull_cox(eta, pars); Cens <- pmin(rexp(n, cens_rate), pars$admin_cens)
      data.frame(time = pmin(Tlat, Cens), status = as.numeric(Tlat <= Cens))
    })
  cbind(data.frame(trt = factor(trt, levels = arms)), X, out)
}

#' 基线表（未舍入）：每个协变量给均值、SD、Q1/中位数/Q3、最小/最大、切点以上比例（可另给错位切点的比例）
baseline_table <- function(X, cut = PARS$cut, mis_cut = NULL) {
  cont <- function(v) {
    ct <- list(mean = mean(v), sd = sd(v), q = unname(quantile(v, c(0.25, 0.5, 0.75))), range = range(v), n = length(v),
               p_cut = mean(v > cut), cut = cut)
    if (!is.null(mis_cut)) { ct$p_mis <- mean(v > mis_cut); ct$mis_cut <- mis_cut }
    ct
  }
  tb <- lapply(COVS, function(v) cont(X[[v]])); names(tb) <- COVS
  tb$n <- nrow(X)
  tb
}

#' 把 BC 试验的 IPD 压缩成发表物：效应估计（B vs C，边际）+ 基线表
#' @param unanchored TRUE 时只有 B 臂；效应项改为 B 臂的结局汇总（事件比例 / RMST）
publish_bc <- function(ipd_bc, endpoint, cut = PARS$cut, mis_cut = NULL, unanchored = FALSE, tau = 24) {
  X <- ipd_bc[, COVS]
  eff <- if (!unanchored) {
    switch(endpoint,
      binary     = { f <- glm(y ~ trt, family = binomial, data = ipd_bc); c(est = unname(coef(f)[2]), se = unname(sqrt(vcov(f)[2, 2]))) },
      continuous = { f <- lm(y ~ trt, data = ipd_bc);                    c(est = unname(coef(f)[2]), se = unname(sqrt(vcov(f)[2, 2]))) },
      tte        = { f <- survival::coxph(survival::Surv(time, status) ~ trt, data = ipd_bc)
                     c(est = unname(coef(f)[1]), se = unname(sqrt(vcov(f)[1, 1]))) })
  } else {
    switch(endpoint,
      binary     = { p <- mean(ipd_bc$y); c(est = p, se = sqrt(p * (1 - p) / nrow(ipd_bc))) },
      continuous = c(est = mean(ipd_bc$y), se = sd(ipd_bc$y) / sqrt(nrow(ipd_bc))),
      tte        = { km <- survival::survfit(survival::Surv(time, status) ~ 1, data = ipd_bc)   # KM 的 RMST(τ)
                     s  <- summary(km, rmean = tau)$table; c(est = unname(s["rmean"]), se = unname(s["se(rmean)"])) })
  }
  # 次要估计目标的发表值：锚定二分类给 B vs C 的风险差；锚定生存给 B − C 的 RMST(τ)；其余 NA
  eff_sec <- c(est = NA_real_, se = NA_real_)
  if (!unanchored && endpoint == "binary") {
    pB <- mean(ipd_bc$y[ipd_bc$trt == "B"]); pC <- mean(ipd_bc$y[ipd_bc$trt == "C"]); nB <- sum(ipd_bc$trt == "B"); nC <- sum(ipd_bc$trt == "C")
    eff_sec <- c(est = pB - pC, se = sqrt(pB * (1 - pB) / nB + pC * (1 - pC) / nC))
  }
  if (!unanchored && endpoint == "tte") {
    km <- survival::survfit(survival::Surv(time, status) ~ trt, data = ipd_bc); tb <- summary(km, rmean = tau)$table   # 行序 = C, B
    eff_sec <- c(est = unname(tb[2, "rmean"] - tb[1, "rmean"]), se = unname(sqrt(tb[2, "se(rmean)"]^2 + tb[1, "se(rmean)"]^2)))
  }
  list(effect = eff, effect_sec = eff_sec, table = baseline_table(X, cut, mis_cut), endpoint = endpoint, unanchored = unanchored)
}

#' 可得性：按报告情景把完整基线表变成分析者"看得到"的汇总量（含舍入）；三个协变量按同一方式报告
#' R1 均值(SD)；R2 中位数[IQR]；R2r 中位数(范围)；R3 均值(SD)+n(%)>切点；R3m 同 R3 但切点错位；
#' R4 均值(SD)+中位数[IQR]；R5 全部。未报告的项为 NULL。
#' @param prop_at "true" 按真实切点报告 n (%)；"mis" 按错位切点报告（R3m 总是错位；阶段 1 的切点扫描用）
visible_table <- function(tb, R = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"), round_dp = 1, prop_at = c("true", "mis")) {
  R <- match.arg(R); prop_at <- match.arg(prop_at); rd <- function(v) round(v, round_dp)
  use_mis <- R == "R3m" || prop_at == "mis"
  vis_cont <- function(ct) {
    out <- list(n = ct$n)
    if (R %in% c("R1", "R3", "R3m", "R4", "R5")) { out$mean <- rd(ct$mean); out$sd <- rd(ct$sd) }
    if (R %in% c("R2", "R2r"))                   out$q <- c(NA, rd(ct$q[2]), NA)      # 只有中位数
    if (R %in% c("R2", "R4", "R5"))              out$q <- rd(ct$q)                     # Q1 / 中位数 / Q3
    if (R %in% c("R2r", "R5"))                   out$range <- rd(ct$range)             # 范围（同样舍入到 1 位小数）
    if (R %in% c("R3", "R3m", "R5")) {                                                 # 切点以上比例，报告到 0.1%
      if (use_mis) { out$p_cut <- round(ct$p_mis, round_dp + 2); out$cut <- ct$mis_cut }
      else         { out$p_cut <- round(ct$p_cut, round_dp + 2); out$cut <- ct$cut }
    }
    out
  }
  vis <- lapply(COVS, function(v) vis_cont(tb[[v]])); names(vis) <- COVS
  vis$n <- tb$n; vis$R <- R
  vis
}

# ---- 真值 ---------------------------------------------------------------------
#' 分块计算边际生存函数与风险（大样本时避免 时间 × 个体 的巨矩阵）
marginal_surv_chunked <- function(X, trt, shape, b_em, tt, pars = PARS, chunk = 2e4) {
  n <- nrow(X); S <- numeric(length(tt)); Hn <- numeric(length(tt))
  for (i in seq(1, n, by = chunk)) {
    idx <- i:min(i + chunk - 1, n)
    sm <- surv_matrices(X[idx, ], trt, shape, b_em, tt, pars)
    S <- S + rowSums(sm$S); Hn <- Hn + rowSums(sm$H)
  }
  list(S = S / n, h = Hn / S)
}

#' 在给定的协变量集合（权重相等）上计算 A vs B 的边际估计目标
#' 锚定生存：两个 vs-C 的 Cox 极限之差（共同删失分布）；另给直接 A vs B 极限（仅记录，用于局限性说明）
estimands_on <- function(X, shape, b_em, endpoint, cens_rate = NULL, pars = PARS,
                         tt = seq(0.05, pars$admin_cens, by = 0.05), tau = 24) {
  w <- rep(1, nrow(X))
  switch(endpoint,
    binary     = { ab <- marginal_binary(X, w, "A", "B", shape, b_em, pars); c(logOR = unname(ab["logOR"]), RD = unname(ab["RD"])) },
    continuous = c(meandiff = marginal_meandiff(X, w, "A", "B", shape, b_em, pars)),
    tte        = {
      G <- G_cens(tt, cens_rate, pars$admin_cens)
      sA <- marginal_surv_chunked(X, "A", shape, b_em, tt, pars); sB <- marginal_surv_chunked(X, "B", shape, b_em, tt, pars)
      sC <- marginal_surv_chunked(X, "C", shape, b_em, tt, pars)
      c(logHR_diff   = cox_limit(sA$S, sA$h, sC$S, sC$h, G) - cox_limit(sB$S, sB$h, sC$S, sC$h, G),
        logHR_direct = cox_limit(sA$S, sA$h, sB$S, sB$h, G),
        RMST_diff    = rmst(sA$S, tt, tau) - rmst(sB$S, tt, tau))
    })
}

#' 样本真值：在实际生成的 BC 样本（两臂合并的协变量）上算 A vs B 的边际效应
sample_truth <- function(ipd_bc, shape, b_em, endpoint, cens_rate = NULL, pars = PARS)
  estimands_on(ipd_bc[, COVS], shape, b_em, endpoint, cens_rate, pars)

#' DGM（超总体）真值：在 BC 总体上大样本分批计算，给出均值与蒙特卡洛标准误
#' @param spec dgm_spec() 的输出；@param n 总抽样量（二分类默认 2e6；生存默认 2e5，因每批要解 Cox 极限）
dgm_truth <- function(spec, shape, b_em, endpoint, cens_rate = NULL, n = NULL, batches = 10, seed = 20261010, pars = NULL) {
  if (is.null(pars)) pars <- pars_for(spec)
  if (is.null(n)) n <- if (endpoint == "tte") 2e5 else 2e6
  set.seed(seed)
  est <- t(sapply(seq_len(batches), function(b) {
    X <- population_draw(spec$bc, n %/% batches)
    estimands_on(X, shape, b_em, endpoint, cens_rate, pars)
  }))
  list(truth = colMeans(est), mcse = apply(est, 2, sd) / sqrt(batches), n = n, batches = batches)
}

# ---- 一次重复 -----------------------------------------------------------------
#' 生成一次重复：AC 试验 IPD、BC 试验 IPD（内部保留，用于样本真值）、BC 发表物
#' @param spec dgm_spec() 的输出；@param anchored FALSE 时 AC 试验只有 A 臂，BC 试验只有 B 臂
#' @param mis_cut R3m 的错位切点（默认 60 = 真实切点 65 − 0.5 SD）
gen_replicate <- function(spec, shape, b_em, endpoint, cens_rate = NULL, anchored = TRUE,
                          N_AC = 300, N_BC = 300, mis_cut = 60, pars = NULL) {
  if (is.null(pars)) pars <- pars_for(spec)
  ipd_ac <- gen_trial(spec$ac, N_AC, if (anchored) c("A", "C") else "A", shape, b_em, endpoint, cens_rate, pars)
  ipd_bc <- gen_trial(spec$bc, N_BC, if (anchored) c("B", "C") else "B", shape, b_em, endpoint, cens_rate, pars)
  if (anchored) ipd_bc$trt <- factor(ipd_bc$trt, levels = c("C", "B"))   # 让 glm/coxph 的系数是 B vs C
  list(ipd_ac = ipd_ac, ipd_bc = ipd_bc,
       pub_bc = publish_bc(ipd_bc, endpoint, pars$cut, mis_cut, unanchored = !anchored),
       truth_sample = sample_truth(ipd_bc, shape, b_em, endpoint, cens_rate, pars))
}
