# 样本层面的 MAIC 族估计量（protocol v0.2 第 5 节：M1 全部精确、M1m 只用均值、M7 sMAIC、O1/O2 oracle）
# 流程：可见的基线表 → 约束矩阵 → 指数倾斜权重（含不可行时的回退规则）→ AC 试验的加权 A vs C 对比
#       → 与发表的 B vs C 相减（锚定）→ 三种方差：(a) 条件于目标、含权重估计的 sandwich；
#       (b) bootstrap + 独立重抽目标；(c) bootstrap + 联合重抽目标（从重建的边际分布抽伪样本）。
# 依赖：R/01–05（由调用脚本 source）。可打包代码：不用 setwd()，不写绝对路径。

# ---- 1. 由可见的基线表构造多协变量约束 -------------------------------------------
#' 把一个连续协变量的可见汇总量翻译成 (H, targets, tol, se)；range 返回截取掩码
#' @param use 用到的汇总量名（"mean","sd","quantiles","median","prop","range"）；@param soft 其中软匹配的
cov_constraints <- function(x, vis, use, soft = character(0), kappa = 1.96, round_dp = 1) {
  H <- NULL; tg <- NULL; se <- NULL; nm <- NULL; r_half <- 10^(-round_dp) / 2; n <- vis$n
  add <- function(h, t, s, name) { H <<- cbind(H, h); tg <<- c(tg, t); se <<- c(se, s); nm <<- c(nm, name) }
  if ("mean" %in% use && !is.null(vis$mean)) add(x, vis$mean, vis$sd / sqrt(n), "mean")
  if ("sd" %in% use && !is.null(vis$sd))                                    # SD → 二阶原点矩
    add(x^2, vis$sd^2 + vis$mean^2, vis$sd^2 * sqrt(2 / (n - 1)), "sd")
  if ("quantiles" %in% use && !is.null(vis$q))                              # 可见的分位数 → 指示函数
    for (j in which(!is.na(vis$q))) { p <- c(0.25, 0.5, 0.75)[j]; add(as.numeric(x <= vis$q[j]), p, sqrt(p * (1 - p) / n), paste0("q", j)) }
  if ("median" %in% use && !is.null(vis$q) && !is.na(vis$q[2])) add(as.numeric(x <= vis$q[2]), 0.5, sqrt(0.25 / n), "q2")
  if ("prop" %in% use && !is.null(vis$p_cut)) { p <- vis$p_cut; add(as.numeric(x > vis$cut), p, sqrt(p * (1 - p) / n), "prop") }
  tol <- ifelse(nm %in% c("mean", "sd", "q1", "q2", "q3", "prop") & nm %in% soft_names(soft), kappa * se + r_half * c(mean = 1, sd = 2 * (if (is.null(vis$sd)) 1 else vis$sd), q1 = 0, q2 = 0, q3 = 0, prop = 0.01)[nm], 0)
  mask <- if ("range" %in% use && !is.null(vis$range)) x >= vis$range[1] & x <= vis$range[2] else rep(TRUE, length(x))
  list(H = H, targets = tg, tol = unname(tol), se = se, names = nm, mask = mask)
}
soft_names <- function(soft) unlist(list(sd = "sd", quantiles = c("q1", "q2", "q3"), mean = "mean", prop = "prop")[soft])

#' 多协变量约束：三个连续协变量按同一配方翻译，再把约束矩阵拼起来（范围掩码取交集）
#' @param recipe list(use = ..., soft = ..., convert = NULL/"iqr"/"range")：convert 时把中位数换算成均值/SD
build_constraints_multi <- function(X, vis, recipe, kappa = 1.96, round_dp = 1) {
  parts <- lapply(COVS, function(v) {
    vc <- vis[[v]]
    if (!is.null(recipe$convert) && is.null(vc$mean)) {                     # 现行做法：没报均值时把中位数换算成均值/SD
      cv <- if (!is.null(vc$q) && !any(is.na(vc$q)))
              estmeansd::qe.mean.sd(q1.val = vc$q[1], med.val = vc$q[2], q3.val = vc$q[3], n = vc$n)
            else if (!is.null(vc$range) && !is.null(vc$q) && !is.na(vc$q[2]))
              estmeansd::qe.mean.sd(min.val = vc$range[1], med.val = vc$q[2], max.val = vc$range[2], n = vc$n)
            else NULL
      if (!is.null(cv)) { vc$mean <- cv$est.mean; vc$sd <- cv$est.sd }
    }
    cov_constraints(X[[v]], vc, recipe$use, recipe$soft, kappa, round_dp)
  })
  names(parts) <- COVS
  H <- do.call(cbind, lapply(parts, `[[`, "H")); tg <- unlist(lapply(parts, `[[`, "targets"))
  tol <- unlist(lapply(parts, `[[`, "tol")); se <- unlist(lapply(parts, `[[`, "se"))
  nm <- unlist(lapply(COVS, function(v) paste0(v, "_", parts[[v]]$names)))
  mask <- Reduce(`&`, lapply(parts, `[[`, "mask"))
  list(H = H, targets = unname(tg), tol = unname(tol), se = unname(se), names = nm, mask = mask)
}

#' 提取策略配方（样本层面）。use 中的项只在基线表可见时生效
MAIC_RECIPES <- list(
  M1  = list(use = c("mean", "sd", "quantiles", "prop", "range")),                     # 全部可得汇总量精确匹配
  M1m = list(use = "mean", convert = "iqr"),                                          # 只用均值（中位数换算）
  M1c = list(use = c("mean", "sd"), convert = "iqr"),                                  # 换算后的均值 + SD
  M1q = list(use = c("mean", "quantiles", "prop", "range")),                          # 决策 24：均值 + 分位数 + n(%)（+ 范围截取），表上有 SD 也不用
  M7  = list(use = c("mean", "sd", "quantiles", "prop", "range"), soft = c("sd", "quantiles")),  # sMAIC（κ = 1.96）
  M7k1  = list(use = c("mean", "sd", "quantiles", "prop", "range"), soft = c("sd", "quantiles"), kappa = 1),    # sMAIC，κ = 1（决策 24 的容差敏感性）
  M7k05 = list(use = c("mean", "sd", "quantiles", "prop", "range"), soft = c("sd", "quantiles"), kappa = 0.5)   # sMAIC，κ = 0.5
)

# ---- 2. 权重与回退规则 ---------------------------------------------------------
#' 按配方拟合权重；M1 不可行时依次去掉范围 → n(%) → 分位数 → SD（记录回退层级与不可行类型）
#' 2026-10-10（决策 20）：n(%) 排在分位数之前。主运行诊断发现 R5 下 13–22% 的不可行几乎全是"n(%) 的切点（60）落在
#' 某个四分位数的舍入距离内"造成的自相矛盾（P(X ≤ 60.0) = 0.75 与 P(X > 60) = 0.247 不能同时成立）；
#' 先丢掉这 1 个冗余的 n(%) 比丢掉 3 个分位数损失的信息少。fallback 字段：1 = 去范围，2 = 再去 n(%)，3 = 再去分位数，4 = 再去 SD
fit_maic_weights <- function(X, vis, recipe, kappa = 1.96, round_dp = 1) {
  if (!is.null(recipe$kappa)) kappa <- recipe$kappa                          # 配方自带 κ（sMAIC 容差敏感性）时覆盖默认值
  drop_order <- list(character(0), "range", c("range", "prop"), c("range", "prop", "quantiles"), c("range", "prop", "quantiles", "sd"))
  for (lvl in seq_along(drop_order)) {
    rc <- recipe; rc$use <- setdiff(rc$use, drop_order[[lvl]])
    cs <- build_constraints_multi(X, vis, rc, kappa, round_dp)
    if (is.null(cs$H)) break                                                # 回退后一个约束都不剩（如 R2r 下只有中位数 + 范围）→ 判不可行
    base <- as.numeric(cs$mask); if (sum(base) < 10) next
    tw <- tilt_weights(cs$H, cs$targets, cs$tol, base = base / sum(base))
    if (tw$feasible) {
      w <- tw$w / sum(tw$w) * nrow(X)                                       # 归一化到和为 n（截掉的个体权重为 0）
      return(list(w = w, ess = ess(w), feasible = TRUE, fallback = lvl - 1L, constraints = cs, lambda = tw$lambda,
                  infeasible_type = if (lvl == 1) "none" else "inconsistent_or_hull"))
    }
    if (length(recipe$use) == 1) break                                      # 只用均值的配方没有回退余地
  }
  list(w = rep(1, nrow(X)), ess = nrow(X), feasible = FALSE, fallback = NA_integer_, constraints = cs,
       lambda = NULL, infeasible_type = "out_of_hull")
}

# ---- 3. AC 试验的加权 A vs C 对比 ------------------------------------------------
#' @return c(est, se_fixedw)：se_fixedw 把权重当已知（Hájek 线性化 / 加权 Cox 的稳健 SE）
weighted_contrast <- function(ipd, w, endpoint, arms = c("A", "C"), tau = 24) {
  a <- ipd$trt == arms[1]; c0 <- ipd$trt == arms[2]
  switch(endpoint,
    binary = {
      pA <- sum(w[a] * ipd$y[a]) / sum(w[a]); pC <- sum(w[c0] * ipd$y[c0]) / sum(w[c0])
      vA <- sum((w[a] * (ipd$y[a] - pA))^2) / sum(w[a])^2; vC <- sum((w[c0] * (ipd$y[c0] - pC))^2) / sum(w[c0])^2
      c(est = qlogis(pA) - qlogis(pC), se = sqrt(vA / (pA * (1 - pA))^2 + vC / (pC * (1 - pC))^2),
        RD = pA - pC, RD_se = sqrt(vA + vC))
    },
    continuous = {
      mA <- sum(w[a] * ipd$y[a]) / sum(w[a]); mC <- sum(w[c0] * ipd$y[c0]) / sum(w[c0])
      vA <- sum((w[a] * (ipd$y[a] - mA))^2) / sum(w[a])^2; vC <- sum((w[c0] * (ipd$y[c0] - mC))^2) / sum(w[c0])^2
      c(est = mA - mC, se = sqrt(vA + vC))
    },
    tte = {
      d <- ipd[w > 0, ]; ww <- w[w > 0]; d$trt <- factor(d$trt, levels = arms[2:1])   # 系数 = arms[1] vs arms[2]
      f <- survival::coxph(survival::Surv(time, status) ~ trt, data = d, weights = ww, robust = TRUE)
      km <- survival::survfit(survival::Surv(time, status) ~ trt, data = d, weights = ww)
      tb <- summary(km, rmean = tau)$table
      c(est = unname(coef(f)[1]), se = unname(sqrt(f$var[1, 1])),
        RMST = unname(tb[2, "rmean"] - tb[1, "rmean"]), RMST_se = unname(sqrt(tb[2, "se(rmean)"]^2 + tb[1, "se(rmean)"]^2)))
    })
}

# ---- 4. 方差 (a)：条件于目标、含权重估计的 M-估计 sandwich（二分类 / 连续） -----------
#' 参数 θ = (λ, μ_A, μ_C)；估计方程：Σ w_i(λ)(h_i − s) = 0，Σ w_i T_A (y − μ_A) = 0，Σ w_i T_C (y − μ_C) = 0
#' 软约束中起作用的项按其达到的值当作精确约束，不起作用的项（λ_k = 0）剔除
sandwich_maic <- function(ipd, fit, endpoint, arms = c("A", "C")) {
  if (!endpoint %in% c("binary", "continuous") || is.null(fit$lambda)) return(NA_real_)
  cs <- fit$constraints; keep <- fit$constraints$mask
  H <- cs$H[keep, , drop = FALSE]; y <- ipd$y[keep]; tA <- ipd$trt[keep] == arms[1]; tC <- ipd$trt[keep] == arms[2]
  act <- abs(fit$lambda) > 1e-8 | cs$tol == 0                               # 精确约束 + 起作用的软约束
  if (!any(act)) return(NA_real_)
  Hs <- H[, act, drop = FALSE]; s_eff <- as.numeric(crossprod(Hs, fit$w[keep])) / sum(fit$w[keep])   # 达到的值
  Z <- sweep(Hs, 2, s_eff); K <- ncol(Z); n <- nrow(Z)
  wfun <- function(l) exp(as.numeric(Z %*% l))
  score_mean <- function(th) { w <- wfun(th[1:K]); muA <- th[K + 1]; muC <- th[K + 2]
    c(colMeans(w * Z), mean(w * tA * (y - muA)), mean(w * tC * (y - muC))) }
  w0 <- fit$w[keep]; muA <- sum(w0 * tA * y) / sum(w0 * tA); muC <- sum(w0 * tC * y) / sum(w0 * tC)
  th <- c(fit$lambda[act], muA, muC)
  A <- -numDeriv::jacobian(score_mean, th)
  w <- wfun(th[1:K]); Psi <- cbind(w * Z, w * tA * (y - muA), w * tC * (y - muC))
  B <- crossprod(Psi) / n
  V <- tryCatch(solve(A) %*% B %*% t(solve(A)) / n, error = function(e) matrix(NA, K + 2, K + 2))
  g <- c(rep(0, K), if (endpoint == "binary") c(1 / (muA * (1 - muA)), -1 / (muC * (1 - muC))) else c(1, -1))
  sqrt(as.numeric(t(g) %*% V %*% g) * n / (n - K - 2))                    # HC1 型小样本校正（自由度）
}

#' 无锚定：单臂（A）的加权结局汇总；二分类给比例，生存给 RMST(τ)
weighted_arm <- function(ipd, w, endpoint, arm = "A", tau = 24) {
  a <- ipd$trt == arm
  switch(endpoint,
    binary = { p <- sum(w[a] * ipd$y[a]) / sum(w[a]); c(est = p, se = sqrt(sum((w[a] * (ipd$y[a] - p))^2)) / sum(w[a])) },
    continuous = { m <- sum(w[a] * ipd$y[a]) / sum(w[a]); c(est = m, se = sqrt(sum((w[a] * (ipd$y[a] - m))^2)) / sum(w[a])) },
    tte = { d <- ipd[a & w > 0, ]; ww <- w[a & w > 0]
      km <- survival::survfit(survival::Surv(time, status) ~ 1, data = d, weights = ww); tb <- summary(km, rmean = tau)$table
      c(est = unname(tb["rmean"]), se = unname(tb["se(rmean)"])) })
}

#' 加权 KM 的 RMST(τ) 及其无穷小刀切（IJ）标准误，权重视为固定（2026-10-13 决策 39）
#' 原因：weighted_arm() 取的 summary(survfit(weights = w))$table 的 se(rmean) 把权重当频数权重，不随 ESS 下降而增大
#'（一次重复上 ESS 177：0.56，而 IJ 0.82、权重固定的 bootstrap 0.83；不加权时 IJ 与 Greenwood 公式相同）。
#' IJ：survfit(influence = TRUE) 给出每个个体对 S(t) 的影响值，在 [0, τ] 的阶梯上积分得 RMST 的影响值，平方和开根即标准误。
rmst_ij <- function(time, status, w, tau = 24) {
  km <- survival::survfit(survival::Surv(time, status) ~ 1, weights = w, influence = TRUE)
  keep <- km$time <= tau
  grid <- c(0, km$time[keep], tau); dt <- diff(grid)                    # 每段 [grid_k, grid_k+1) 上 S 为常数；第一段 S = 1
  infl <- km$influence.surv; if (is.null(infl)) infl <- km$influence
  IF <- as.numeric(as.matrix(infl)[, keep, drop = FALSE] %*% dt[-1])
  c(est = sum(c(1, km$surv[keep]) * dt), se = sqrt(sum(IF^2)))
}

# ---- 5. 方差 (b)(c)：bootstrap，独立 / 联合重抽目标 ------------------------------
#' 独立重抽：每个可见汇总量各自加 N(0, SE²) 噪声
redraw_table_indep <- function(vis) {
  for (v in COVS) { vc <- vis[[v]]; n <- vc$n
    if (!is.null(vc$mean)) vc$mean <- vc$mean + rnorm(1, 0, vc$sd / sqrt(n))
    if (!is.null(vc$sd))   vc$sd   <- sqrt(max(vc$sd^2 + rnorm(1, 0, vc$sd^2 * sqrt(2 / (n - 1))), 1e-6))
    if (!is.null(vc$q))    { ok <- !is.na(vc$q); vc$q[ok] <- vc$q[ok] + rnorm(sum(ok), 0, 1.25 * (if (!is.null(vc$sd)) vc$sd else diff(vc$q[c(1, 3)]) / 1.35) / sqrt(n)) }
    if (!is.null(vc$p_cut)) vc$p_cut <- min(max(vc$p_cut + rnorm(1, 0, sqrt(vc$p_cut * (1 - vc$p_cut) / n)), 0.001), 0.999)
    vis[[v]] <- vc }
  vis
}
#' 联合重抽：从按可见汇总量重建的边际分布抽一个 N_BC 的伪样本，重新算整张表（各汇总量一起抖动）
#' @param Qs 预先重建好的分位数函数 list(x1 = Q, x2 = Q, x3 = Q)（每次重复只拟合一次，bootstrap 里反复用）
redraw_table_joint <- function(vis, Qs, round_dp = 1, cut = PARS$cut) {
  n <- vis$n
  pseudo <- as.data.frame(lapply(Qs[COVS], function(Q) Q(runif(n))))      # 伪样本的三个协变量独立抽（边际重抽）
  tb <- baseline_table(pseudo, cut = cut, mis_cut = if (!is.null(vis$x1$cut)) vis$x1$cut else cut)   # R3m 时 vis 里的切点是错位切点
  visible_table(tb, vis$R, round_dp)
}

#' bootstrap：重抽 AC（按臂分层）+ 按 mode 重抽目标；返回 A vs C 估计的 SD
boot_maic <- function(ipd, vis, recipe, endpoint, B = 200, mode = c("fixed", "indep", "joint"), kappa = 1.96, round_dp = 1) {
  mode <- match.arg(mode); X <- ipd[, COVS]
  idx_by_arm <- split(seq_len(nrow(ipd)), ipd$trt)
  Qs <- if (mode == "joint") { q <- lapply(COVS, function(v) reconstruct_marginal(vis[[v]])); names(q) <- COVS; q } else NULL
  ests <- vapply(seq_len(B), function(b) {
    idx <- unlist(lapply(idx_by_arm, function(i) sample(i, replace = TRUE)))
    vb <- switch(mode, fixed = vis, indep = redraw_table_indep(vis), joint = redraw_table_joint(vis, Qs, round_dp))
    fb <- fit_maic_weights(X[idx, ], vb, recipe, kappa, round_dp)
    if (!fb$feasible) return(NA_real_)
    weighted_contrast(ipd[idx, ], fb$w, endpoint)[["est"]]
  }, numeric(1))
  c(se = sd(ests, na.rm = TRUE), fail = mean(is.na(ests)))
}

# ---- 6. 一次重复里跑一个 MAIC 族方法 ----------------------------------------------
#' 锚定：A vs C（加权）减 B vs C（发表）；无锚定：A 臂加权结局 vs B 臂发表结局
#' @return 一行 data.frame：主估计（est_ab）、次要估计（est_sec：风险差 / RMST）、四种 SE、ESS、可行性、回退层级
run_maic_method <- function(rep, vis, method, endpoint, B = 200, kappa = 1.96, round_dp = 1, do_boot = TRUE, anchored = TRUE) {
  X <- rep$ipd_ac[, COVS]
  fit <- if (method == "M0")                                                 # M0：未调整（Bucher / 朴素比较），权重全为 1
    list(w = rep(1, nrow(X)), ess = nrow(X), feasible = TRUE, fallback = 0L, infeasible_type = "none", lambda = NULL)
  else fit_maic_weights(X, vis, MAIC_RECIPES[[method]], kappa, round_dp)
  se_pub <- rep$pub_bc$effect[["se"]]; pub <- rep$pub_bc$effect[["est"]]
  if (anchored) {
    wc <- weighted_contrast(rep$ipd_ac, fit$w, endpoint)
    est_ac <- wc[["est"]]; est_ab <- est_ac - pub; se_fw <- wc[["se"]]
    se_a <- if (endpoint %in% c("binary", "continuous") && method != "M0") sandwich_maic(rep$ipd_ac, fit, endpoint) else wc[["se"]]
    sec_name <- switch(endpoint, binary = "RD", tte = "RMST", NA); ps <- rep$pub_bc$effect_sec
    est_sec <- if (!is.na(sec_name)) wc[[sec_name]] - ps[["est"]] else NA_real_
    se_sec  <- if (!is.na(sec_name)) sqrt(wc[[paste0(sec_name, "_se")]]^2 + ps[["se"]]^2) else NA_real_
  } else {                                                                  # 无锚定：A 臂 vs 发表的 B 臂
    wa <- weighted_arm(rep$ipd_ac, fit$w, endpoint)
    if (endpoint == "binary") {                                             # 主：logOR；次：风险差
      pA <- wa[["est"]]; pB <- pub
      est_ac <- qlogis(pA); est_ab <- qlogis(pA) - qlogis(pB)
      se_fw <- sqrt((wa[["se"]] / (pA * (1 - pA)))^2 + (se_pub / (pB * (1 - pB)))^2)
      est_sec <- pA - pB; se_sec <- sqrt(wa[["se"]]^2 + se_pub^2); se_pub <- 0   # 发表 SE 已并入 se_fw
    } else {                                                                # 生存：主 RMST 差；连续：均值差
      est_ac <- wa[["est"]]; est_ab <- est_ac - pub; se_fw <- sqrt(wa[["se"]]^2 + se_pub^2); se_pub <- 0
      est_sec <- NA_real_; se_sec <- NA_real_
    }
    se_a <- NA_real_
  }
  out <- data.frame(method = method, est_ac = est_ac, est_ab = est_ab, est_sec = est_sec, se_sec = se_sec, ess = fit$ess,
                    feasible = fit$feasible, fallback = fit$fallback, infeasible_type = fit$infeasible_type,
                    se_fixedw = if (anchored) sqrt(se_fw^2 + se_pub^2) else se_fw,
                    se_sandwich = if (anchored) sqrt(se_a^2 + se_pub^2) else NA_real_,
                    se_boot_indep = NA_real_, se_boot_joint = NA_real_, boot_fail = NA_real_)
  if (do_boot && fit$feasible && anchored && method != "M0") {
    bi <- boot_maic(rep$ipd_ac, vis, MAIC_RECIPES[[method]], endpoint, B, "indep", kappa, round_dp)
    bj <- boot_maic(rep$ipd_ac, vis, MAIC_RECIPES[[method]], endpoint, B, "joint", kappa, round_dp)
    out$se_boot_indep <- sqrt(bi[["se"]]^2 + se_pub^2); out$se_boot_joint <- sqrt(bj[["se"]]^2 + se_pub^2)
    out$boot_fail <- max(bi[["fail"]], bj[["fail"]])
  }
  out
}
