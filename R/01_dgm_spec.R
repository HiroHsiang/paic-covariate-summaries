# 协变量分布规格（protocol v1.2 第 3.1–3.2 节；决策 18、19）。
# 三个连续协变量 X1–X3，抽象的"评分型"尺度，不赋临床含义；三者边际分布相同，用高斯 copula 相关（rho）。
# 两个试验共享同一分布族（同一疾病里同一指标的"生物学形状"相同）；BC 相对 AC 只有
#   位置（+0.3 SD）、尺度（×1.2）、以及入排标准造成的截断（D4）这三类差别。
#   D1：两边正态；D2：两边右偏（偏度 +0.8）；D3：两边左偏（−0.8）；D4：AC 正态、BC 在 70（AC 均值 + 2 SD）处上截断。
# AC：均值 50、SD 10；BC：均值 53、SD 12（D4 为截断后均值 53、截断前 SD 12）。
# 可打包代码：不用 setwd()，不写绝对路径。

SCALE   <- list(mean = 50, sd = 10)                                    # 标准化基准：z = (x − 50) / 10（AC 的均值 / SD）
AC_BASE <- list(x = list(type = "normal", mean = 50, sd = 10), rho = 0.35)   # D1 / D4 的 AC 规格；也是 PARS 默认的中心化基准

#' 一个 DGM 的两试验协变量规格
#' @param D "D1"（正态）/"D2"（右偏）/"D3"（左偏）/"D4"（BC 上截断）
#' @param delta_sd    BC 均值相对 AC 的平移（以 AC 的 SD 为单位，默认 0.3 → 53）
#' @param scale_ratio BC 的 SD / AC 的 SD（默认 1.2 → 12）；D4 为截断前的 SD。
#'   取 1.2 而不是 1.3：三个协变量同时平移、放大后，总体重叠是三维的，1.3 时 D2 的 oracle ESS 只剩 15%，
#'   1.2 时 D1/D2/D4 分别约 0.56/0.43/0.50，与 v1.1 单变量设计的重叠程度相当（见 sims/01 表 1b）
#' @param skew        D2 / D3 的偏度（两试验相同；D2 取 +skew，D3 取 −skew）
#' @param comp_ratio  混合成分 SD 与总 SD 之比（默认 0.75 → 左成分占 90.5%，右成分在 +2 SD 处形成"肩部"；
#'   NULL 时由 choose_comp_ratio() 按 oracle ESS 自动选，但那会选出 1–3% 的远端小簇，形状不自然，故固定）
#' @param trunc_sd    D4 的截断点（AC 均值 + trunc_sd × AC 的 SD，默认 70）
#' @param rho         高斯 copula 相关系数（两试验相同）
#' @return list(D, ac = list(x, rho), bc = list(x, rho), ...)：x 是一个边际分布规格，三个协变量共用
dgm_spec <- function(D = c("D1", "D2", "D3", "D4"), delta_sd = 0.3, scale_ratio = 1.2, skew = 0.8,
                     comp_ratio = 0.75, trunc_sd = 2, rho = 0.35) {
  D <- match.arg(D)
  m_ac <- SCALE$mean; s_ac <- SCALE$sd
  m_bc <- m_ac + delta_sd * s_ac; s_bc <- s_ac * scale_ratio; trunc_at <- m_ac + trunc_sd * s_ac
  sk <- switch(D, D1 = 0, D2 = +abs(skew), D3 = -abs(skew), D4 = 0)
  if (D %in% c("D2", "D3")) {
    if (is.null(comp_ratio)) comp_ratio <- choose_comp_ratio(sk, delta_sd, scale_ratio)
    ac_x <- mix_spec(m_ac, s_ac, sk, comp_ratio)
    bc_x <- mix_spec(m_bc, s_bc, sk, comp_ratio)          # 同一成分比例与偏度 → BC 恰是 AC 的位置–尺度变换
  } else {
    ac_x <- list(type = "normal", mean = m_ac, sd = s_ac)
    bc_x <- if (D == "D1") list(type = "normal", mean = m_bc, sd = s_bc)
            else list(type = "truncnorm", mean = trunc_mean_pre(m_bc, s_bc, trunc_at), sd = s_bc, upper = trunc_at,
                      mean_target = m_bc)                  # 截断前均值经反解，使截断后的均值仍为 53（与 D1 同位置）
  }
  list(D = D, ac = list(x = ac_x, rho = rho), bc = list(x = bc_x, rho = rho),
       delta_sd = delta_sd, scale_ratio = scale_ratio, skew = sk, comp_ratio = comp_ratio,
       trunc_at = if (D == "D4") trunc_at else NA)
}

#' 两成分等方差正态混合的规格（由均值 / SD / 偏度 / 成分 SD 比例确定）
mix_spec <- function(mean, sd, skew, comp_ratio)
  c(list(type = "mixture", skew_target = skew, mean = mean, sd = sd, comp_ratio = comp_ratio),
    solve_mixture(mean, sd, skew, comp_ratio * sd))

#' 上截断正态：解截断前的均值 μ*，使 N(μ*, sd²) 在 upper 处截断后的均值等于 target
trunc_mean_pre <- function(target, sd, upper) {
  tmean <- function(mu) { a <- (upper - mu) / sd; mu - sd * dnorm(a) / pnorm(a) }   # 上截断正态的均值
  uniroot(function(mu) tmean(mu) - target, c(target, target + 4 * sd), tol = 1e-8)$root
}

# ---- 两成分等方差正态混合：由均值 / SD / 偏度解参数 ------------------------
# 记 d = m2 − m1，p 为左成分的权重，则
#   均值 = p·m1 + (1−p)·m2；方差 = s² + p(1−p)d²；三阶中心矩 = d³·p(1−p)(2p−1)
# p > 0.5 且 d > 0 时为右偏（大成分在左、小成分拖出右尾）；左偏取镜像。
# 解出的 p 只依赖 sd_comp / sd 与偏度，所以同比例、同偏度的两个混合互为位置–尺度变换。
solve_mixture <- function(mean, sd, skew, sd_comp) {
  stopifnot(sd_comp < sd)
  d_of <- function(p) sqrt((sd^2 - sd_comp^2) / (p * (1 - p)))
  f <- function(p) d_of(p)^3 * p * (1 - p) * (2 * p - 1) / sd^3 - abs(skew)
  p <- uniroot(f, c(0.5 + 1e-6, 0.995), tol = 1e-10)$root
  d <- d_of(p)
  if (skew >= 0) list(p = p, m1 = mean - (1 - p) * d, m2 = mean + p * d, s = sd_comp)
  else           list(p = p, m1 = mean + (1 - p) * d, m2 = mean - p * d, s = sd_comp)
}

#' 在成分 SD 比例的候选值中选使"AC → BC"oracle 权重的 ESS 最大的那个（一维数值积分）
#' 比例越小成分分得越开（尾部越重、权重越极端）；比例越大越接近正态，但同样偏度需要更小的远端成分
choose_comp_ratio <- function(skew, delta_sd = 0.3, scale_ratio = 1.2, grid = seq(0.35, 0.95, by = 0.025),
                              x = seq(-20, 140, by = 0.01)) {
  best <- NULL
  for (r in grid) {
    f_ac <- dmix(x, mix_spec(SCALE$mean, SCALE$sd, skew, r)); f_ac <- f_ac / sum(f_ac)
    f_bc <- dmix(x, mix_spec(SCALE$mean + delta_sd * SCALE$sd, SCALE$sd * scale_ratio, skew, r)); f_bc <- f_bc / sum(f_bc)
    ess <- 1 / sum(f_bc^2 / f_ac)                      # oracle 权重的 ESS/N 极限（一维）
    if (is.null(best) || ess > best$ess) best <- list(r = r, ess = ess)
  }
  best$r
}

dmix <- function(x, m) m$p * dnorm(x, m$m1, m$s) + (1 - m$p) * dnorm(x, m$m2, m$s)
pmix <- function(x, m) m$p * pnorm(x, m$m1, m$s) + (1 - m$p) * pnorm(x, m$m2, m$s)
qmix <- function(prob, m) {                       # 分位数：细网格上反插值
  xs <- seq(m$m1 - 8 * m$s - abs(m$m2 - m$m1), m$m2 + 8 * m$s + abs(m$m2 - m$m1), length.out = 20001)
  F <- pmix(xs, m); keep <- c(TRUE, diff(F) > 0)
  approx(F[keep], xs[keep], xout = prob, rule = 2)$y
}

# ---- 统一的边际分布接口：d（密度）、p（分布函数）、q（分位数）、r（随机数） ----
marg_dist <- function(spec_x) {
  switch(spec_x$type,
    normal = list(
      d = function(x) dnorm(x, spec_x$mean, spec_x$sd),
      p = function(x) pnorm(x, spec_x$mean, spec_x$sd),
      q = function(p) qnorm(p, spec_x$mean, spec_x$sd),
      r = function(n) rnorm(n, spec_x$mean, spec_x$sd)),
    mixture = list(
      d = function(x) dmix(x, spec_x),
      p = function(x) pmix(x, spec_x),
      q = function(p) qmix(p, spec_x),
      r = function(n) { k <- rbinom(n, 1, spec_x$p)
        rnorm(n, ifelse(k == 1, spec_x$m1, spec_x$m2), spec_x$s) }),
    truncnorm = {                                   # 上截断正态
      Fu <- pnorm(spec_x$upper, spec_x$mean, spec_x$sd)
      list(
        d = function(x) ifelse(x <= spec_x$upper, dnorm(x, spec_x$mean, spec_x$sd) / Fu, 0),
        p = function(x) pmin(pnorm(x, spec_x$mean, spec_x$sd) / Fu, 1),
        q = function(p) qnorm(p * Fu, spec_x$mean, spec_x$sd),
        r = function(n) qnorm(runif(n) * Fu, spec_x$mean, spec_x$sd))
    })
}

#' 边际分布的矩与分位数（数值积分），用于核对与报告
marg_moments <- function(spec_x, cut = 65, x = seq(-20, 140, by = 0.01)) {
  f <- marg_dist(spec_x)$d(x); f <- f / sum(f)
  m <- sum(x * f); s <- sqrt(sum((x - m)^2 * f))
  list(mean = m, sd = s, skew = sum(((x - m) / s)^3 * f),
       q = marg_dist(spec_x)$q(c(0.25, 0.5, 0.75)), p_cut = sum(f[x > cut]))
}

# ---- 从 copula 生成一个（大）总体样本 ------------------------------------------
#' @param spec_trial list(x = 边际规格, rho)（即 dgm_spec()$ac 或 $bc）
#' @param n 样本量；@param seed 种子（总体真值用固定种子）
#' @return data.frame(x1, x2, x3)：三个协变量同一边际，经高斯 copula 相关
population_draw <- function(spec_trial, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  R <- matrix(spec_trial$rho, 3, 3); diag(R) <- 1
  Z <- matrix(rnorm(n * 3), n, 3) %*% chol(R)              # 相关的标准正态
  Q <- marg_dist(spec_trial$x)$q
  data.frame(x1 = Q(pnorm(Z[, 1])), x2 = Q(pnorm(Z[, 2])), x3 = Q(pnorm(Z[, 3])))   # 经 copula 变换到边际分布
}

# ---- 三维密度比（oracle 权重）与其 ESS：检查"AC → BC"的权重是否有界 ----------------
#' f_BC(x) / f_AC(x)，x 为 data.frame(x1, x2, x3)：边际密度比的乘积 × 高斯 copula 密度比
density_ratio_3d <- function(spec, X) {
  dac <- marg_dist(spec$ac$x); dbc <- marg_dist(spec$bc$x)
  R <- matrix(spec$ac$rho, 3, 3); diag(R) <- 1; A <- solve(R) - diag(3)
  clamp <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)
  lr <- 0; Qac <- NULL; Qbc <- NULL
  for (v in c("x1", "x2", "x3")) {
    x <- X[[v]]
    lr <- lr + log(dbc$d(x)) - log(dac$d(x))
    Qac <- cbind(Qac, qnorm(clamp(dac$p(x)))); Qbc <- cbind(Qbc, qnorm(clamp(dbc$p(x))))
  }
  lr <- lr - 0.5 * (rowSums((Qbc %*% A) * Qbc) - rowSums((Qac %*% A) * Qac))   # copula 密度比（两试验 rho 相同）
  w <- exp(lr); w[!is.finite(w)] <- 0; w
}

#' oracle 权重（真实密度比）在 AC 总体上的 ESS/N，以及 AC 中央区域内的最大权重
oracle_ess_3d <- function(spec, n = 2e5, seed = 20261001) {
  X <- population_draw(spec$ac, n, seed); w <- density_ratio_3d(spec, X)
  c(ess_oracle_3d = mean(w)^2 / mean(w^2), max_w_central = unname(quantile(w, 0.999)))
}
