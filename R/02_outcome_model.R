# 结局模型（protocol v1.2 第 3.3–3.4 节）：效应修饰形状、线性预测、二分类 / 生存 / 连续结局，
# 以及总体层面的边际效应（用于真值和校准）。可打包代码：不用 setwd()，不写绝对路径。
# 决策 18：三个协变量都是预后因子，也都是效应修饰因子；同一形状 g、同一强度 b_EM 作用于 X1、X2、X3。

# ---- 固定参数 ---------------------------------------------------------------
PARS <- list(
  beta        = -log(0.67),   # 每个（标准化）协变量的预后效应
  theta_A     = log(0.5),     # A vs C 的条件效应（EM 项取 AC 平均时）
  theta_B     = log(0.7),     # B vs C 的条件效应
  alpha       = qlogis(0.3),  # 二分类结局的截距（C 臂基准风险 30%）
  weib_shape  = 1.3,          # Weibull 形状参数
  weib_median = 12,           # 基线中位生存（月）
  admin_cens  = 36,           # 行政删失（月）
  cut         = 65,           # 阈值型 EM 的切点（AC 均值 + 1.5 SD；决策 21：离 BC 的四分位数（≈ 45 / 53 / 61）足够远，
                              #   避免报告的 n (%) 与舍入后的四分位数约束互相矛盾；60 时 25–43% 的样本会冲突）
  g_center    = NULL          # 各形状下 g 在 AC 人群的均值（由 pars_for() 按 DGM 填入；默认为正态 AC）
)

# ---- 效应修饰形状 g(X) ------------------------------------------------------
z_std <- function(x) (x - SCALE$mean) / SCALE$sd    # 标准化（以 AC 的均值 / SD 为准）

#' @param shape "linear" / "quad" / "threshold"
g_em <- function(x, shape, cut = PARS$cut) {
  switch(shape,
    linear    = z_std(x),
    quad      = z_std(x)^2,
    threshold = as.numeric(x > cut),
    stop("未知的 EM 形状：", shape))
}

#' g(X) 在 AC 人群上的均值（三种形状），用于把 EM 项中心化：
#' 中心化后 θ_A、θ_B 是"AC 平均病人"的条件效应，EM 只改变效应在病人之间的差异，不改变基准水平。
#' AC 的形状随 D 变（决策 19），所以按 AC 的边际分布数值积分
g_center_of <- function(spec_x, cut = PARS$cut, shapes = c("linear", "quad", "threshold"), x = seq(-20, 140, by = 0.01)) {
  f <- marg_dist(spec_x)$d(x); f <- f / sum(f)
  vapply(shapes, function(s) sum(g_em(x, s, cut) * f), numeric(1))
}
PARS$g_center <- g_center_of(AC_BASE$x)              # 默认：正态 AC（linear 0；quad 1；threshold 1 − Φ(1)）

#' 给定 DGM 规格，返回带有该 DGM 中心化常数的参数表（所有生成 / 真值函数都应传入这个 pars）
pars_for <- function(spec, pars = PARS) { pars$g_center <- g_center_of(spec$ac$x, pars$cut); pars }

# ---- 线性预测 ---------------------------------------------------------------
#' @param X data.frame(x1, x2, x3)
#' @param trt "A" / "B" / "C"（可为向量）
#' @param b_em 效应修饰强度（三个协变量共用，按形状校准）
#' @param endpoint "binary" 含截距 alpha；"tte" / "continuous" 不含
lin_pred <- function(X, trt, shape, b_em, endpoint = "binary", pars = PARS, cut = pars$cut) {
  gc <- pars$g_center[[shape]]                                            # EM 项在 AC 人群上中心化
  em  <- b_em * ((g_em(X$x1, shape, cut) - gc) + (g_em(X$x2, shape, cut) - gc) + (g_em(X$x3, shape, cut) - gc))
  eta <- pars$beta * (z_std(X$x1) + z_std(X$x2) + z_std(X$x3)) +
         (trt == "A") * (pars$theta_A + em) + (trt == "B") * (pars$theta_B + em)
  if (endpoint == "binary") eta <- eta + pars$alpha
  eta
}

# ---- Weibull 比例风险：基线函数 ----------------------------------------------
weib_scale0 <- function(pars = PARS) pars$weib_median / log(2)^(1 / pars$weib_shape)
H0 <- function(t, pars = PARS) (t / weib_scale0(pars))^pars$weib_shape          # 累积基线风险
h0 <- function(t, pars = PARS) pars$weib_shape / weib_scale0(pars) * (t / weib_scale0(pars))^(pars$weib_shape - 1)
S_cond <- function(t, eta, pars = PARS) exp(-H0(t, pars) * exp(eta))             # 条件生存函数

# ---- 总体层面的边际效应（给定一个总体样本或一组权重） ------------------------
#' 边际 logOR / 风险差：对 a 臂与 b 臂分别取加权平均风险后作对比
marginal_binary <- function(X, w, trt_a, trt_b, shape, b_em, pars = PARS) {
  w <- w / sum(w)
  pa <- sum(w * plogis(lin_pred(X, trt_a, shape, b_em, "binary", pars)))
  pb <- sum(w * plogis(lin_pred(X, trt_b, shape, b_em, "binary", pars)))
  c(logOR = qlogis(pa) - qlogis(pb), RD = pa - pb)
}

#' 边际均值差（连续结局）
marginal_meandiff <- function(X, w, trt_a, trt_b, shape, b_em, pars = PARS) {
  w <- w / sum(w)
  sum(w * (lin_pred(X, trt_a, shape, b_em, "continuous", pars) -
           lin_pred(X, trt_b, shape, b_em, "continuous", pars)))
}

#' 边际生存函数与边际风险（时间网格 tt 上），用于 Cox 极限与 RMST
#' 分块矩阵运算（每块 chunk 个个体 × 全部时间点），权重为 0 的个体直接跳过；边际风险 h(t) = h0(t)·Σ w e^η S_i(t) / Σ w S_i(t)
marginal_surv <- function(X, w, trt, shape, b_em, tt, pars = PARS, chunk = 5000) {
  keep <- w > 0; w <- w[keep] / sum(w[keep])
  ex  <- exp(lin_pred(X[keep, , drop = FALSE], trt, shape, b_em, "tte", pars))
  H0t <- H0(tt, pars); h0t <- h0(tt, pars); n <- length(ex)
  S <- numeric(length(tt)); Hn <- numeric(length(tt))
  for (i in seq(1, n, by = chunk)) {
    idx <- i:min(i + chunk - 1, n)
    E <- exp(-outer(H0t, ex[idx]))                              # 行 = 时间，列 = 个体：条件生存
    S <- S + as.numeric(E %*% w[idx]); Hn <- Hn + as.numeric(E %*% (ex[idx] * w[idx]))
  }
  list(S = S, h = h0t * Hn / S)
}

#' 条件生存矩阵（时间 × 个体）及其风险加权版本，供同一（不太大的）总体上反复换权重时使用：
#' 边际生存 S(t) = Smat %*% w；边际风险 h(t) = (Hmat %*% w) / S(t)
surv_matrices <- function(X, trt, shape, b_em, tt, pars = PARS) {
  ex   <- exp(lin_pred(X, trt, shape, b_em, "tte", pars))
  Smat <- exp(-outer(H0(tt, pars), ex))                      # 行 = 时间，列 = 个体
  list(S = Smat, H = Smat * outer(h0(tt, pars), ex))
}

#' 由生存矩阵与权重得到边际 S、h
marginal_surv_from <- function(sm, w) {
  w <- w / sum(w); S <- as.numeric(sm$S %*% w); h <- as.numeric(sm$H %*% w) / S
  list(S = S, h = h)
}

#' 两臂边际数据上 Cox 回归的概率极限（1:1 分配，删失生存函数 G）
#' 解总体偏似然得分方程：∫ G(t)[S1 h1 − E(t)(S1 h1 + S0 h0)] dt = 0，E(t) = e^b S1 / (e^b S1 + S0)
cox_limit <- function(S1, h1, S0, h0m, G) {
  score <- function(b) { E <- exp(b) * S1 / (exp(b) * S1 + S0)
    sum(G * (S1 * h1 - E * (S1 * h1 + S0 * h0m))) }
  uniroot(score, c(-4, 4), tol = 1e-10)$root
}

#' 删失生存函数：指数随机删失 + 行政删失
G_cens <- function(tt, rate, admin = PARS$admin_cens) exp(-rate * tt) * (tt < admin)

#' 校准随机删失率，使 AC 试验（A、C 各半）总删失比例 ≈ target
calibrate_cens_rate <- function(X_ac, shape, b_em, target = 0.30, pars = PARS,
                                tt = seq(0.01, pars$admin_cens, by = 0.01)) {
  w <- rep(1, nrow(X_ac))
  SA <- marginal_surv(X_ac, w, "A", shape, b_em, tt, pars)$S
  SC <- marginal_surv(X_ac, w, "C", shape, b_em, tt, pars)$S
  S  <- (SA + SC) / 2
  # 删失比例 = P(C < T) = ∫ S(t)·(−dG(t))：随机删失部分 + 行政删失部分 S(36)
  cens_frac <- function(rate) {
    dG <- rate * exp(-rate * tt) * (tt < pars$admin_cens)
    sum(S * dG) * 0.01 + S[length(S)] * exp(-rate * pars$admin_cens)
  }
  uniroot(function(r) cens_frac(r) - target, c(1e-5, 1), tol = 1e-8)$root
}

#' RMST(τ)：边际生存曲线在 [0, τ] 上的积分
rmst <- function(S, tt, tau) sum(S[tt <= tau]) * (tt[2] - tt[1])

# ---- 迁移差距与 b_EM 校准（第 3.3 节） --------------------------------------
#' 未调整时的迁移差距：A vs C 的边际 logOR 在 AC 人群与 BC 人群之间的差
transport_gap <- function(shape, b_em, X_ac, X_bc, pars = PARS) {
  unname(marginal_binary(X_ac, rep(1, nrow(X_ac)), "A", "C", shape, b_em, pars)["logOR"] -
         marginal_binary(X_bc, rep(1, nrow(X_bc)), "A", "C", shape, b_em, pars)["logOR"])
}

#' 解 b_EM（三个协变量共用）使 |迁移差距| = gap（默认 0.20 logOR）
#' 注意：边际尺度上的差距随 b_EM 先增后减（不可压缩性使异质性越大、边际效应越缩水），
#' 所以目标不能定得太高，并只在单调区间内求根。
calibrate_b_em <- function(shape, X_ac, X_bc, gap = 0.20, pars = PARS, interval = c(0.005, 1.5)) {
  f <- function(b) abs(transport_gap(shape, b, X_ac, X_bc, pars)) - gap
  if (f(interval[1]) * f(interval[2]) > 0)
    stop("形状 ", shape, "：在区间内迁移差距达不到 ", gap, "，请降低 gap")
  uniroot(f, interval, tol = 1e-6)$root
}
