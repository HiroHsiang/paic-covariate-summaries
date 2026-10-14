# 平衡权重（MAIC = 指数倾斜 / 熵平衡），支持精确约束和带容差的软约束（sMAIC）。
# 同一函数既用于样本（基权重均匀），也用于总体极限计算（基权重 = 范围截取后的均匀权重）。
# 由基线表构造约束矩阵的函数在 R/06_maic.R（cov_constraints / build_constraints_multi）。
# 2026-10-10 修订：原先把精确约束也写成 λ = λ⁺ − λ⁻ 交给 L-BFGS-B，λ⁺、λ⁻ 同时增大时目标函数不变（退化方向），
# 偶尔会漂移到 1e34 并报"不可行"（主运行中 R5 的回退率被抬到 16–25%）。现在：精确约束的 λ 自由、用 Newton 法求解；
# 软约束才用 λ⁺ − λ⁻ 的 L1 惩罚形式（惩罚项使该方向不再退化），交给 L-BFGS-B，并以 Newton 法的结果作起点。

#' 指数倾斜权重
#' 对偶问题：min_λ  log Σ_i b_i exp(λ'(h_i − s)) + Σ_k t_k |λ_k|；t_k = 0 即精确约束。
#' @param H n×K 的约束函数矩阵（如 cbind(x, x^2, x <= q)）
#' @param targets 长度 K 的目标值（基线表上的数或其变换）
#' @param tol 长度 K 的容差（0 = 精确）
#' @param base 基权重（默认均匀）；总体极限时传入（范围截取后的）基权重
#' @return list(w = 和为 n 的权重（样本）或和为 1（给了 base）, achieved, feasible, ess, lambda, convergence)
tilt_weights <- function(H, targets, tol = 0, base = NULL, maxit = 5000, feas_tol = 1e-4) {
  H <- as.matrix(H); K <- ncol(H); n <- nrow(H)
  if (length(tol) == 1) tol <- rep(tol, K)
  uniform <- is.null(base) || all(base == base[1])
  if (is.null(base)) base <- rep(1 / n, n) else base <- base / sum(base)
  Z <- sweep(H, 2, targets)                                   # 中心化到目标
  sc <- sqrt(colSums(base * Z^2)); sc[sc < 1e-12] <- 1        # 按列缩放，改善条件数
  Zs <- sweep(Z, 2, sc, "/"); ts <- tol / sc
  ex <- ts <= 0; K_e <- sum(ex); K_s <- K - K_e                # 精确 / 软约束的个数
  fail <- function(code) list(w = if (uniform) rep(1, n) else base, achieved = rep(NA_real_, K), feasible = FALSE,
                              ess = if (uniform) n else 1, lambda = rep(NA_real_, K), convergence = code)
  probs <- function(l) { v <- as.numeric(Zs %*% l); p <- base * exp(v - max(v)); p / sum(p) }
  lse   <- function(l) { v <- as.numeric(Zs %*% l); m <- max(v); m + log(sum(base * exp(v - m))) }

  # ---- 第一步：Newton 法解精确约束部分（软约束的 λ 先置 0） --------------------------------
  l <- rep(0, K); conv <- 0L
  if (K_e > 0) {
    ok <- FALSE
    for (it in 1:200) {
      p <- probs(l); g <- as.numeric(crossprod(Zs, p))[ex]                 # 梯度 = 加权均值与目标之差（缩放后）
      if (max(abs(g)) < 1e-10) { ok <- TRUE; break }
      Ze <- Zs[, ex, drop = FALSE]; M <- crossprod(Ze * sqrt(p)); M <- M - tcrossprod(as.numeric(crossprod(Ze, p)))   # 加权协方差
      step <- tryCatch(solve(M + diag(1e-10, K_e), g), error = function(e) NULL)
      if (is.null(step) || !all(is.finite(step))) break
      f0 <- lse(l); a <- 1; moved <- FALSE
      for (h in 1:30) { l2 <- l; l2[ex] <- l[ex] - a * step; f2 <- lse(l2)   # 步长减半的线搜索
        if (is.finite(f2) && f2 <= f0 - 1e-12 * a * sum(g * step) + 1e-14) { l <- l2; moved <- TRUE; break }; a <- a / 2 }
      if (!moved) break
    }
    if (!ok) { p <- probs(l); g <- as.numeric(crossprod(Zs, p))[ex]; if (!all(is.finite(g)) || max(abs(g)) > feas_tol) conv <- 1L }
  }

  # ---- 第二步：有软约束时，用 L-BFGS-B 解完整的 L1 惩罚对偶（精确 λ 自由，软 λ = λ⁺ − λ⁻ ≥ 0） ----
  if (K_s > 0) {
    ie <- which(ex); is <- which(!ex); tss <- ts[is]
    unpack <- function(u) { lam <- numeric(K); lam[ie] <- u[seq_len(K_e)]
      lam[is] <- u[K_e + seq_len(K_s)] - u[K_e + K_s + seq_len(K_s)]; lam }
    fn <- function(u) lse(unpack(u)) + sum(tss * (u[K_e + seq_len(K_s)] + u[K_e + K_s + seq_len(K_s)]))
    gr <- function(u) { g <- as.numeric(crossprod(Zs, probs(unpack(u)))); c(g[ie], g[is] + tss, -g[is] + tss) }
    u0 <- c(l[ie], rep(0, 2 * K_s))
    op <- tryCatch(optim(u0, fn, gr, method = "L-BFGS-B", lower = c(rep(-Inf, K_e), rep(0, 2 * K_s)),
                         control = list(maxit = maxit, factr = 1e5)), error = function(e) NULL)
    if (is.null(op) || !all(is.finite(op$par))) return(fail(99L))
    l <- unpack(op$par); conv <- as.integer(op$convergence)
  }

  w <- probs(l); lambda <- l / sc
  achieved <- as.numeric(crossprod(H, w))
  gap <- abs(achieved - targets) - tol                          # 超出容差的部分
  feasible <- all(is.finite(gap)) && all(gap <= feas_tol * pmax(1, abs(targets)))   # 相对容差判可行（Newton 的精确部分远小于此）
  if (!all(is.finite(w))) return(fail(99L))
  list(w = if (uniform) w * n else w, achieved = achieved, feasible = feasible,
       ess = ess_frac(w, base) * (if (uniform) n else 1), lambda = lambda, convergence = conv)
}

#' Kish 有效样本量（样本权重）
ess <- function(w) sum(w)^2 / sum(w^2)

#' ESS 占比的极限形式：权重 w 相对基权重 base 的 (Σ b r)² / Σ b r²，r = w/b
ess_frac <- function(w, base) { r <- w / base; sum(base * r)^2 / sum(base * r^2) }
