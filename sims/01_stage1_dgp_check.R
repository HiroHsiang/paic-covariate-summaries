# 阶段 1 · 第一步：定两个试验的"病人画像"，并校准效应修饰强度与删失率（protocol v1.2 第 3.2–3.4 节）
# 决策 18：三个连续协变量、同形状同强度的 EM；决策 19：两个试验共享分布族（D1 正态、D2 右偏、D3 左偏、D4 BC 截断）。
# 全部是确定性计算（数值积分 + 固定种子的大样本），运行约 2 分钟。
# 用法：在项目根目录运行  Rscript sims/01_stage1_dgp_check.R
# 产出：output/poplimit/dgm_params.rds（后续脚本读取）、output/poplimit/dgm_check.md（人读）

source("R/01_dgm_spec.R"); source("R/02_outcome_model.R"); source("R/03_balancing.R")
dir.create("output/poplimit", recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()
Ds <- c("D1", "D2", "D3", "D4")

# ---- 1. 四个 DGM 的两试验边际分布：矩、分位数、切点以上比例、oracle 权重的 ESS ------------------
specs <- lapply(Ds, dgm_spec); names(specs) <- Ds
x <- seq(-20, 140, by = 0.01)
dist_tab <- do.call(rbind, lapply(Ds, function(D) {
  sp <- specs[[D]]
  rows <- lapply(c("AC", "BC"), function(tr) {
    sx <- if (tr == "AC") sp$ac$x else sp$bc$x; mom <- marg_moments(sx, PARS$cut, x)
    data.frame(D = D, trial = tr, type = sx$type, mean = mom$mean, sd = mom$sd, skew = mom$skew,
               Q1 = mom$q[1], median = mom$q[2], Q3 = mom$q[3], p_gt_cut = mom$p_cut,
               comp_ratio = if (sx$type == "mixture") sx$comp_ratio else NA, mix_p = if (sx$type == "mixture") sx$p else NA,
               mix_m1 = if (sx$type == "mixture") sx$m1 else NA, mix_m2 = if (sx$type == "mixture") sx$m2 else NA,
               mix_s = if (sx$type == "mixture") sx$s else NA)
  })
  do.call(rbind, rows)
}))
# 一维与三维 oracle ESS（AC → BC 的真实密度比权重）；三维用 2e5 的 AC 样本
ess_tab <- do.call(rbind, lapply(Ds, function(D) {
  sp <- specs[[D]]
  f_ac <- marg_dist(sp$ac$x)$d(x); f_ac <- f_ac / sum(f_ac); f_bc <- marg_dist(sp$bc$x)$d(x); f_bc <- f_bc / sum(f_bc)
  ok <- f_ac > 1e-12
  o3 <- oracle_ess_3d(sp)
  data.frame(D = D, ess_oracle_1d = 1 / sum(f_bc[ok]^2 / f_ac[ok]), ess_oracle_3d = o3[["ess_oracle_3d"]],
             w_q999_3d = o3[["max_w_central"]])
}))
dist_tab[, -(1:3)] <- round(dist_tab[, -(1:3)], 3); ess_tab[, -1] <- round(ess_tab[, -1], 3)
cat("== 表 1：两试验的协变量边际分布（三个协变量相同；AC 均值 50 / SD 10，BC 均值 53 / SD 12） ==\n")
print(dist_tab, row.names = FALSE)
cat("\n== 表 1b：oracle 权重（真实密度比）的 ESS/N：一维（单个协变量）与三维（三个协变量，rho = 0.35） ==\n")
print(ess_tab, row.names = FALSE)
cat("\n合格判定（三维 oracle ESS ≥ 0.35）：", paste(ess_tab$D, ifelse(ess_tab$ess_oracle_3d >= 0.35, "OK", "偏低"), collapse = "; "), "\n")

# ---- 2. 校准 b_EM（三个协变量共用）：D1 下三种形状的迁移差距都等于 0.20 logOR ----------------
n_pop <- 4e5                                              # 大样本近似总体（MC 误差约 1/√n，可忽略）
pops <- lapply(Ds, function(D) list(ac = population_draw(specs[[D]]$ac, n_pop, seed = 20261008),
                                    bc = population_draw(specs[[D]]$bc, n_pop, seed = 20261009)))
names(pops) <- Ds
shapes <- c("linear", "quad", "threshold")
b_em <- vapply(shapes, function(s) calibrate_b_em(s, pops$D1$ac, pops$D1$bc, gap = 0.20, pars = pars_for(specs$D1)), numeric(1))
gap_tab <- do.call(rbind, lapply(shapes, function(s) {
  g <- vapply(Ds, function(D) transport_gap(s, b_em[[s]], pops[[D]]$ac, pops[[D]]$bc, pars_for(specs[[D]])), numeric(1))
  g2 <- vapply(Ds, function(D) transport_gap(s, 2 * b_em[[s]], pops[[D]]$ac, pops[[D]]$bc, pars_for(specs[[D]])), numeric(1))
  data.frame(shape = s, b_em = round(b_em[[s]], 4), t(round(g, 3)), strong_b_em = round(2 * b_em[[s]], 4), t(round(g2, 3)))
}))
names(gap_tab) <- c("shape", "b_em", paste0("gap_", Ds), "b_em_x2", paste0("gap2_", Ds))
cat("\n== 表 2：b_EM（中等档）与各 D 下未调整的迁移差距（A vs C 的边际 logOR，AC 人群 − BC 人群）；gap2 = 强档（×2） ==\n")
print(gap_tab, row.names = FALSE)

# 中心化常数与条件效应的量级（供理解）
gc_tab <- do.call(rbind, lapply(Ds, function(D) data.frame(D = D, t(round(pars_for(specs[[D]])$g_center, 4)))))
cat("\n== 表 3：g 在各 D 的 AC 人群上的均值（EM 项的中心化常数） ==\n"); print(gc_tab, row.names = FALSE)
cat("\n条件效应修饰量级（logOR，每个协变量）：线性型每 1 SD =", round(b_em[["linear"]], 3),
    "；U 形 +1 SD (z=1) 相对均值 =", round(b_em[["quad"]], 3),
    "；阈值型 > 65 的跳变 =", round(b_em[["threshold"]], 3),
    "；A vs C 的基准条件效应 =", round(PARS$theta_A, 3), "\n")

# ---- 3. 校准随机删失率（D1 的 AC 试验、线性 EM 下总删失 ≈ 30%，所有 DGM 共用） -------------
cens_rate <- calibrate_cens_rate(pops$D1$ac[sample(n_pop, 5e4), ], "linear", b_em[["linear"]], target = 0.30, pars = pars_for(specs$D1))
cat("\n随机删失率 =", signif(cens_rate, 4), "（另有 36 个月行政删失）\n")

# ---- 4. 保存 --------------------------------------------------------------------
for (D in Ds) specs[[D]]$g_center <- pars_for(specs[[D]])$g_center       # 记录用；生成函数仍通过 pars_for() 取
params <- list(specs = specs, b_em = b_em, cens_rate = cens_rate, dist_tab = dist_tab, ess_tab = ess_tab, gap_tab = gap_tab,
               g_center = gc_tab, n_pop = n_pop, seeds = c(ac = 20261008, bc = 20261009), version = "v1.2", date = Sys.Date())
saveRDS(params, "output/poplimit/dgm_params.rds")
md_table <- function(df) c(paste0("| ", paste(names(df), collapse = " | "), " |"),   # data.frame → markdown 表
                           paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
                           apply(df, 1, function(r) paste0("| ", paste(trimws(r), collapse = " | "), " |")))
md <- c("# DGM check (Stage 1, step 1; protocol v1.2)", paste0("Generated ", Sys.time()), "",
        "## Table 1. Marginal distributions of X1-X3 in the two trials (AC mean 50 / SD 10; BC mean 53 / SD 12; D4 truncated at 70)", "",
        md_table(dist_tab), "", "## Table 1b. Oracle-weight ESS/N (1-D and 3-D, rho = 0.35)", "", md_table(ess_tab), "",
        "## Table 2. b_EM (shared by X1-X3) and unadjusted transport gap (marginal logOR, AC pop - BC pop); gap2 = strong EM (x2)", "",
        md_table(gap_tab), "", "## Table 3. Centring constants E_AC[g(X)]", "", md_table(gc_tab), "",
        paste0("Censoring rate (exponential, shared by both trials): ", signif(cens_rate, 4)))
writeLines(md, "output/poplimit/dgm_check.md")
cat("\n已保存 output/poplimit/dgm_params.rds 与 dgm_check.md；用时", round(difftime(Sys.time(), t0, units = "secs")), "秒\n")
