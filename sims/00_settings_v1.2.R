# 00_settings_v1.2.R
# ------------------------------------------------------------
# CN: 本模拟研究（protocol v1.2，2026-10-10；决策 18、19）的全局设置与场景网格 pc；只生成设置对象，供审阅与下游脚本读取
# EN: Global settings and scenario grid for the simulation study (protocol v1.2)
# 风格参照 ITC-sim2/unanchored_settings.R。可打包代码：不写绝对路径；校准值从 output/poplimit/dgm_params.rds 读取。
# 用法：Rscript sims/00_settings_v1.2.R   → 打印设置摘要，保存 output/settings/settings_v1.2.rds 与 pc_v1.2.csv
# ------------------------------------------------------------
for (f in list.files("R", pattern = "^0[1-2]_", full.names = TRUE)) source(f)
pr <- readRDS("output/poplimit/dgm_params.rds")          # sims/01 校准产物：D1–D4 两试验规格、b_EM、删失率、中心化常数
dir.create("output/settings", recursive = TRUE, showWarnings = FALSE)

## 1) 全局 | Global
replicates  <- 1000                                       # 主运行每个 DGM 的重复次数
replicates_mlnmr <- 500                                   # ML-NMR 子集
replicates_cov   <- 300                                   # G-computation / G-MAIC 区间（bootstrap B = 100）
endpoints   <- c("binary", "tte")                         # 主研究；continuous 仅附录
anchoring   <- c(anchored = TRUE, unanchored = FALSE)
N_AC <- 300; N_BC <- 300; allocation <- 0.5               # 两试验样本量与分配比
tau_rmst    <- 24                                         # RMST 的 τ（月）
seed_main   <- 20261100                                   # 主运行种子 = seed_main + DGM 行号（doRNG 并行流）
ncores      <- 12

## 2) 协变量 | Covariates（p = 3，全部连续、抽象"评分型"尺度；决策 18：三个都是预后因子 + 效应修饰因子）
p        <- 3
EM_idx   <- 1:3
covariates <- data.frame(
  var   = c("x1", "x2", "x3"),
  type  = "continuous (score-type, no clinical label)",
  AC    = "mean 50, SD 10; shape by D",
  BC    = "mean 53 (+0.3 SD), SD 12 (x1.2); same shape as AC (D4: truncated at 70)",
  prognostic = sprintf("beta = %.4f per SD", PARS$beta),
  EM_shape = "same g for all three (linear / U-shape / threshold at 65 = AC mean + 1.5 SD)",
  reported = "per availability R (all three reported the same way)")
rho      <- 0.35                                          # 高斯 copula 相关系数（两试验相同）

## 3) 分布情景 D | Distribution scenarios（决策 19：两试验共享形状；BC = AC 的位置–尺度变换，D4 另加截断）
D_scenarios <- pr$dist_tab[, c("D", "trial", "type", "mean", "sd", "skew", "Q1", "median", "Q3", "p_gt_cut", "mix_p", "mix_m1", "mix_m2", "mix_s")]
D_overlap   <- pr$ess_tab                                 # oracle 权重的 ESS/N（一维 / 三维）
D_meaning   <- c(D1 = "both normal (control: mean + SD exact)", D2 = "both right-skewed, skew +0.8 (90.5% + 9.5% shoulder)",
                 D3 = "both left-skewed, skew -0.8 (mirror)", D4 = "both normal; BC truncated at 70 by eligibility (post-truncation mean 53)")

## 4) 效应修饰 | Effect modification（三个协变量同形状、同强度；EM 项在各 D 的 AC 人群上中心化）
EM_shapes <- c(linear = "z", quad = "z^2", threshold = "1{X > 65}")   # z = (X − 50)/10；切点 = AC 均值 + 1.5 SD（决策 21）
EM_center <- pr$g_center                                  # E_AC[g]：按 D（AC 形状）数值积分
b_EM      <- pr$b_em                                      # 中等档：D1 下迁移差距 0.20 logOR
em_mult   <- c(moderate = 1, strong = 2)                  # 强档 = 系数翻倍
gap_table <- pr$gap_tab                                   # 各 D 下的实际迁移差距（logOR），中等与强档

## 5) 结局生成参数 | Endpoint generator parameters（η = α + β·(z1 + z2 + z3) + T_A(θ_A + EM) + T_B(θ_B + EM)，EM = b_EM Σ_j [g(z_j) − c_D]）
bin_params  <- list(alpha = PARS$alpha, beta = PARS$beta, theta_A = PARS$theta_A, theta_B = PARS$theta_B,
                    link = "logit", baseline_risk_C = plogis(PARS$alpha),
                    published_BC = "logOR (glm) with SE; RD with SE (secondary)",
                    published_B_unanchored = "event proportion with binomial SE",
                    estimand = "marginal logOR A vs B in BC superpopulation; RD secondary")
surv_params <- list(dist = "Weibull PH", shape = PARS$weib_shape, baseline_median = PARS$weib_median,
                    scale0 = weib_scale0(PARS), beta = PARS$beta, theta_A = PARS$theta_A, theta_B = PARS$theta_B,
                    cens_rate = pr$cens_rate, cens_target_AC = 0.30, admin_cens = PARS$admin_cens, tau = tau_rmst,
                    published_BC = "logHR (Cox) with SE; RMST(24) difference with SE (secondary)",
                    published_B_unanchored = "KM-based RMST(24) with SE",
                    estimand_anchored = "difference of the two vs-C marginal Cox limits under the common censoring law; RMST secondary",
                    estimand_unanchored = "RMST(24) difference")
cont_params <- list(sigma = 1, beta = PARS$beta, theta_A = PARS$theta_A, theta_B = PARS$theta_B,
                    note = "appendix only; unit test (D1 x quad x mean+SD must be exact)")

## 6) 报告方式（可得性） | Availability scenarios（三个协变量同一方式报告）
availability <- data.frame(
  R = c("R1", "R2", "R2r", "R3", "R3m", "R4", "R5"),
  reported = c("mean (SD)", "median [IQR]", "median (range)", "mean (SD) + n(%) > 65", "mean (SD) + n(%) > 60 (misaligned by 0.5 SD)",
               "mean (SD) + median [IQR]", "mean (SD) + median [IQR] + range + n(%) > 65"),
  constraints_R5_per_cov = c(2, 3, 1, 3, 3, 5, 6), rounding = "1 dp (proportions 0.1%)")

## 7) 方法 | Methods
methods <- data.frame(
  code = c("M0", "M1m", "M1c", "M1", "M7", "M2", "M3a", "M3b", "M5", "M5s", "M9"),
  method = c("unadjusted (Bucher / naive)", "MAIC means only (median converted)", "MAIC converted mean+SD", "MAIC all exact (+fallback)",
             "sMAIC (means/proportions exact; SD, quantiles soft, kappa 1.96)", "plug-in STC (conditional)",
             "G-computation, normal marginals", "G-computation, reconstructed marginals", "G-MAIC (exact weights)",
             "soft G-MAIC (sMAIC weights)", "ML-NMR (multinma, anchored binary only)"),
  run = c("supp", "main", "main", "main", "main", "main", "main", "main", "main", "supp", "mlnmr"))

## 8) 场景网格 pc | Scenario grid（与 sims/10 的顺序一致：二分类在前；种子 = seed_main + k）
pc <- expand.grid(endpoint = endpoints, anchored = c(TRUE, FALSE), em_mult = em_mult,
                  shape = names(EM_shapes), D = names(pr$specs), stringsAsFactors = FALSE)
pc <- pc[order(pc$endpoint != "binary", pc$D, pc$shape, pc$em_mult, !pc$anchored), ]
pc$k <- seq_len(nrow(pc)); pc$seed <- seed_main + pc$k
pc$N_AC <- N_AC; pc$N_BC <- N_BC; pc$rho <- rho
pc$b_EM <- b_EM[pc$shape] * pc$em_mult
pc$cens_rate <- ifelse(pc$endpoint == "tte", pr$cens_rate, NA)
pc$file_id <- sprintf("%s_%s_em%d_%s_%s", pc$D, pc$shape, pc$em_mult, pc$endpoint, ifelse(pc$anchored, "anch", "unanch"))
rownames(pc) <- NULL

## 9) 汇总设置对象 s | Settings object
s <- list(version = "v1.2", date = "2026-10-10", replicates = replicates, replicates_mlnmr = replicates_mlnmr,
          replicates_cov = replicates_cov, endpoints = endpoints, N_AC = N_AC, N_BC = N_BC, allocation = allocation,
          tau_rmst = tau_rmst, seed_main = seed_main, ncores = ncores,
          p = p, EM_idx = EM_idx, covariates = covariates, rho = rho, D_scenarios = D_scenarios, D_overlap = D_overlap, D_meaning = D_meaning,
          EM_shapes = EM_shapes, EM_center = EM_center, b_EM = b_EM, em_mult = em_mult, gap_table = gap_table,
          bin_params = bin_params, surv_params = surv_params, cont_params = cont_params,
          availability = availability, methods = methods, pc = pc)
saveRDS(s, "output/settings/settings_v1.2.rds"); write.csv(pc, "output/settings/pc_v1.2.csv", row.names = FALSE)

## 10) 打印摘要 | Print
options(width = 170)
cat("== settings", s$version, "==\n")
cat("replicates:", replicates, "| N_AC/N_BC:", N_AC, "/", N_BC, "| rho:", rho, "| endpoints:", paste(endpoints, collapse = ", "), "| DGMs:", nrow(pc), "\n\n")
cat("-- covariates --\n"); print(t(covariates[1, -1]))
cat("\n-- distribution scenarios D (both trials; three covariates identical) --\n"); print(D_scenarios, row.names = FALSE, digits = 3)
cat("\n-- meaning --\n"); for (D in names(D_meaning)) cat(" ", D, ":", D_meaning[[D]], "\n")
cat("\n-- overlap: oracle-weight ESS/N --\n"); print(D_overlap, row.names = FALSE)
cat("\n-- effect modification --\n"); print(data.frame(shape = names(EM_shapes), g = EM_shapes, b_EM_moderate = round(b_EM, 4), b_EM_strong = round(2 * b_EM, 4), row.names = NULL))
cat("centring constants E_AC[g] by D:\n"); print(EM_center, row.names = FALSE)
cat("realised transport gap (logOR; moderate and strong):\n"); print(gap_table, row.names = FALSE)
cat("\n-- binary --\n"); str(bin_params, give.attr = FALSE)
cat("\n-- time-to-event --\n"); str(surv_params, give.attr = FALSE)
cat("\n-- continuous (appendix) --\n"); str(cont_params, give.attr = FALSE)
cat("\n-- availability --\n"); print(availability, row.names = FALSE)
cat("\n-- methods --\n"); print(methods, row.names = FALSE)
cat("\n-- pc (first 8 of", nrow(pc), ") --\n"); print(head(pc, 8), row.names = FALSE)
