# 示例：用一张发表的基线表做人群调整间接比较，比较只用均值、均值 + SD 与使用全部汇总量的方法
# ------------------------------------------------------------------------------------------
# 用法（仓库根目录）：Rscript examples/apply_to_published_table.R      （约 1–2 分钟；首次运行要缓存偏正态分位数表）
#
# 数据：用本研究的数据生成器造一次“申报”——IPD 试验（A vs C，各 150 例）与对照试验（B vs C）的发表物。
# 情景：协变量右偏、阈值型效应修饰（65 岁以上获益改变）、强效应修饰、二分类结局、锚定比较。
# 换成自己的数据时：
#   1. ipd_ac：data.frame，列 trt（"A"/"C"）、x1、x2、x3（三个连续协变量）、y（0/1 结局）；
#   2. pub_bc：发表的 B vs C 对数比值比及其标准误；
#   3. vis：按下面“第 2 步”的格式手工填写发表的基线表。
# 目前的代码固定为三个连续协变量、列名 x1–x3（与论文的模拟设计一致）。

for (f in list.files("R", pattern = "^0[0-8]_", full.names = TRUE)) source(f)
suppressPackageStartupMessages({ library(survival); library(sn); library(estmeansd); library(numDeriv) })
set.seed(2026)

# ---- 第 1 步：造一次申报（真实应用中这一步换成自己的 IPD 与发表物） --------------------------
pr <- readRDS("output/poplimit/dgm_params.rds")                         # 研究中校准好的设计参数
D <- "D2"; shape <- "threshold"; em_mult <- 2; endpoint <- "binary"
spec <- pr$specs[[D]]; b_em <- pr$b_em[[shape]] * em_mult; pars <- pars_for(spec)
rep <- gen_replicate(spec, shape, b_em, endpoint, pr$cens_rate, anchored = TRUE, N_AC = 300, N_BC = 300, pars = pars)
ipd_ac <- rep$ipd_ac                                                      # 分析者手里的 IPD
cat("IPD 试验：", nrow(ipd_ac), "例；前几行：\n"); print(head(ipd_ac, 3))
cat("\n发表的 B vs C 对数比值比：", round(rep$pub_bc$effect[["est"]], 3), "（SE", round(rep$pub_bc$effect[["se"]], 3), "）\n")

# ---- 第 2 步：发表的基线表（舍入到一位小数） -----------------------------------------------
# vis_all：表上报告的全部汇总量；vis_msd：只看均值 (SD)。真实应用中按下面的格式手工填写，例如：
#   vis_all <- list(
#     x1 = list(n = 300, mean = 53.0, sd = 12.0, q = c(45.0, 53.1, 61.3), range = c(21.4, 88.9), p_cut = 0.157, cut = 65),
#     x2 = list(...), x3 = list(...), n = 300, R = "R5")
#   其中 q = 第 1 四分位数、中位数、第 3 四分位数；p_cut = 高于切点 cut 的患者比例（0–1）；没有报告的项删掉即可。
vis_all <- visible_table(rep$pub_bc$table, "R5")
vis_msd <- visible_table(rep$pub_bc$table, "R1")
cat("\n基线表中 x1 的汇总量：\n"); str(vis_all$x1)

# ---- 第 3 步：各方法的估计（A vs B 的边际对数比值比） ----------------------------------------
# 加权方法的区间用把权重视为固定的稳健标准误；G 类方法用 IPD 的 bootstrap（这里 B = 100）
B <- 100
res <- rbind(
  cbind(label = "MAIC (means)",                    run_maic_method(rep, vis_msd, "M1m", endpoint, do_boot = FALSE)),
  cbind(label = "MAIC (means and SDs)",            run_maic_method(rep, vis_msd, "M1",  endpoint, do_boot = FALSE)),
  cbind(label = "MAIC (all reported summaries)",   run_maic_method(rep, vis_all, "M1",  endpoint, do_boot = FALSE)),
  cbind(label = "sMAIC",                           run_maic_method(rep, vis_all, "M7",  endpoint, do_boot = FALSE)),
  cbind(label = "STC",                             run_stc_method(rep, vis_msd, endpoint)),
  cbind(label = "G-computation",                   run_gcomp_method(rep, vis_msd, "M3a", endpoint, B = B)),
  cbind(label = "G-MAIC (means and SDs)",          run_gcomp_method(rep, vis_msd, "M5",  endpoint, B = B)),
  cbind(label = "G-MAIC",                          run_gcomp_method(rep, vis_all, "M5",  endpoint, B = B)),
  cbind(label = "sG-MAIC",                         run_gcomp_method(rep, vis_all, "M5s", endpoint, B = B)))
res$se <- ifelse(is.na(res$se_boot_joint), res$se_fixedw, res$se_boot_joint)

# 该情景的真值（对照试验人群中的边际对数比值比；大样本模拟，约 10 秒）
truth <- dgm_truth(spec, shape, b_em, endpoint, pr$cens_rate, n = 2e5)$truth[["logOR"]]

out <- data.frame(method = res$label, estimate = round(res$est_ab, 3), se = round(res$se, 3),
                  lower = round(res$est_ab - 1.96 * res$se, 3), upper = round(res$est_ab + 1.96 * res$se, 3),
                  ESS = round(res$ess), fallback = res$fallback)
cat("\nA vs B 的边际对数比值比（真值", round(truth, 3), "）：\n"); print(out, row.names = FALSE)
cat("\nESS：加权方法的有效样本量（IPD 共 300 例）；fallback：约束不可行时依次去掉范围、比例、四分位数、SD 的层级（0 = 未回退）。\n",
    "STC 估计的是条件效应，与边际真值不可直接比较。单次重复的估计与真值之差主要来自抽样误差（包括各方法共用的发表 B vs C 估计）；\n",
    "方法的表现须看 1000 次重复的结果，见论文与 output/summary/。\n", sep = "")
