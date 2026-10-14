# ESM 表 S1（2026-10-13，决策 32）：逐情景、逐方法的偏倚与 RMSE。原正文表 3（逐情景 RMSE 最低方法）按用户意见撤出正文，
# 并入本表：每格 = 偏倚 / RMSE；RMSE 最低者加粗，与之无法区分者（同一批重复的配对 95% MC 区间含 0）用斜体；
# 末列 = 若 RMSE 最低者用到了均值和标准差以外的汇总量，最多只用均值和标准差的方法中最低 RMSE 高出的百分比。
# 每种设计 × 结局一个分表（a–d），行 = 24 个情景（顺序与图 4 相同），列 = 表 2 的 11 种方法。
# 输入：output/summary/rmse_rank_by_scenario.csv、best_by_scenario.csv（sims/35）、perf_long.csv（sims/14）。
# 产出：output/figures/tableS1_bias_rmse.md（英文）、tableS1_bias_rmse_zh.md（中文表头，数字相同）、.csv
# 用法：Rscript sims/36_tableS1_bias_rmse.R

rk <- read.csv("output/summary/rmse_rank_by_scenario.csv"); b <- read.csv("output/summary/best_by_scenario.csv")
p <- read.csv("output/summary/perf_long.csv"); p$key <- paste(p$method, p$R)
rk <- merge(rk, p[, c("dgm_id", "key", "bias")], by = c("dgm_id", "key"))
meth <- c("Unadjusted (Bucher)", "MAIC (means)", "MAIC (means and SDs)", "MAIC (all reported summaries)", "MAIC (means, quartiles and proportions)",
          "sMAIC", "STC", "G-computation", "G-MAIC (means and SDs)", "G-MAIC", "sG-MAIC")
D_lab <- c(D1 = "Normal", D2 = "Right-skewed", D3 = "Left-skewed", D4 = "Truncated comparator")
D_zh <- c(D1 = "正态", D2 = "右偏", D3 = "左偏", D4 = "对照试验截断")
S_lab <- c(linear = "Linear", quad = "U-shaped", threshold = "Threshold"); S_zh <- c(linear = "线性", quad = "U 形", threshold = "阈值型")
designs <- data.frame(ep = c("binary", "binary", "tte", "tte"), an = c(TRUE, FALSE, TRUE, FALSE), part = c("a", "b", "c", "d"),
  en = c("Anchored comparison, binary outcome (marginal log odds ratio)", "Unanchored comparison, binary outcome (marginal log odds ratio)",
         "Anchored comparison, time-to-event outcome (marginal log hazard ratio)", "Unanchored comparison, time-to-event outcome (difference in 24-month RMST, months)"),
  zh = c("锚定比较，二分类结局（边际对数比值比）", "无锚定比较，二分类结局（边际对数比值比）", "锚定比较，生存结局（边际对数风险比）", "无锚定比较，生存结局（24 个月 RMST 差，月）"))
fm <- function(x) { x[abs(x) < 0.005] <- 0; gsub("-", "−", sprintf("%.2f", x)) }              # 负号用 U+2212；舍入为 0 的不带负号
grid <- expand.grid(shape = names(S_lab), em_mult = c(1, 2), D = names(D_lab), stringsAsFactors = FALSE)
grid <- grid[order(match(grid$D, names(D_lab)), match(grid$shape, names(S_lab)), grid$em_mult), ]

build <- function(lang) {
  md <- character(0); csv <- list()
  for (j in seq_len(nrow(designs))) {
    hdr_m <- meth; if (lang == "zh") hdr_m[1] <- "Unadjusted (Bucher)"
    md <- c(md, if (lang == "en") sprintf("**Table S1%s** %s: bias / RMSE of each method in each scenario", designs$part[j], designs$en[j])
                else sprintf("**表 S1%s** %s：各方法在每个情景中的偏倚 / RMSE", designs$part[j], designs$zh[j]), "",
            paste0("| ", if (lang == "en") "Covariate distribution | Effect modification | " else "协变量分布 | 效应修饰 | ", paste(hdr_m, collapse = " | "),
                   if (lang == "en") " | Lowest RMSE with means and SDs only |" else " | 只用均值和标准差时的最低 RMSE |"),
            paste0("|", paste(rep("---|", length(meth) + 3), collapse = "")))
    for (i in seq_len(nrow(grid))) {
      z <- rk[rk$D == grid$D[i] & rk$shape == grid$shape[i] & rk$em_mult == grid$em_mult[i] & rk$endpoint == designs$ep[j] & rk$anchored == designs$an[j], ]
      bb <- b[b$D == grid$D[i] & b$shape == grid$shape[i] & b$em_mult == grid$em_mult[i] & b$endpoint == designs$ep[j] & b$anchored == designs$an[j], ]
      cells <- sapply(meth, function(m) { y <- z[z$label == m, ]; if (!nrow(y)) return("—")
        r <- fm(y$rmse); if (y$is_best) r <- paste0("**", r, "**") else if (y$indist) r <- paste0("*", r, "*")
        paste0(fm(y$bias), " / ", r) })
      pen <- if (bb$best_msd_ratio <= 1 + 1e-12) "0%" else if (100 * (bb$best_msd_ratio - 1) < 1) "<1%" else sprintf("+%.0f%%", 100 * (bb$best_msd_ratio - 1))
      dist <- if (i %% 6 == 1) (if (lang == "en") D_lab[grid$D[i]] else D_zh[grid$D[i]]) else ""
      sc <- if (lang == "en") sprintf("%s, %s", S_lab[grid$shape[i]], ifelse(grid$em_mult[i] == 1, "moderate", "strong")) else sprintf("%s，%s", S_zh[grid$shape[i]], ifelse(grid$em_mult[i] == 1, "中等", "强"))
      md <- c(md, sprintf("| %s | %s | %s | %s |", dist, sc, paste(cells, collapse = " | "), pen))
      if (lang == "en") csv[[length(csv) + 1]] <- data.frame(design = designs$en[j], distribution = D_lab[grid$D[i]], scenario = sc, method = meth,
        bias = sapply(meth, function(m) { y <- z[z$label == m, ]; if (nrow(y)) y$bias else NA }), rmse = sapply(meth, function(m) { y <- z[z$label == m, ]; if (nrow(y)) y$rmse else NA }),
        lowest = sapply(meth, function(m) { y <- z[z$label == m, ]; if (nrow(y)) y$is_best else NA }), indistinguishable = sapply(meth, function(m) { y <- z[z$label == m, ]; if (nrow(y)) y$indist else NA }),
        msd_only_penalty = bb$best_msd_ratio - 1, row.names = NULL)
    }
    md <- c(md, "")
  }
  note_en <- paste0("Each cell gives bias / RMSE (300 patients per trial; 1000 replicates per scenario). Bold: lowest RMSE among the methods in the scenario; italic: RMSE not distinguishable from the lowest ",
    "(paired 95% Monte Carlo interval of the difference in mean squared error, computed on the same replicates, included zero). The last column gives the amount by which the lowest RMSE among methods that use at most the means and SDs ",
    "(Unadjusted (Bucher), MAIC (means), MAIC (means and SDs), STC, G-computation, G-MAIC (means and SDs)) exceeded the lowest RMSE; 0% when the method with the lowest RMSE used at most the means and SDs. ",
    "Monte Carlo standard errors: bias at most 0.017 on the log scales and 0.032 months. MAIC, matching-adjusted indirect comparison; RMSE, root mean squared error; RMST, restricted mean survival time; SD, standard deviation; STC, simulated treatment comparison.")
  note_zh <- paste0("每格为偏倚 / RMSE（每个试验 300 例；每个情景 1000 次重复）。加粗：该情景中 RMSE 最低者；斜体：与最低者无法区分（在同一批重复上计算的均方误差之差，其配对 95% 蒙特卡洛区间包含 0）。",
    "末列：最多只用均值和标准差的方法（Unadjusted (Bucher)、MAIC (means)、MAIC (means and SDs)、STC、G-computation、G-MAIC (means and SDs)）中最低的 RMSE 比全部方法中最低的 RMSE 高出的百分比；RMSE 最低者本身只用均值和标准差时为 0%。",
    "蒙特卡洛标准误：偏倚在对数尺度上不超过 0.017，RMST 差不超过 0.032 个月。")
  list(md = c(md, if (lang == "en") note_en else note_zh), csv = if (lang == "en") do.call(rbind, csv) else NULL)
}
en <- build("en"); zh <- build("zh")
writeLines(en$md, "output/figures/tableS1_bias_rmse.md"); writeLines(zh$md, "output/figures/tableS1_bias_rmse_zh.md")
write.csv(en$csv, "output/figures/tableS1_bias_rmse.csv", row.names = FALSE)
cat(en$md[1:8], sep = "\n")
