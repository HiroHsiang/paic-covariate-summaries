# 逐情景的 RMSE 排序（2026-10-12，用户：每个情景是独立的，要回答“什么情景下哪个方法的 RMSE 更低”，不做跨情景汇总）。
# 对每个 DGM（96 个）：在正文方法（表 2 的 11 种）中找 RMSE 最低者；其余方法与最低者在同一批重复上做配对比较
# （逐次重复的平方误差之差，均值与 MCSE），配对 95% MC 区间含 0 的方法视为“与最低者无法区分”。
# 另记录：只用均值 (SD) 的方法中 RMSE 最低者（对应基线表只报均值 (SD) 的情形）。
# 数据与 sims/14 相同（主运行 + 各补丁 + 补充运行），误差 = 估计 − 总体真值。
# 产出：output/summary/rmse_rank_by_scenario.csv（逐 DGM × 方法：RMSE、相对最低者的比值、配对差及其 MCSE、是否可区分）
#       output/summary/best_by_scenario.csv（逐 DGM：最低者、无法区分者、只用均值 (SD) 时的最低者）
# 用法：Rscript sims/35_best_rmse_by_scenario.R（约 1–2 分钟）

truths <- readRDS("output/main/truths.rds")
prim <- function(ep, an) if (ep == "binary") "logOR" else if (an) "logHR_diff" else "RMST_diff"
keep_cols <- c("method", "R", "rep", "est_ab", "D", "shape", "endpoint", "anchored", "em_mult", "dgm_id")
rd <- function(f) { d <- readRDS(f); d[, intersect(keep_cols, names(d))] }
ids <- sub("[.]rds$", "", list.files("output/main", pattern = "^D.*[.]rds$"))
assemble <- function(id) {                                                    # 与 sims/14 的 assemble() 相同
  d <- rd(file.path("output/main", paste0(id, ".rds")))
  pf <- file.path("output/main_patch", paste0(id, ".rds"))
  if (file.exists(pf)) { p <- rd(pf); sel <- d$method %in% c("M3b", "M5") & d$R %in% c("R2", "R2r", "R4", "R5"); d <- rbind(d[!sel, ], p[, names(d)]) }
  pf <- file.path("output/main_patch_stc", paste0(id, ".rds"))
  if (file.exists(pf)) { p <- rd(pf); sel <- d$method == "M2"; d <- rbind(d[!sel, ], p[, names(d)]) }
  for (dir in c("output/supp", "output/extra_m1q", "output/full_kappa1")) { f <- file.path(dir, paste0(id, ".rds")); if (file.exists(f)) { x <- rd(f); d <- rbind(d, x[, names(d)]) } }
  d$dgm_id <- id; d
}

# 正文方法（表 2）；msd = 只用均值 (SD) 的方法
M <- data.frame(key = c("M0 R1", "M1m R1", "M1 R1", "M1 R5", "M1q R5", "M7 R5", "M2 R1", "M3a R1", "M5 R1", "M5 R5", "M5s R5"),
                label = c("Unadjusted (Bucher)", "MAIC (means)", "MAIC (means and SDs)", "MAIC (all reported summaries)", "MAIC (means, quartiles and proportions)",
                          "sMAIC", "STC", "G-computation", "G-MAIC (means and SDs)", "G-MAIC", "sG-MAIC"),
                msd = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)

rank_rows <- list(); best_rows <- list()
for (id in ids) {
  d <- assemble(id); ep <- d$endpoint[1]; an <- d$anchored[1]
  key <- sprintf("%s_%s_em%d_%s", d$D[1], d$shape[1], d$em_mult[1], ep)
  tr <- truths$truth[truths$key == key & truths$estimand == prim(ep, an)]
  d$key <- paste(d$method, d$R); d <- d[d$key %in% M$key, ]; d$err <- d$est_ab - tr
  w <- reshape(d[, c("rep", "key", "err")], idvar = "rep", timevar = "key", direction = "wide")   # 每行一次重复，每列一种方法的误差
  names(w) <- sub("^err[.]", "", names(w))
  ks <- intersect(M$key, names(w)); ks <- ks[sapply(ks, function(k) sum(is.finite(w[[k]])) > 100)]
  rmse <- sapply(ks, function(k) sqrt(mean(w[[k]]^2, na.rm = TRUE)))
  kb <- names(which.min(rmse))
  # 与最低者的配对比较：逐次重复平方误差之差
  cmp <- do.call(rbind, lapply(ks, function(k) {
    v <- w[[k]]^2 - w[[kb]]^2; ok <- is.finite(v)
    data.frame(dgm_id = id, D = d$D[1], shape = d$shape[1], em_mult = d$em_mult[1], endpoint = ep, anchored = an, key = k, label = M$label[M$key == k],
               msd = M$msd[M$key == k], rmse = rmse[[k]], ratio_to_best = rmse[[k]] / rmse[[kb]], dmse = mean(v[ok]), dmse_mcse = sd(v[ok]) / sqrt(sum(ok)),
               n = sum(ok), stringsAsFactors = FALSE)
  }))
  cmp$is_best <- cmp$key == kb
  cmp$indist <- !cmp$is_best & (cmp$dmse - 1.96 * cmp$dmse_mcse <= 0)            # 配对 95% MC 区间含 0（dmse ≥ 0 恒成立，只看下界）
  rank_rows[[id]] <- cmp
  ms <- cmp[cmp$msd, ]; kbm <- ms$key[which.min(ms$rmse)]
  best_rows[[id]] <- data.frame(dgm_id = id, D = d$D[1], shape = d$shape[1], em_mult = d$em_mult[1], endpoint = ep, anchored = an,
    best = M$label[M$key == kb], best_rmse = rmse[[kb]],
    indist = paste(cmp$label[cmp$indist][order(cmp$rmse[cmp$indist])], collapse = "; "),
    indist_ratio = paste(sprintf("%.3f", cmp$ratio_to_best[cmp$indist][order(cmp$rmse[cmp$indist])]), collapse = "; "),
    second = cmp$label[order(cmp$rmse)][2], second_ratio = sort(cmp$ratio_to_best)[2],
    best_msd = M$label[M$key == kbm], best_msd_ratio = rmse[[kbm]] / rmse[[kb]],
    ratio_means = rmse[["M1m R1"]] / rmse[[kb]], stringsAsFactors = FALSE)
}
rank <- do.call(rbind, rank_rows); best <- do.call(rbind, best_rows)
write.csv(rank, "output/summary/rmse_rank_by_scenario.csv", row.names = FALSE)
write.csv(best, "output/summary/best_by_scenario.csv", row.names = FALSE)
cat("written rmse_rank_by_scenario.csv and best_by_scenario.csv:", nrow(best), "DGMs\n")
