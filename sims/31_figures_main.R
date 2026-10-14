# 正文三张主图（决策 29 之后的最终方案，2026-10-12）：
#   图 2 / 图 3：锚定 / 无锚定的偏倚点图网格（Gao 2025、Phillippo 2020 的样式）：分面 = 分布 D × EM 形状，纵轴方法，横轴偏倚。
#             同一行上画两个点：空心圈 = 中等档 EM，实心点 = 强档 EM（同一条水平线，便于比较）；点后的短粗线 = 95% MC 区间；
#             自 0 出发的细线 = 偏倚的“茎”；超出横轴范围的点（只会是 STC）截到边界并在旁边标出数值（中等 / 强）。a = 二分类，b = 生存。
#   图 4：配对的 RMSE 比瓦片图：各方法相对“只均值 MAIC”在同一批重复内的 RMSE 比，96 个 DGM × 7 种方法；配对 95% MC 区间含 1 的格子灰色；覆盖率 < 0.925 加 ▲。
# 配色：Nature（NPG）色板（ggsci::pal_npg 的十六进制值），并用点形区分；字体 Arial。输出 output/figures/fig{2,3,4}*.png（图 1 示意图见 sims/33）（300 dpi）与 .pdf。
# 数据：output/summary/perf_long.csv、perf_paired.csv（sims/14）。用法：Rscript sims/31_figures_main.R

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
perf <- read.csv("output/summary/perf_long.csv"); paired <- read.csv("output/summary/perf_paired.csv")
dir.create("output/figures", showWarnings = FALSE)
npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# ---- 方法定义（正文 9 种；各自“规范”的报告方式）。2026-10-13 用户确认加入 G-MAIC (means and SDs)（M5 @R1），
#      配色取 NPG 浅青（与 G-MAIC 的深青同族） -------------------------------------------------------------------
M <- data.frame(method = c("M1m", "M1", "M1", "M7", "M2", "M3a", "M5", "M5", "M5s"),
                R = c("R1", "R1", "R5", "R5", "R1", "R1", "R1", "R5", "R5"),
                label = c("MAIC (means)", "MAIC (means and SDs)", "MAIC (all reported summaries)", "sMAIC",
                          "STC", "G-computation", "G-MAIC (means and SDs)", "G-MAIC", "sG-MAIC"), stringsAsFactors = FALSE)
M$key <- paste(M$method, M$R); cols <- setNames(npg[c(10, 2, 4, 1, 9, 6, 7, 3, 8)], M$label)
shapes_em <- c(moderate = 1, strong = 16)                                   # 空心 = 中等档，实心 = 强档
D_lab <- c(D1 = "Normal", D2 = "Right-skewed", D3 = "Left-skewed", D4 = "Truncated comparator")
S_lab <- c(linear = "Linear effect modification", quad = "U-shaped effect modification", threshold = "Threshold effect modification")
S_short <- c(linear = "Linear", quad = "U-shaped", threshold = "Threshold")
pick <- function(df) { df$key <- paste(df$method, df$R); df <- df[df$key %in% M$key, ]; df$label <- factor(M$label[match(df$key, M$key)], levels = rev(M$label)); df }
theme_pub <- theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8), legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        axis.text = element_text(size = 7.5), axis.title = element_text(size = 9), plot.tag = element_text(face = "bold", size = 11))

# ---- 图 1 / 图 2：偏倚点图网格 ------------------------------------------------------------------------------
dotgrid <- function(ep, an, xlab, tag) {
  d <- pick(perf[perf$endpoint == ep & perf$anchored == an, ])
  d$em <- factor(ifelse(d$em_mult == 1, "moderate", "strong"), levels = c("moderate", "strong"))
  d$lo <- d$bias - 1.96 * d$bias_mcse; d$hi <- d$bias + 1.96 * d$bias_mcse
  d$D <- factor(D_lab[d$D], levels = D_lab); d$shape <- factor(S_lab[d$shape], levels = S_lab)
  # 未运行的方法（本研究中只有无锚定生存的 STC）：整行 NA，在图上标 "not run"
  na_rows <- unique(d[is.na(d$bias), c("D", "shape", "label")]); d <- d[!is.na(d$bias), ]
  # 横轴范围由 STC 以外的方法决定；超出的点截到边界并标数值
  core <- d[d$method != "M2", ]; lim <- range(c(core$lo, core$hi)); lim <- lim + c(-1, 1) * 0.08 * diff(lim); lim <- c(min(lim[1], -0.02), max(lim[2], 0.02))
  d$off <- d$bias < lim[1] | d$bias > lim[2]; d$bias_plot <- pmin(pmax(d$bias, lim[1]), lim[2])
  d$lo_plot <- ifelse(d$off, NA, pmax(d$lo, lim[1])); d$hi_plot <- ifelse(d$off, NA, pmin(d$hi, lim[2]))
  # 每行一条“茎”：从 0 到两个点中最远的一个（两点异号时覆盖两边）
  stem <- do.call(rbind, lapply(split(d, list(d$D, d$shape, d$label), drop = TRUE), function(s)
    data.frame(D = s$D[1], shape = s$shape[1], label = s$label[1], x0 = min(0, s$bias_plot), x1 = max(0, s$bias_plot))))
  # 超出范围的点：数值写在该行的另一端（那一端通常是空的，不与行内的点相撞），前面带点形符号（○ 中等 / ● 强）
  offlab <- NULL
  if (any(d$off)) offlab <- do.call(rbind, lapply(split(d[d$off, ], list(d$D[d$off], d$shape[d$off], d$label[d$off], d$bias[d$off] < lim[1]), drop = TRUE), function(s) {
    s <- s[order(s$em), ]; left <- s$bias[1] < lim[1]
    data.frame(D = s$D[1], shape = s$shape[1], label = s$label[1], x = if (left) lim[2] else lim[1], hj = if (left) 1.05 else -0.05,
               txt = paste(sprintf("%s %.2f", ifelse(s$em == "moderate", "○", "●"), s$bias), collapse = "  ")) }))
  p <- ggplot(d, aes(y = label, colour = label)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_segment(data = stem, aes(x = x0, xend = x1, yend = label), linewidth = 0.35, alpha = 0.6) +
    geom_segment(aes(x = lo_plot, xend = hi_plot, yend = label), linewidth = 1.1, alpha = 0.35, na.rm = TRUE, lineend = "butt") +
    geom_point(data = d[d$em == "strong", ], aes(x = bias_plot, shape = em), size = 1.7) +
    geom_point(data = d[d$em == "moderate", ], aes(x = bias_plot, shape = em), size = 2.5, stroke = 0.65) +
    facet_grid(D ~ shape) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_shape_manual(values = shapes_em, name = "Effect-modification strength", labels = c("moderate", "strong"), drop = FALSE) +
    guides(shape = guide_legend(override.aes = list(size = c(2.5, 1.7), stroke = c(0.65, 0.5)))) +
    scale_y_discrete(drop = FALSE) +                                       # 保留整行 NA 的方法（STC）的位置
    coord_cartesian(xlim = lim) +
    labs(x = xlab, y = NULL, tag = tag) + theme_pub
  if (!is.null(offlab)) p <- p + geom_text(data = offlab, aes(x = x, label = txt, hjust = hj), size = 2.0, family = "Arial", colour = "grey20")
  if (nrow(na_rows)) p <- p + geom_text(data = na_rows, aes(x = 0, label = "not run"), hjust = -0.3, size = 2.1, family = "Arial", colour = "grey45", fontface = "italic")
  p
}
f1a <- dotgrid("binary", TRUE, "Bias in the marginal log odds ratio", "a")
f1b <- dotgrid("tte", TRUE, "Bias in the marginal log hazard ratio", "b")
f1 <- f1a / f1b + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("output/figures/fig2_anchored_bias_dotgrid.png", f1, width = 7.2, height = 12.6, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave("output/figures/fig2_anchored_bias_dotgrid.pdf", f1, width = 7.2, height = 12.6, device = cairo_pdf)
f2a <- dotgrid("binary", FALSE, "Bias in the marginal log odds ratio", "a")
f2b <- dotgrid("tte", FALSE, "Bias in the 24-month RMST difference (months)", "b")
f2 <- f2a / f2b + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("output/figures/fig3_unanchored_bias_dotgrid.png", f2, width = 7.2, height = 12.6, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave("output/figures/fig3_unanchored_bias_dotgrid.pdf", f2, width = 7.2, height = 12.6, device = cairo_pdf)

# ---- 图 3：配对 RMSE 比瓦片图 ---------------------------------------------------------------------------------
tiles <- function(ep, tag) {
  d <- pick(perf[perf$endpoint == ep & perf$method != "M1m", ])
  ref <- perf[perf$endpoint == ep & perf$method == "M1m" & perf$R == "R1", c("dgm_id", "rmse")]; names(ref)[2] <- "rmse_ref"
  d <- merge(d, ref, by = "dgm_id")
  pr <- paired[paired$endpoint == ep, c("dgm_id", "R", "method", "dmse_vs_means", "dmse_vs_means_mcse")]
  d <- merge(d, pr, by = c("dgm_id", "R", "method"))
  d$ratio <- d$rmse / d$rmse_ref
  d$ns <- abs(d$dmse_vs_means) < 1.96 * d$dmse_vs_means_mcse                  # 配对 MSE 差的 95% MC 区间含 0 ⇔ RMSE 比的区间含 1
  d$lowcov <- !is.na(d$cover) & d$cover < 0.925                         # 决策 39：无锚定 RMST 差改用 IJ 稳健标准误（sims/41）后恢复显示
  row_lv <- as.vector(outer(c("moderate", "strong"), S_short, function(e, s) sprintf("%s, %s", s, e)))   # 行 = 形状 × 强度（形状为主）
  d$row <- factor(sprintf("%s, %s", S_short[d$shape], ifelse(d$em_mult == 1, "moderate", "strong")), levels = rev(row_lv))
  d$D <- factor(D_lab[d$D], levels = D_lab)                                                              # 分面行 = 协变量分布
  d$design <- factor(ifelse(d$anchored, "Anchored", "Unanchored"), levels = c("Anchored", "Unanchored"))
  d$label <- factor(as.character(d$label), levels = M$label[-1])
  d$fill_val <- ifelse(d$ns | is.na(d$ratio), NA, pmin(pmax(d$ratio, 0.6), 1.6))
  d$txt <- ifelse(is.na(d$ratio), "not run", sprintf("%.2f", d$ratio))
  ggplot(d, aes(x = label, y = row)) +
    geom_tile(aes(fill = fill_val), colour = "white", linewidth = 0.5) +
    geom_text(aes(label = txt), size = 2.0, family = "Arial", colour = "grey15") +
    geom_point(data = d[d$lowcov, ], aes(x = label, y = row), shape = 17, size = 1.3, colour = "black", position = position_nudge(x = 0.38, y = 0.3)) +
    facet_grid(D ~ design) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", midpoint = 1, limits = c(0.6, 1.6), na.value = "grey85",
                         name = "RMSE relative to MAIC (means)") +
    scale_x_discrete(position = "top") +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    labs(x = NULL, y = NULL, tag = tag) + theme_pub +
    theme(axis.text.x.top = element_text(angle = 35, hjust = 0, size = 7), axis.text.y = element_text(size = 7), panel.grid = element_blank(), strip.text.y = element_text(size = 7.5),
          legend.key.width = unit(1.6, "cm"), plot.margin = margin(5, 28, 5, 5))
}
f3 <- tiles("binary", "a") / tiles("tte", "b") + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("output/figures/fig4_paired_rmse_ratio_tiles.png", f3, width = 7.6, height = 12.5, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave("output/figures/fig4_paired_rmse_ratio_tiles.pdf", f3, width = 7.6, height = 12.5, device = cairo_pdf)
cat("figures written to output/figures\n")
