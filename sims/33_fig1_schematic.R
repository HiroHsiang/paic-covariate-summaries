# 图 1（示意图，2026-10-12）：发表的基线表报告了什么、每种方法从中取什么。全文的“地图”与词汇表。
#   a：单个连续协变量在对照试验基线表里的五种报告方式（模拟“表 1”中的几行，数值固定：N = 300，
#      均值 53.0（SD 12.0），中位数 53.1 [45.0, 61.3]，范围 21.4–88.9，切点 65 以上 47 (15.7%)）
#   b：在协变量分布（正态，均值 53、SD 12，与上述数值一致）上，各汇总量分别刻画什么：
#      均值 = 位置（峰顶上方的倒三角），SD = 离散度（拐点到中心的尺寸线，长 1 SD），四分位数 = 横轴上的三个切点，
#      切点以上比例 = 阴影尾部，范围 = 支撑的两端（尺寸线）。
#      注意：53 + 12 = 65 恰好等于切点，若画 ±1 SD 双箭头，右端会落在尾部的角上、像是“指向”比例，故只画左侧一段。
#   c：方法 × 汇总量网格。实心圆 = 精确匹配；空心圆 = 在容差内匹配；浅绿底格 = 用于重建对照协变量分布；短横 = 未使用。
#      底格与圆可以叠加（G-MAIC、Soft G-MAIC：既用于重建，又被权重匹配）。
# 样式同 sims/31_figures_main.R：ggplot2 + patchwork，theme_bw(9, Arial)，NPG 色板，条带 grey92，标签小写粗体。
# a、c 两个面板的数据坐标按英寸设定（面板实宽：a ≈ 3.6 in，c ≈ 6.85 in，由 7.2 in 总宽与边距决定），便于按磅值放文字。
# 所有文字 ≥ 7 pt。输出：output/figures/fig1_table_to_constraints.png（300 dpi）与 .pdf。
# 用法（项目根目录）：Rscript sims/33_fig1_schematic.R

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
dir.create("output/figures", showWarnings = FALSE)

npg <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
ink <- "grey15"                     # 文字与线条
col_exact <- npg[4]                 # 深蓝：匹配（精确 / 容差）
col_recon_fill <- npg[7]            # 浅绿：用于重建
col_recon_line <- npg[3]
col_tail <- npg[1]                  # 尾部阴影
fs <- 7.5 / .pt                     # 注释字号 7.5 pt
fam <- "Arial"

theme_pub <- theme_bw(base_size = 9, base_family = fam) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8), legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        axis.text = element_text(size = 7.5), axis.title = element_text(size = 9), plot.tag = element_text(face = "bold", size = 11))
# 示意图面板：去掉坐标轴、网格与面板边框，只保留条带作为面板标题
theme_schem <- theme_pub + theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
                                 panel.grid = element_blank(), panel.border = element_blank(), panel.background = element_blank(),
                                 strip.text = element_text(size = 8, hjust = 0, margin = margin(3, 4, 3, 4)),
                                 plot.margin = margin(2, 6, 2, 2))
base <- function(title) ggplot(data.frame(ttl = title)) + facet_wrap(~ ttl)   # 单面板 + 条带标题
txt <- function(...) annotate("text", ..., size = fs, family = fam, colour = ink)
arw <- function(ends) arrow(ends = ends, length = unit(0.045, "in"), angle = 22, type = "closed")

# ---- a：基线表中的五种报告方式 ------------------------------------------------------------------------------
W_a <- 3.6                                                          # 面板实宽（英寸）
tab <- data.frame(
  pattern = c("Mean (SD)", "Median [Q1, Q3]", "Median (range)", "Mean (SD), n (%) above cut-point", "All of the above"),
  value = c("53.0 (12.0)", "53.1 [45.0, 61.3]", "53.1 (21.4–88.9)", "53.0 (12.0); 47 (15.7%) above 65",
            "53.0 (12.0); 53.1 [45.0, 61.3];\nrange 21.4–88.9; 47 (15.7%) above 65"),
  y = c(-0.37, -0.545, -0.72, -0.895, -1.125))
x_val <- 1.76                                                       # 数值列起点（最长行标签 1.55 in）
rules <- data.frame(y = c(-0.24, -1.29), lw = c(0.3, 0.45))         # 表头下线、底线（条带已充当顶线，不再画顶线）
pa <- base("Ways a trial publication may report one covariate") +
  geom_segment(data = rules, aes(x = 0, xend = W_a, y = y, yend = y, linewidth = lw), colour = ink) +
  scale_linewidth_identity() +
  annotate("text", x = c(0.04, x_val), y = -0.12, label = c("Reporting pattern", "Comparator trial (N = 300)"),
           hjust = 0, size = fs, family = fam, colour = ink, fontface = "bold") +
  geom_text(data = tab, aes(x = 0.04, y = y, label = pattern), hjust = 0, size = fs, family = fam, colour = ink) +
  geom_text(data = tab, aes(x = x_val, y = y, label = value), hjust = 0, size = fs, family = fam, colour = ink, lineheight = 0.95) +
  coord_cartesian(xlim = c(0, W_a), ylim = c(-1.45, 0.04), expand = FALSE, clip = "off") +
  labs(tag = "a") + theme_schem

# ---- b：每个汇总量在分布上刻画什么 -----------------------------------------------------------------------------
mu <- 53; s <- 12; lo <- 21.4; hi <- 88.9; qs <- c(45.0, 53.1, 61.3); cut <- 65
f <- function(x) exp(-0.5 * ((x - mu) / s)^2)                       # 正态密度，峰高归一为 1
xs <- seq(lo, hi, length.out = 500); xt <- seq(cut, hi, length.out = 300)
dens <- data.frame(x = c(lo, xs, hi), y = c(0, f(xs), 0))
tail_df <- data.frame(x = c(cut, xt, hi), y = c(0, f(xt), 0))
y_q <- -0.15; y_r <- -0.34                                          # 横轴下方两行：四分位数标签、范围尺寸线
pb <- base("Feature of the distribution captured by each summary") +
  geom_polygon(data = dens, aes(x, y), fill = "grey93", colour = NA) +
  geom_polygon(data = tail_df, aes(x, y), fill = col_tail, alpha = 0.45, colour = NA) +            # 切点以上比例
  geom_line(data = dens[-c(1, nrow(dens)), ], aes(x, y), colour = ink, linewidth = 0.45) +
  annotate("segment", x = 12, xend = 98, y = 0, yend = 0, colour = ink, linewidth = 0.35) +      # 横轴（协变量取值）
  # 四分位数：分布内部的虚线 + 横轴上的刻度
  annotate("segment", x = qs, xend = qs, y = 0, yend = f(qs), colour = "grey35", linewidth = 0.35, linetype = "22") +
  annotate("segment", x = qs, xend = qs, y = -0.06, yend = 0, colour = ink, linewidth = 0.45) +
  txt(x = qs[2], y = y_q, label = "Quartiles") +
  # 均值：峰顶上方的倒三角（位置）
  annotate("point", x = mu, y = 1.08, shape = 25, size = 1.9, fill = ink, colour = ink) +
  txt(x = mu, y = 1.22, label = "Mean") +
  # SD：从左拐点（均值 − 1 SD）到中心的尺寸线（离散度）
  annotate("segment", x = mu - s, xend = mu, y = f(mu - s), yend = f(mu - s), colour = ink, linewidth = 0.4, arrow = arw("both")) +
  txt(x = mu - s - 1.5, y = f(mu - s), label = "SD", hjust = 1) +
  # 切点与尾部
  annotate("segment", x = cut, xend = cut, y = 0, yend = 0.86, colour = ink, linewidth = 0.35, linetype = "42") +
  txt(x = cut + 1.2, y = 0.9, label = "Proportion above\ncut-point", hjust = 0, vjust = 1, lineheight = 0.95) +
  # 范围：支撑两端的延伸线 + 尺寸线（中间留空写标签）
  annotate("segment", x = c(lo, hi), xend = c(lo, hi), y = y_r - 0.05, yend = 0.08, colour = "grey35", linewidth = 0.3) +
  annotate("segment", x = lo, xend = 28.5, y = y_r, yend = y_r, colour = ink, linewidth = 0.4, arrow = arw("first")) +
  annotate("segment", x = 37.5, xend = hi, y = y_r, yend = y_r, colour = ink, linewidth = 0.4, arrow = arw("last")) +
  txt(x = 33, y = y_r, label = "Range") +
  coord_cartesian(xlim = c(12, 98), ylim = c(-0.42, 1.3), expand = FALSE, clip = "off") +
  labs(tag = "b") + theme_schem

# ---- c：每种方法使用哪些汇总量、怎样使用 -------------------------------------------------------------------------
W_c <- 6.85                                                         # 面板实宽（英寸）
meth <- c("Unadjusted", "MAIC (means)", "MAIC (means and SDs)", "MAIC (all reported summaries)",
          "MAIC (means, quartiles and proportions)", "sMAIC", "STC", "G-computation", "G-MAIC (means and SDs)", "G-MAIC", "sG-MAIC")
summ <- c("Mean", "SD", "Quartiles", "Range", "Proportion\nabove\ncut-point")
# 权重 / 中心化如何使用各汇总量："E" = 精确匹配，"T" = 在容差内匹配，"" = 不用于匹配
# （STC 的均值用于把结局模型的协变量中心化，按“精确”记；范围把 IPD 截到报告的支撑内，按“精确”记——sMAIC 与 soft G-MAIC 也是如此，
#  它们只有 SD 和四分位数带容差（R/06_maic.R：soft = c("sd", "quantiles")）。2026-10-13 决策 37 据代码改正，并加入 G-MAIC (means and SDs)）
use <- rbind(c("",  "",  "",  "",  ""),     # Unadjusted
             c("E", "",  "",  "",  ""),     # MAIC (means)
             c("E", "E", "",  "",  ""),     # MAIC (means and SDs)
             c("E", "E", "E", "E", "E"),    # MAIC (all reported summaries)
             c("E", "",  "E", "E", "E"),    # MAIC (means, quartiles and proportions)
             c("E", "T", "T", "E", "E"),    # sMAIC（范围为精确截取）
             c("E", "",  "",  "",  ""),     # STC
             c("",  "",  "",  "",  ""),     # G-computation
             c("E", "E", "",  "",  ""),     # G-MAIC (means and SDs)
             c("E", "E", "E", "E", "E"),    # G-MAIC
             c("E", "T", "T", "E", "E"))    # Soft G-MAIC（范围为精确截取）
# 是否用于重建对照试验的协变量分布
# 重建只用均值、SD、四分位数和范围（R/07 reconstruct_marginal），不用切点以上比例；G-computation 与 G-MAIC (means and SDs) 用正态（均值、SD）
recon <- rbind(matrix(FALSE, 7, 5), c(TRUE, TRUE, FALSE, FALSE, FALSE), c(TRUE, TRUE, FALSE, FALSE, FALSE),
               matrix(c(TRUE, TRUE, TRUE, TRUE, FALSE), 2, 5, byrow = TRUE))
pitch <- 0.19; gap <- 0.07                                          # 行距与组间额外间隔
y_row <- numeric(length(meth)); cur <- 0
for (i in seq_along(meth)) { y_row[i] <- cur; cur <- cur - pitch - if (i %in% c(1, 6, 8)) gap else 0 }
x_lab <- 1.98; x_col <- 2.33 + 0.5 * 0:4                            # 行标签右端（最长 1.87 in）与五列中心
cells <- expand.grid(i = seq_along(meth), j = seq_along(summ))
cells$x <- x_col[cells$j]; cells$y <- y_row[cells$i]
cells$use <- use[cbind(cells$i, cells$j)]; cells$recon <- recon[cbind(cells$i, cells$j)]
tile_w <- 0.34; tile_h <- 0.15
y_head <- 0.2                                                       # 表头底边
y_top <- y_head + 0.45                                              # 三行表头顶端（用于纵轴范围）
sep_y <- (y_row[c(1, 6, 8)] + y_row[c(2, 7, 9)]) / 2                # 组间分隔线
x_end <- x_col[5] + 0.33
# 图例（面板内右侧，顶端与第一行方法对齐）
x_lg <- x_end + 0.36
y_lg <- y_row[1] - c(0, 0.24, 0.52, 0.82)
pc <- base("Summaries used by each method") +
  annotate("segment", x = 0, xend = x_end, y = c(y_head - 0.08, min(y_row) - 0.13),            # 表头下线、底线
           yend = c(y_head - 0.08, min(y_row) - 0.13), colour = ink, linewidth = c(0.3, 0.45)) +
  annotate("segment", x = 0, xend = x_end, y = sep_y, yend = sep_y, colour = "grey80", linewidth = 0.3) +
  annotate("text", x = x_col, y = y_head, label = summ, vjust = 0, size = fs, family = fam, colour = ink, fontface = "bold", lineheight = 0.95) +
  annotate("text", x = x_lab, y = y_head, label = "Method", hjust = 1, vjust = 0, size = fs, family = fam, colour = ink, fontface = "bold") +
  annotate("text", x = x_lab, y = y_row, label = meth, hjust = 1, size = fs, family = fam, colour = ink) +
  geom_tile(data = cells[cells$recon, ], aes(x, y), width = tile_w, height = tile_h, fill = col_recon_fill, colour = col_recon_line,
            linewidth = 0.25, alpha = 0.55) +
  geom_point(data = cells[cells$use == "E", ], aes(x, y), shape = 21, size = 2.4, fill = col_exact, colour = col_exact, stroke = 0.6) +
  geom_point(data = cells[cells$use == "T", ], aes(x, y), shape = 21, size = 2.4, fill = "white", colour = col_exact, stroke = 0.75) +
  geom_text(data = cells[cells$use == "" & !cells$recon, ], aes(x, y), label = "–", size = fs, family = fam, colour = "grey55") +
  # 图例
  annotate("point", x = x_lg, y = y_lg[1], shape = 21, size = 2.4, fill = col_exact, colour = col_exact, stroke = 0.6) +
  annotate("point", x = x_lg, y = y_lg[2], shape = 21, size = 2.4, fill = "white", colour = col_exact, stroke = 0.75) +
  annotate("tile", x = x_lg, y = y_lg[3], width = tile_w, height = tile_h, fill = col_recon_fill, colour = col_recon_line,
           linewidth = 0.25, alpha = 0.55) +
  annotate("text", x = x_lg, y = y_lg[4], label = "–", size = fs, family = fam, colour = "grey55") +
  annotate("text", x = x_lg + 0.25, y = y_lg, vjust = c(0.5, 0.5, 0.78, 0.5), hjust = 0, size = fs, family = fam, colour = ink, lineheight = 0.95,
           label = c("Matched exactly", "Matched within a tolerance", "Used to reconstruct the comparator\ncovariate distribution", "Not used")) +
  coord_cartesian(xlim = c(0, W_c), ylim = c(min(y_row) - 0.2, y_top + 0.04), expand = FALSE, clip = "off") +
  labs(tag = "c") + theme_schem

# ---- 拼图与输出 -------------------------------------------------------------------------------------------------
top <- pa + pb + plot_layout(widths = c(3.6, 3.04))
fig <- top / pc + plot_layout(heights = c(1.5, 2.94))                # c 面板比原来多一行（+0.19 in）
H <- 4.94
ggsave("output/figures/fig1_table_to_constraints.png", fig, width = 7.2, height = H, dpi = 300, bg = "white", device = ragg::agg_png)
ggsave("output/figures/fig1_table_to_constraints.pdf", fig, width = 7.2, height = H, device = cairo_pdf)
cat("written output/figures/fig1_table_to_constraints.{png,pdf}\n")
