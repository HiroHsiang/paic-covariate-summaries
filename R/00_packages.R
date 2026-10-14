# 本项目（本模拟研究，protocol v0.2 第 7 节）所需的 R 包清单。
# 可打包代码：不用 setwd()，不写绝对路径。换机器时 source 本文件后调用 ensure_packages()；
# 正式运行前调用 write_env_manifest() 记录各包版本，供论文投稿时的可重复性审查使用。

# 核心：数据生成、估计量、性能评价 ------------------------------------------
pkgs_core <- c(
  "survival",     # Cox / Weibull 拟合，RMST
  "flexsurv",     # 参数生存模型，RMST 辅助函数
  "copula",       # 高斯 copula，生成相关的协变量
  "sn",           # 偏正态分布，用于由汇总量重建边际分布
  "rmetalog",     # metalog 分位数参数化分布，重建边际分布的另一选择
  "estmeansd",    # 中位数/IQR/范围 → 均值/SD 换算（McGrath 2020）
  "sandwich",     # 稳健方差的组件
  "numDeriv",     # sandwich 估计量所需的数值雅可比矩阵
  "matrixStats",  # 行列汇总的快速实现
  "rsimsum",      # 性能指标及其蒙特卡洛标准误
  "foreach", "doParallel", "doRNG",  # 可重复的并行重复模拟
  "data.table", "dplyr", "ggplot2"
)

# 对照方法的现成实现与交叉核对 ----------------------------------------------
pkgs_methods <- c(
  "maicplus",     # MAIC 交叉核对（hta-pharma）
  "outstandR",    # G-computation 交叉核对
  "optweight"     # 带容差的平衡权重（sMAIC 交叉核对；注意 tols 默认按 SD 单位，需 std.cont = FALSE）
)

# 开发 ------------------------------------------------------------------------
pkgs_dev <- c("testthat")

pkgs_all <- unique(c(pkgs_core, pkgs_methods, pkgs_dev))

#' 安装 `pkgs` 中尚未安装的包（Windows 下用 CRAN 二进制版，不需要 Rtools）。
ensure_packages <- function(pkgs = pkgs_all, repos = "https://cloud.r-project.org") {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) {
    message("正在安装：", paste(missing, collapse = ", "))
    install.packages(missing, repos = repos,
                     type = if (.Platform$OS.type == "windows") "binary" else getOption("pkgType"))
  }
  still <- setdiff(pkgs, rownames(installed.packages()))
  if (length(still)) warning("未能安装：", paste(still, collapse = ", "))
  invisible(setdiff(pkgs, still))
}

#' 把 R 版本、平台、核数和各包版本写到 `path`（作为运行清单的一部分）。
write_env_manifest <- function(path, pkgs = pkgs_all) {
  ver <- vapply(pkgs, function(p)
    if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "MISSING",
    character(1))
  lines <- c(
    paste0("# Environment manifest (", format(Sys.time(), "%Y-%m-%d %H:%M"), ")"),
    "", R.version.string, paste("Platform:", R.version$platform),
    paste("Cores:", parallel::detectCores()), "",
    "| package | version |", "|---|---|",
    sprintf("| %s | %s |", names(ver), ver)
  )
  writeLines(lines, path)
  invisible(ver)
}
