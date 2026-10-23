#==============================================#
#       生成 Gamma DGP 下的 ELBH/TBH 结果     #
#==============================================#
rm(list = ls())

if (!requireNamespace("emplik", quietly = TRUE)) {
  stop("需要安装 R 包 'emplik'。", call. = FALSE)
}

#==============================================#
#                 辅助函数                     #
#==============================================#
# 上侧 p 值：EL 用经验似然比统计量，失败返回 NA
el_pvalue <- function(z, delta0) {
  if (mean(z) <= delta0) {
    return(1)
  }
  fit <- tryCatch(
    suppressWarnings(emplik::el.test(z, mu = delta0)),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$`-2LLR`)) {
    return(NA_real_)
  }
  0.5 * stats::pchisq(max(0, as.numeric(fit$`-2LLR`)), df = 1, lower.tail = FALSE)
}

# 上侧 p 值：T 用单样本 t 检验
t_pvalue <- function(z, delta0) {
  n <- length(z)
  stats::pt((mean(z) - delta0) / (stats::sd(z) / sqrt(n)), df = n - 1, lower.tail = FALSE)
}

# Benjamini–Hochberg 拒绝集（返回 logical 向量）
bh_reject <- function(p, alpha) {
  m <- length(p)
  ord <- order(p)
  k <- which(p[ord] <= seq_len(m) * alpha / m)
  if (length(k) == 0L) {
    return(rep(FALSE, m))
  }
  p <= p[ord][max(k)]
}

#==============================================#
#                 数据保存                     #
#==============================================#
save_simulation_data <- function(results, output_dir = ".") {
  ms_table <- sort(unique(results$m))
  paths <- character(0)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (m in ms_table) {
    table_data <- results[results$m == m, , drop = FALSE]
    file_path <- file.path(output_dir, paste0("Table2-Gamma-FFR-Power-m-", m, ".csv"))
    utils::write.csv(table_data, file_path, row.names = FALSE)
    paths <- c(paths, file_path)
    message("数据已保存：", file_path)
  }

  invisible(paths)
}

#==============================================#
#                 结果作图                     #
#==============================================#
open_eps_device <- function(file_path, width, height) {
  if (tolower(tools::file_ext(file_path)) != "eps") {
    stop("图片仅保存为 EPS。", call. = FALSE)
  }
  grDevices::postscript(file_path, width = width, height = height, paper = "special",
                        horizontal = FALSE, onefile = FALSE, family = "Helvetica")
}

# 每个 pi1 一个面板：pi1=0 只画 FFR；pi1=0.5 与 pi1=1 分别并列 FFR 与 Power
draw_ffr_power <- function(results, m, alpha = 0.05) {
  plot_data <- results[results$m == m, , drop = FALSE]
  if (nrow(plot_data) == 0L) {
    stop("没有与 m 匹配的作图数据。", call. = FALSE)
  }
  methods_plot <- intersect(c("ELBH", "TBH"), unique(plot_data$method))
  if (length(methods_plot) < 2L) {
    stop("同图对比需要同时包含 ELBH 和 TBH 的结果。", call. = FALSE)
  }

  pi1_values <- sort(unique(plot_data$target_pi1))
  ns_plot <- sort(unique(plot_data$n))
  panel_specs <- do.call(rbind, lapply(pi1_values, function(pi1) {
    pi1_data <- plot_data[plot_data$target_pi1 == pi1, , drop = FALSE]
    # m0 = m（即 m1 = 0，全部为原假设）时只画 FFR；否则 FFR 与 Power 并列
    is_null <- all(pi1_data$m1 == 0L)
    metrics <- if (is_null) "ffr" else c("ffr", "power")
    data.frame(target_pi1 = pi1, metric = metrics, stringsAsFactors = FALSE)
  }))
  method_colors <- c(ELBH = "red", TBH = "blue")
  method_points <- c(ELBH = 16, TBH = 17)
  n_line_types <- rep(c(1, 2, 3, 4, 5, 6), length.out = length(ns_plot))
  style_grid <- expand.grid(
    method = methods_plot, n = ns_plot,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  style_grid$color <- unname(method_colors[style_grid$method])
  style_grid$point <- unname(method_points[style_grid$method])
  style_grid$line_type <- n_line_types[match(style_grid$n, ns_plot)]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1L, nrow(panel_specs)), oma = c(1.1, 1.3, 2.2, 0.4),
    mar = c(3.2, 3.5, 2.2, 0.8), mgp = c(2.1, 0.6, 0), las = 1,
    cex.main = 1.05, cex.lab = 1, cex.axis = 0.9
  )

  for (panel_index in seq_len(nrow(panel_specs))) {
    pi1 <- panel_specs$target_pi1[panel_index]
    metric <- panel_specs$metric[panel_index]
    panel <- plot_data[plot_data$target_pi1 == pi1, , drop = FALSE]
    x_range <- range(panel$effect_size)
    if (diff(x_range) == 0) {
      x_range <- x_range + c(-0.01, 0.01)
    }
    y_limit <- c(0, 1)
    y_label <- if (metric == "ffr") "False flagging rate" else "Power"

    graphics::plot(NA, xlim = x_range, ylim = y_limit, xaxs = "i", yaxs = "i",
                   xaxt = "n", yaxt = "n", xlab = "", ylab = "",
                   main = sprintf("target pi1 = %.2f", pi1))
    graphics::axis(1, at = sort(unique(panel$effect_size)), font = 2, lwd = 1.5, tcl = 0.35)
    graphics::axis(2, at = pretty(y_limit, n = 5), font = 2, lwd = 1.5, tcl = 0.35)
    graphics::axis(3, labels = FALSE, lwd = 1.5, tcl = 0)
    graphics::axis(4, labels = FALSE, lwd = 1.5, tcl = 0)
    graphics::title(xlab = expression(delta), line = 2,   font.lab = 2)
    graphics::title(ylab = y_label,           line = 2.2, font.lab = 2)
    graphics::grid(col = "gray88", lty = 1)
    if (metric == "ffr") {
      graphics::abline(h = alpha, col = "gray35", lty = 2, lwd = 1.5)
    }

    for (j in seq_len(nrow(style_grid))) {
      curve <- panel[
        panel$method == style_grid$method[j] & panel$n == style_grid$n[j],
        , drop = FALSE
      ]
      curve <- curve[order(curve$effect_size), , drop = FALSE]
      graphics::lines(curve$effect_size, curve[[metric]], col = style_grid$color[j],
                      lwd = 2, lty = style_grid$line_type[j])
      graphics::points(curve$effect_size, curve[[metric]], col = style_grid$color[j],
                       pch = style_grid$point[j], cex = 0.9)
    }

    if (panel_index == 1L) {
      legend_labels <- paste0(style_grid$method, ", n = ", style_grid$n)
      legend_colors <- style_grid$color
      legend_lty <- style_grid$line_type
      legend_pch <- style_grid$point
      if (metric == "ffr") {
        legend_labels <- c(legend_labels, substitute(alpha == value, list(value = alpha)))
        legend_colors <- c(legend_colors, "gray35")
        legend_lty <- c(legend_lty, 2)
        legend_pch <- c(legend_pch, NA)
      }
      graphics::legend("topright", legend = legend_labels, col = legend_colors,
                       lty = legend_lty, pch = legend_pch, lwd = 1.8,
                       bg = "white", cex = 0.76)
    }
  }

  graphics::mtext(sprintf("ELBH vs TBH under Gamma DGP (m = %d)", m),
                  side = 3, outer = TRUE, line = 0.5, font = 2, cex = 1.15)
  invisible(NULL)
}

save_ffr_power_plot <- function(results, m, alpha, file_path,
                                width = 13, height = 4.8) {
  open_eps_device(file_path = file_path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_ffr_power(results = results, m = m, alpha = alpha)
  invisible(file_path)
}

#==============================================#
#                 模拟主函数                   #
#==============================================#
run_simulation <- function(ns, m, pi1s, effect_sizes, nsim, shape = 2,
                           delta0 = 1.25, alpha = 0.05,
                           methods = c("EL", "T"),
                           save_data = FALSE, save_plot = FALSE, plot_dir = ".") {
  grid <- expand.grid(n = ns, m = m, pi1 = pi1s,
                      effect_size = effect_sizes,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)

  all_rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    s <- grid[i, ]
    message(sprintf(
      "[%d/%d] n=%d, m=%d, pi1=%.2f, delta=%.2f",
      i, nrow(grid), s$n, s$m, s$pi1, s$effect_size
    ))

    # ---- 构建设置：备择组位置与各组 disparity ----
    m1 <- as.integer(round(s$m * s$pi1))
    alternative <- rep(FALSE, s$m)
    if (m1 > 0L) {
      # 将备择组均匀放在分割区间中，避免总固定在 X 较小的一端
      alternative[unique(as.integer(round(seq(1, s$m, length.out = m1))))] <- TRUE
    }
    disparity <- rep(delta0, s$m)
    disparity[alternative] <- delta0 + s$effect_size
    null_true <- !alternative

    # ---- 重复模拟，记录每次有效重复的指标 ----
    ffr <- matrix(NA_real_, nrow = nsim, ncol = length(methods), dimnames = list(NULL, methods))
    power <- matrix(NA_real_, nrow = nsim, ncol = length(methods), dimnames = list(NULL, methods))

    completed <- 0L
    progress_every <- max(1L, floor(nsim / 10L))

    while (completed < nsim) {
      # 生成数据
      x <- stats::runif(s$n)
      group <- pmin.int(floor(s$m * x) + 1L, s$m)
      # E(M | G_j) = disparity_j，且 Var(M | G_j) = disparity_j^2 / shape > 0
      performance <- stats::rgamma(s$n, shape = shape, scale = disparity[group] / shape)

      # 按组计算各方法的上侧 p 值
      p_values <- matrix(NA_real_, length(methods), s$m, dimnames = list(methods, NULL))
      ok <- TRUE
      for (g in seq_len(s$m)) {
        z <- performance[group == g]
        if (length(z) < 2L || !all(is.finite(z)) || stats::sd(z) <= 0) {
          ok <- FALSE
          break
        }
        for (method in methods) {
          p_values[method, g] <- switch(method, EL = el_pvalue(z, delta0), T = t_pvalue(z, delta0))
        }
        if (any(!is.finite(p_values[, g]))) {
          ok <- FALSE
          break
        }
      }
      if (!ok) {
        next
      }

      # 对每种方法做 BH，并记录 FFR / 功效
      completed <- completed + 1L
      for (method in methods) {
        rejected <- bh_reject(p_values[method, ], alpha)
        total <- sum(rejected)
        ffr[completed, method] <- sum(rejected & null_true) / max(total, 1L)
        if (m1 > 0L) {
          power[completed, method] <- sum(rejected & !null_true) / m1
        }
      }

      if (nsim >= 100L && completed %% progress_every == 0L) {
        message(sprintf("  已完成 %d/%d 次有效重复", completed, nsim))
      }
    }

    # ---- 把当前设置汇总成每个方法一行 ----
    all_rows[[i]] <- do.call(rbind, lapply(methods, function(method) {
      method_power <- power[, method]
      data.frame(
        method = paste0(method, "BH"),
        effect_size = s$effect_size,
        n = s$n, m1 = m1, m = s$m,
        realized_pi1 = m1 / s$m, target_pi1 = s$pi1,
        ffr = mean(ffr[, method]),
        power = if (any(!is.na(method_power))) mean(method_power) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }

  result <- do.call(rbind, all_rows)
  rownames(result) <- NULL

  # ---- 保存、打印并作图 ----
  if (save_data) {
    save_simulation_data(result, output_dir = plot_dir)
  }

  cat(sprintf("\n delta0 = %.2f, shape = %g, alpha = %.2f, nsim = %d, m = %d",
              delta0, shape, alpha, nsim, m), "\n")
  print(result)

  if (save_plot) {
    file_path <- file.path(plot_dir, paste0("Table2-Gamma-EL-vs-T-m-", m, ".eps"))
    save_ffr_power_plot(results = result, m = m, alpha = alpha, file_path = file_path)
    message("图片已保存：", file_path)
  }

  if (interactive()) {
    draw_ffr_power(results = result, m = m, alpha = alpha)
  }

  result
}

#==============================================#
#                 参数设置                     #
#==============================================#
nsim <- 2000
ns <- c(1000)
ms <- c(10, 30)
pi1s <- c(0, 0.5, 1)               # 备择组占比 pi1：0=全部原假设；0.5=一半为备择组
effect_sizes <- seq(0, 1, 0.05)    # 备择组效应量：disparity 超出 delta0 的量，0 到 1 步长 0.05
shape <- 2                         # shape 越小，Gamma 分布偏态越强
delta0 <- 1.25                     # H0j: Delta_j <= 1.25
alpha <- 0.05
methods <- c("EL", "T")            # ELBH 为主方法，TBH 为偏态数据下的参照方法
save_data <- T                  # TRUE: 保存 CSV；FALSE: 只计算和打印
save_plot <- T                  # TRUE: 额外保存 EPS
plot_dir <- "."                    # CSV 与 EPS 保存目录（默认与脚本同路径）

# 逐个 m 模拟：run_simulation 内部算完一种 m 就保存数据并作图
seed <- 123
set.seed(seed)
for (m in ms) {
  run_simulation(
    ns = ns, m = m, pi1s = pi1s, effect_sizes = effect_sizes,
    nsim = nsim, shape = shape, delta0 = delta0, alpha = alpha,
    methods = methods, save_data = save_data, save_plot = save_plot,
    plot_dir = plot_dir
  )
}
