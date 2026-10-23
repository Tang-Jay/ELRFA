#==============================================#
#       读取模拟数据并生成表格与图形         #
#==============================================#
rm(list = ls())

#==============================================#
#                 数据读取                     #
#==============================================#
read_simulation_data <- function(data_dir, m) {
  file_path <- file.path(data_dir, paste0("Table2-Gamma-FFR-Power-m-", m, ".csv"))
  if (!file.exists(file_path)) {
    stop("找不到数据文件：", file_path, call. = FALSE)
  }

  required_columns <- c(
    "method", "effect_size", "n", "m1", "m", "realized_pi1",
    "target_pi1", "ffr", "power"
  )
  table_data <- utils::read.csv(file_path, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, names(table_data))
  if (length(missing_columns) > 0L) {
    stop(sprintf(
      "%s 缺少列：%s", basename(file_path), paste(missing_columns, collapse = ", ")
    ), call. = FALSE)
  }
  table_data[required_columns]
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
    # 全部为原假设时只画 FFR；全部为备择时只画 Power；混合情形两者并列
    is_null <- all(pi1_data$m1 == 0L)
    is_alternative <- all(pi1_data$m1 == pi1_data$m)
    metrics <- if (is_null) {
      "ffr"
    } else if (is_alternative) {
      "power"
    } else {
      c("ffr", "power")
    }
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
    y_label <- if (metric == "ffr") "FFR" else "Power"

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

  graphics::mtext(sprintf("ELBH vs TBH (m = %d)", m),
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
#                 参数与运行                   #
#==============================================#
process_results <- function(m, results, alpha = 0.05, plot_dir = ".", save_plot = FALSE) {
  # ELDI_SMOKE_TEST=true Rscript Table2.R：读取和制表，但不写入 EPS
  smoke_test <- tolower(Sys.getenv("ELDI_SMOKE_TEST", unset = "false")) == "true"
  if (smoke_test) {
    save_plot <- FALSE
  }

  if (nrow(results) == 0L) {
    stop("没有 m = ", m, " 的模拟数据。", call. = FALSE)
  }

  message("m = ", m)
  print(results)

  file_path <- file.path(plot_dir, paste0("Table2-Gamma-EL-vs-T-m-", m, ".eps"))
  if (save_plot) {
    save_ffr_power_plot(results = results, m = m, alpha = alpha, file_path = file_path)
    message("图片已保存：", file_path)
  }

  if (interactive()) {
    draw_ffr_power(results = results, m = m, alpha = alpha)
  }

  invisible(list(table = results, plot_path = file_path))
}

#==============================================#
#                 参数设置                     #
#==============================================#
ms <- c(10, 30)
alpha <- 0.05
data_dir <- "./data"
plot_dir <- "."
save_plot = T

for (m in ms) {
  m_results <- read_simulation_data(data_dir, m)
  process_results(m = m, results = m_results, alpha = alpha, plot_dir = plot_dir, save_plot = save_plot)
}
