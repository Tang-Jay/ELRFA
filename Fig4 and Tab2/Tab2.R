#==============================================#
#       读取模拟数据并生成表格与图形         #
#==============================================#
rm(list = ls())

#==============================================#
#                 数据读取                     #
#==============================================#
read_simulation_data <- function(data_dir, m, dgp = "Gamma") {
  file_path <- file.path(data_dir, paste0("Table2-", dgp, "-FFR-Power-m-", m, ".csv"))
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
#                 结果制表                     #
#==============================================#
# 每个 pi0（= 1 - target_pi1）一张表：
#   列 = delta（effect_size），行 = 方法 × 指标（ffr / power）
build_ffr_power_tables <- function(results, m, methods = c("ELBH", "TBH"), delta = NULL) {
  plot_data <- results[results$m == m, , drop = FALSE]
  if (nrow(plot_data) == 0L) {
    stop("没有 m = ", m, " 的模拟数据。", call. = FALSE)
  }
  methods_use <- intersect(methods, unique(plot_data$method))
  if (length(methods_use) == 0L) {
    stop("数据中不存在所需方法：", paste(methods, collapse = ", "), call. = FALSE)
  }
  method_label <- c(ELBH = "EL", TBH = "T")

  all_delta <- sort(unique(plot_data$effect_size))
  if (is.null(delta)) {
    delta_values <- all_delta
  } else {
    # 校验指定 delta 均存在于数据中（容差 1e-8 规避浮点误差），并按指定顺序保留
    matched <- vapply(delta, function(d) any(abs(all_delta - d) < 1e-8), logical(1))
    if (any(!matched)) {
      stop("数据中找不到 delta = ",
           paste(delta[!matched], collapse = ", "), call. = FALSE)
    }
    delta_values <- delta
  }
  # pi0 升序（0 -> 0.5 -> 1）：即 target_pi1 降序
  pi1_values <- sort(unique(plot_data$target_pi1), decreasing = TRUE)

  tables <- list()
  for (pi1 in pi1_values) {
    pi0 <- 1 - pi1
    panel <- plot_data[plot_data$target_pi1 == pi1, , drop = FALSE]
    # pi1 = 0（全部原假设）时没有备择组，power 为 NA，只输出 ffr
    metrics <- if (all(panel$m1 == 0L)) "ffr" else c("ffr", "power")

    row_specs <- expand.grid(
      method = methods_use, metric = metrics,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    # 行顺序：按指标（ffr -> power）再按方法（EL -> T）
    row_specs <- row_specs[order(match(row_specs$metric, metrics),
                                 match(row_specs$method, methods_use)), , drop = FALSE]

    mat <- matrix(NA_real_, nrow = nrow(row_specs), ncol = length(delta_values))
    for (r in seq_len(nrow(row_specs))) {
      for (c in seq_along(delta_values)) {
        hit <- panel[panel$method == row_specs$method[r] &
                       abs(panel$effect_size - delta_values[c]) < 1e-8, , drop = FALSE]
        if (nrow(hit) == 1L) {
          mat[r, c] <- hit[[row_specs$metric[r]]]
        }
      }
    }

    tbl <- data.frame(
      Method = paste0(method_label[row_specs$method], " ", row_specs$metric),
      mat,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    names(tbl)[-1L] <- sprintf("%.2f", delta_values)
    tables[[sprintf("%.2f", pi0)]] <- tbl
  }
  tables
}

print_ffr_power_tables <- function(tables, m, digits = 4) {
  for (pi0 in names(tables)) {
    cat("\n", strrep("=", 80), "\n", sep = "")
    cat(sprintf("\n pi0 = %s   (m = %d)\n", pi0, m))
    tbl <- tables[[pi0]]
    num_cols <- setdiff(names(tbl), "Method")
    formatted <- tbl
    formatted[num_cols] <- lapply(tbl[num_cols], function(x) sprintf("%.*f", digits, x))
    print(formatted, row.names = FALSE)
  }
  invisible(tables)
}

save_ffr_power_tables <- function(tables, m, table_dir = ".", dgp = "Gamma") {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character(0)
  for (pi0 in names(tables)) {
    file_path <- file.path(table_dir, sprintf("Table2-%s-pi0-%s-m-%d.csv", dgp, pi0, m))
    utils::write.csv(tables[[pi0]], file_path, row.names = FALSE)
    paths <- c(paths, file_path)
    message("表格已保存：", file_path)
  }
  invisible(paths)
}

# 把单个 (m, pi0) 的 ffr/power 表渲染成图形（纯 base 绘图，无需外部包）
draw_ffr_power_table <- function(tbl, m, pi0, digits = 4) {
  delta_names <- names(tbl)[-1L]
  method_names <- tbl[[1L]]
  n_col <- length(delta_names)
  n_row <- nrow(tbl)

  # 数值列格式化为字符串（NA -> "—"）
  num_chr <- vapply(seq_len(n_col), function(j) {
    x <- tbl[[j + 1L]]
    ifelse(is.na(x), "—", sprintf("%.*f", digits, x))
  }, character(n_row))

  # 完整网格：第 1 行表头（Method + delta），其后为数据行
  grid <- rbind(c("Method", delta_names),
                cbind(method_names, num_chr))
  n_row_total <- nrow(grid)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(0.6, 0.6, 2.6, 0.6))

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, n_col + 2.2), ylim = c(0, n_row_total + 0.6))

  metric_txt <- paste(toupper(unique(sub("^[^ ]+ ", "", method_names))), collapse = " & ")
  graphics::title(main = sprintf("%s   (m = %d,  pi0 = %s)", metric_txt, m, pi0),
                  line = 0.6, cex.main = 1.15, font.main = 2)

  ys <- (n_row_total:1L) - 0.5
  for (i in seq_len(n_row_total)) {
    fnt <- if (i == 1L) 2 else 1
    graphics::text(0.05, ys[i], labels = grid[i, 1L], adj = 0, cex = 0.95, font = fnt)
    for (k in seq_len(n_col)) {
      graphics::text(k + 1.5, ys[i], labels = grid[i, k + 1L], adj = 1, cex = 0.95, font = fnt)
    }
  }

  # 表头下划线
  graphics::lines(c(0, n_col + 2.2), rep(ys[1L] - 0.55, 2), lwd = 1.5)

  invisible(NULL)
}

# 每个 (m, pi0) 一张 EPS 表图
save_ffr_power_table_plots <- function(tables, m, plot_dir = ".", dgp = "Gamma", digits = 4) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- character(0)
  for (pi0 in names(tables)) {
    tbl <- tables[[pi0]]
    n_col <- ncol(tbl) - 1L
    n_row <- nrow(tbl)
    file_path <- file.path(plot_dir,
                           sprintf("Table2-%s-FFR-Power-m-%d-pi0-%s.eps", dgp, m, pi0))
    width  <- max(7, n_col * 0.6 + 3)
    height <- n_row * 0.42 + 1.8

    open_eps_device(file_path, width = width, height = height)
    ok <- tryCatch(
      { draw_ffr_power_table(tbl, m, pi0, digits = digits); TRUE },
      error = function(e) { message("制表图失败：", conditionMessage(e)); FALSE }
    )
    grDevices::dev.off()
    if (ok) {
      paths <- c(paths, file_path)
      message("表格图片已保存：", file_path)
    }
  }
  invisible(paths)
}

#==============================================#
#                 参数与运行                   #
#==============================================#
process_results <- function(m, results, alpha = 0.05, plot_dir = ".", save_plot = FALSE,
                            save_tables = FALSE, save_table_plot = FALSE, table_dir = ".",
                            dgp = "Gamma", digits = 4, delta = NULL) {
  # ELDI_SMOKE_TEST=true Rscript Table2.R：读取和制表，但不写入 EPS
  smoke_test <- tolower(Sys.getenv("ELDI_SMOKE_TEST", unset = "false")) == "true"
  if (smoke_test) {
    save_plot <- FALSE
  }

  if (nrow(results) == 0L) {
    stop("没有 m = ", m, " 的模拟数据。", call. = FALSE)
  }

  tables <- build_ffr_power_tables(results, m, delta = delta)
  print_ffr_power_tables(tables, m, digits = digits)
  if (save_tables) {
    save_ffr_power_tables(tables, m, table_dir = table_dir, dgp = dgp)
  }
  if (save_table_plot) {
    save_ffr_power_table_plots(tables, m, plot_dir = plot_dir, dgp = dgp, digits = digits)
  }

  file_path <- file.path(plot_dir, paste0("Table2-", dgp, "-EL-vs-T-m-", m, ".eps"))
  if (save_plot) {
    save_ffr_power_plot(results = results, m = m, alpha = alpha, file_path = file_path)
    message("图片已保存：", file_path)
  }

  if (interactive()) {
    draw_ffr_power(results = results, m = m, alpha = alpha)
  }

  invisible(list(table = results, tables = tables, plot_path = file_path))
}

#==============================================#
#                 参数设置                     #
#==============================================#
ms <- c(10, 30)
alpha <- 0.05
data_dir <- "./data"
plot_dir <- "."
table_dir <- "."
dgp <- "Gamma"
delta <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
save_plot <- F
save_tables <- F
save_table_plot <- F

for (m in ms) {
  m_results <- read_simulation_data(data_dir, m, dgp = dgp)
  process_results(
    m = m, results = m_results, alpha = alpha, plot_dir = plot_dir,
    save_plot = save_plot, save_tables = save_tables, save_table_plot = save_table_plot,
    table_dir = table_dir, dgp = dgp, delta = delta
  )
}
