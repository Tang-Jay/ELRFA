#==============================================#
#          Gamma = 0 的 EL Q-Q 图          #
#==============================================#

# 读取 EL 统计量并绘制 chi-square Q-Q 图
plot_QQ_figure <- function(n, m, data_dir = "data", save_plot = TRUE, output_path = NULL, width = 4, height = 4, res = 300, confidence_level = 0.95, qtime = 100) {
  file_el <- file.path(data_dir, paste0("null-EL-", n, "-", m, ".csv"))
  if (!file.exists(file_el)) stop("Missing null EL data: ", file_el)
  el_data <- read.csv(file_el)
  if (ncol(el_data) != 1L || !is.numeric(el_data[[1L]])) {
    stop("EL file must contain one numeric statistic column: ", file_el)
  }
  el <- el_data[[1L]]
  if (any(!is.finite(el)) || length(el) < qtime) {
    stop("EL must contain at least qtime finite values")
  }

  # 非交互运行且不保存时只返回数据，避免自动创建 Rplots.pdf。
  if (!save_plot && !interactive()) return(invisible(list(EL = el)))

  probabilities <- (seq_len(qtime) - 0.5) / qtime
  theoretical_quantiles <- qchisq(probabilities, df = m)
  critical_value <- qchisq(confidence_level, df = m)

  if (save_plot) {
    if (is.null(output_path)) {
      stop("output_path is required when save_plot = TRUE")
    }
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    extension <- tolower(tools::file_ext(output_path))
    if (extension == "eps") {
      postscript(output_path, width = width, height = height, paper = "special", horizontal = FALSE, onefile = FALSE, family = "Helvetica")
    } else if (extension == "png") {
      png(output_path, width = width, height = height, units = "in", res = res)
    } else {
      stop("output_path must have extension .eps or .png")
    }
  }

  old_par <- par(no.readonly = TRUE)
  device_open <- save_plot
  on.exit({
    par(old_par)
    if (device_open && dev.cur() > 1L) {
      dev.off()
    }
  }, add = TRUE)

  par(mfrow = c(1, 1), mar = c(2.2, 4, 2, 1), oma = c(0, 0, 0, 0))

  empirical_quantiles <- quantile(el, probs = probabilities, names = FALSE, type = 1)
  axis_limit <- max(empirical_quantiles)
  axis_ticks <- seq(0, axis_limit + 5, by = 2)
  identity_line <- seq(0, axis_limit + 5, by = 2)

  plot(identity_line, identity_line, xlim = c(0, axis_limit), xaxs = "i", xaxt = "n", ylim = c(0, axis_limit), yaxs = "i", yaxt = "n", ann = FALSE, type = "l", lwd = 1.5)
  axis(1, at = axis_ticks, mgp = c(4, 0,   0), tcl = 0.5, font = 2, lwd = 1.5, las = 1)
  axis(2, at = axis_ticks, mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 1.5, las = 1)
  axis(3, at = axis_ticks, tck = 0.01, labels = FALSE, tcl = 0, lwd = 1.5)
  axis(4, at = axis_ticks, tck = 0.01, labels = FALSE, tcl = 0, lwd = 1.5)

  title(main = paste0("EL  n = ", n, "  m = ", m))
  title(xlab = "Chi-Square quantile", line = 1, col.lab = 1, font.lab = 2, cex.lab = 1)
  title(ylab = "EL quantile", line = 2, col.lab = 1, font.lab = 2, cex.lab = 1)
  points(theoretical_quantiles, empirical_quantiles, pch = 20, cex = 0.7, col = "blue1")
  points(critical_value, critical_value, pch = 5, cex = 1.1, lwd = 1.3, col = "blue1")

  invisible(list(EL = el))
}

#==============================================#
#                   参数设置                   #
#==============================================#
ns <- c(1000, 2000, 4000, 6000)
ms <- c(2, 5, 8)
data_dir <- "data"
save_plot <- T

for (m in ms) {
  for (n in ns) {
    output_path <- file.path(".", paste0("QQ-EL-", n, "-", m, ".eps"))
    plot_QQ_figure(n = n, m = m, data_dir = data_dir, save_plot = save_plot, output_path = output_path, qtime = 100)
  }
}
