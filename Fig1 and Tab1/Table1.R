#==============================================#
#         Table 1: Wald、EL 的经验拒绝概率       #
#==============================================#
# Wald 读取五列汇总表；EL 从逐次统计量计算拒绝比例。
table1_merge <- function(
  gamma_cases = c("null", "local", "alter"),
  ns = c(2000, 4000, 6000),
  ms = c(2, 5, 8),
  alpha = 0.05,
  data_dir = "data"
) {
  if (length(gamma_cases) == 0L || anyNA(gamma_cases) ||
      anyDuplicated(gamma_cases) ||
      !all(gamma_cases %in% c("null", "local", "alter"))) {
    stop("gamma_cases must contain distinct values from null, local, alter.")
  }
  is_count_vector <- function(x) {
    is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
      all(x >= 1 & x == floor(x)) && !anyDuplicated(x)
  }
  if (!is_count_vector(ns) || !is_count_vector(ms)) {
    stop("ns and ms must contain distinct positive integers.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a finite scalar in (0, 1).")
  }

  # ---- Wald：按参数匹配，不能依赖 CSV 的行顺序 ----
  file_wald <- file.path(data_dir, "Wald-summary.csv")
  if (!file.exists(file_wald)) stop("Missing Wald summary: ", file_wald)
  wald <- read.csv(file_wald)
  required <- c("gamma", "n", "m", "alpha", "Wald")
  if (!all(required %in% names(wald))) {
    stop("Wald-summary.csv must contain gamma, n, m, alpha, Wald.")
  }
  wald <- wald[, required, drop = FALSE]
  if (nrow(wald) == 0L || !all(vapply(wald, is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(wald))) || any(wald$Wald < 0 | wald$Wald > 1)) {
    stop("Invalid values in Wald-summary.csv.")
  }

  rows <- vector("list", length(gamma_cases) * length(ms) * length(ns))
  row_index <- 1L
  for (gamma_case in gamma_cases) {
    for (m in ms) {
      critical_value <- qchisq(1 - alpha, df = m)
      for (n in ns) {
        gamma <- switch(gamma_case, null = 0, local = 1 / sqrt(n), alter = 1)
        # 容忍 CSV 序列化对 gamma、alpha 引入的浮点舍入误差。
        match_wald <- which(
          wald$n == n & wald$m == m &
          abs(wald$gamma - gamma) < 1e-12 &
          abs(wald$alpha - alpha) < 1e-12
        )
        if (length(match_wald) != 1L) {
          stop("Expected exactly one Wald row for ", gamma_case,
               ", n=", n, ", m=", m, ", alpha=", alpha,
               "; found ", length(match_wald), ".")
        }

        # ---- EL：按组数 m 的卡方临界值计算拒绝比例 ----
        file_el <- file.path(data_dir, paste0(gamma_case, "-EL-", n, "-", m, ".csv"))
        if (!file.exists(file_el)) stop("Missing EL data: ", file_el)
        el_data <- read.csv(file_el)
        if (ncol(el_data) != 1L) {
          stop("EL file must contain one statistic column: ", file_el)
        }
        el <- el_data[[1L]]
        if (!is.numeric(el) || length(el) == 0L || any(!is.finite(el))) {
          stop("Invalid EL data: ", file_el)
        }
        rows[[row_index]] <- data.frame(
          gamma = gamma, n = n, m = m, alpha = alpha,
          Wald = wald$Wald[match_wald],
          EL = mean(el > critical_value)
        )
        row_index <- row_index + 1L
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  print(out, row.names = FALSE, digits = 8)
  invisible(out)
}

#==============================================#
#                   参数设置                   #
#==============================================#
gamma_cases <- c("null", "local", "alter")
ns <- c(1000, 2000, 4000, 6000)
ms <- c(2, 5, 8)
alpha <- 0.05
data_dir <- "/Users/jtan1995/zotero-obsidian/02 任务/REL/2 手稿/代码/Tab1 and Fig1/data"

table1_results <- table1_merge(
  gamma_cases = gamma_cases, ns = ns, ms = ms,
  alpha = alpha, data_dir = data_dir
)
