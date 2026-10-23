#==============================================#
#       EL rejection ratios under Gamma DGP     #
#==============================================#
library(emplik)

# ---- 辅助函数：将样本分配到 m 个等宽子群 ----
generate_one_hot_matrix <- function(X, m) {
  breaks <- seq(0, 1, length.out = m + 1)
  group_indices <- findInterval(X, breaks, rightmost.closed = TRUE)
  groups <- matrix(0, nrow = length(X), ncol = m)
  groups[cbind(seq_along(X), group_indices)] <- 1
  groups
}

# 保留现有函数接口；a 为卡方分位水平，显著性水平为 1 - a。
generate_Data <- function(n, m, gamma, a = 0.95, nsim = 2000,
                          save_data = TRUE) {
  is_count <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x == floor(x)
  }
  if (!is_count(n) || !is_count(m) || !is_count(nsim) || n <= m) {
    stop("n, m, and nsim must be positive integers, with n > m.")
  }
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a) || a <= 0 || a >= 1) {
    stop("a must be a finite scalar in (0, 1).")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) ||
      gamma < 0 || gamma > 1) {
    stop("gamma must be a finite scalar in [0, 1].")
  }
  if (!is.logical(save_data) || length(save_data) != 1L || is.na(save_data)) {
    stop("save_data must be TRUE or FALSE.")
  }

  beta0 <- 2
  theta_P <- 1
  # 检验 H0: Delta = 1；gamma 只控制真实值偏离原假设的强度。
  j <- seq_len(m)
  delta_true <- (1 - gamma) + gamma * (2 * j - 1) / (2 * m)
  critical_value <- qchisq(a, df = m)
  null_moments <- numeric(m)
  el_statistics <- numeric(nsim)

  for (i in seq_len(nsim)) {
    # Var(Y | X) = (1 - gamma) + gamma * X。
    x <- runif(n, min = 0, max = 1)
    y <- rnorm(n, mean = beta0 * x, sd = sqrt((1 - gamma) + gamma * x))

    # 无截距最小二乘估计及平方损失，无需显式求逆。
    beta_hat <- as.numeric(crossprod(x, y) / crossprod(x))
    loss <- (y - beta_hat * x)^2

    # 每列对应一个子群的估计方程 (loss - theta_P) * 1_G。
    groups <- generate_one_hot_matrix(x, m)
    moments <- (loss - theta_P) * groups
    el_statistics[i] <- el.test(
      x = moments, mu = null_moments, maxit = 300
    )[["-2LLR"]]
  }

  rejection_ratio <- mean(el_statistics > critical_value)
  cat(gamma, n, m, 1 - a, rejection_ratio, "\n")

  # 保存逐次 EL 统计量，供 Table1.R 读取；输出到当前工作目录 getwd()。
  if (save_data) {
    gamma_type <- if (gamma == 0) "null" else if (gamma < 1) "local" else "alter"
    file_el <- paste0(gamma_type, "-EL-", n, "-", m, ".csv")
    write.csv(data.frame(x = el_statistics), file = file_el, row.names = FALSE)
  }

  invisible(data.frame(
    gamma = gamma, n = n, m = m, alpha = 1 - a,
    EL = rejection_ratio,
    max_delta_deviation = max(abs(delta_true - 1))
  ))
}

#==============================================#
#                   参数设置                   #
#==============================================#
a <- 0.95
nsim <- 2000
save_data <- T
ms <- c(2, 5, 8)
ns <- c(1000)  # 2000, 4000, 6000
gamma_cases <- c("null", "local", "alter")
seed <- 3

cat("gamma n m alpha EL\n")
for (gamma_case in gamma_cases) {
  for (m in ms) {
    for (n in ns) {
      # 保留原有逐场景重置种子的方式，以便核对已有模拟结果。
      set.seed(seed)
      gamma <- switch(gamma_case, null = 0, local = 1 / sqrt(n), alter = 1)
      generate_Data(
        n = n, m = m, gamma = gamma, a = a,
        nsim = nsim, save_data = save_data
      )
    }
  }
}
