#==============================================#
#          Wald empirical rejection rates      #
#==============================================#
# 与 RejectRatio-EL-Data.R 使用相同的数据生成、同样本 OLS 和场景种子。
# 已知 theta_P = 1，检验 Delta = 1；只累计拒绝次数，不保存逐次数据。

# ---- 与 EL 一致的数据与残差损失 ----
generate_wald_data <- function(n, m, gamma) {
  beta0 <- 2
  x <- runif(n, min = 0, max = 1)
  y <- rnorm(n, mean = beta0 * x, sd = sqrt((1 - gamma) + gamma * x))
  X <- matrix(x, ncol = 1)
  Y <- matrix(y, ncol = 1)
  xtx <- t(X) %*% X
  xtx_inv <- solve(xtx)
  xty <- t(X) %*% y
  beta_hat <- as.numeric(xtx_inv %*% xty)
  y_hat <- beta_hat * X
  loss <- as.numeric((Y - y_hat)^2)
  group <- findInterval(X, seq(0, 1, length.out = m + 1L),
                        rightmost.closed = TRUE)
  list(loss = loss, group = group)
}

# ---- Wald = sum_j n_j * (mean_j - 1)^2 / s_j^2 ----
compute_wald <- function(loss, group, m) {
  groups <- split(loss, factor(group, levels = seq_len(m)))
  nj <- lengths(groups)
  if (any(nj < 2L)) return(NA_real_)
  mu_hat <- vapply(groups, mean, numeric(1))
  s2 <- vapply(groups, var, numeric(1))  # 分母 n_j - 1
  if (any(!is.finite(mu_hat)) || any(!is.finite(s2)) || any(s2 <= 0)) {
    return(NA_real_)
  }
  sum(nj * (mu_hat - 1)^2 / s2)
}

#==============================================#
#                  场景模拟                    #
#==============================================#
run_wald_simulation <- function(ns, ms, gamma_cases, nsim, alpha,
                                seed, save_data, output_dir) {
  is_count_vector <- function(x) {
    is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
      all(x >= 1 & x == floor(x)) && !anyDuplicated(x)
  }
  if (!is_count_vector(ns) || !is_count_vector(ms) ||
      min(ns) < 2 || max(ms) > min(ns)) {
    stop("ns and ms must contain distinct positive integers, with 2 <= n and m <= n.")
  }
  if (!is_count_vector(nsim) || length(nsim) != 1L) {
    stop("nsim must be a positive integer.")
  }
  if (length(gamma_cases) == 0L || anyNA(gamma_cases) ||
      anyDuplicated(gamma_cases) ||
      !all(gamma_cases %in% c("null", "local", "alter"))) {
    stop("gamma_cases must contain distinct values from null, local, alter.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be in (0, 1).")
  }
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("seed must be a nonnegative integer supported by set.seed().")
  }

  grid <- expand.grid(n = ns, m = ms, gamma_type = gamma_cases,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  results <- vector("list", nrow(grid))
  cat("gamma n m alpha Wald\n")
  for (s in seq_len(nrow(grid))) {
    n <- grid$n[s]
    m <- grid$m[s]
    gamma <- switch(grid$gamma_type[s], null = 0, local = 1 / sqrt(n), alter = 1)
    critical_value <- qchisq(1 - alpha, df = m)
    n_reject <- 0L

    # 与 EL 对齐：每个场景重新设置种子，场景内连续生成 nsim 份数据。
    set.seed(seed)
    for (b in seq_len(nsim)) {
      sample <- generate_wald_data(n, m, gamma)
      wald <- compute_wald(sample$loss, sample$group, m)
      # 异常时明确停止，避免静默丢弃重复导致拒绝率分母改变。
      if (!is.finite(wald)) {
        stop(sprintf("Invalid Wald statistic: gamma=%g, n=%d, m=%d, replicate=%d",
                     gamma, n, m, b))
      }
      if (wald > critical_value) n_reject <- n_reject + 1L
    }

    rejection_rate <- n_reject / nsim
    results[[s]] <- data.frame(gamma = gamma, n = n, m = m,
                               alpha = alpha, Wald = rejection_rate)
    cat(sprintf("%.8g %d %d %.4g %.4f\n", gamma, n, m, alpha, rejection_rate))
  }

  results <- do.call(rbind, results)
  if (save_data) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Cannot create output_dir: ", output_dir)
    write.csv(results, file.path(output_dir, "Wald-summary.csv"), row.names = FALSE)
  }
  invisible(results)
}

#==============================================#
#                 参数与执行                   #
#==============================================#
ns <- c(1000, 2000, 4000, 6000)
ms <- c(2, 5, 8)
gamma_cases <- c("null", "local", "alter")
nsim <- 2000
alpha <- 0.05
seed <- 3
save_data <- TRUE  # TRUE: 只保存一张 27 行、5 列的汇总表；FALSE: 只输出并返回结果
output_dir <- "/Users/jtan1995/zotero-obsidian/02 任务/REL/2 手稿/代码/Tab1 and Fig1/data"

# 与 EL 使用相同的 R 版本和 RNGkind；本脚本不强行修改 RNGkind。
wald_results <- run_wald_simulation(
  ns = ns, ms = ms, gamma_cases = gamma_cases, nsim = nsim,
  alpha = alpha, seed = seed, save_data = save_data, output_dir = output_dir
)
