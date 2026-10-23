# =============================================== #
#       Compas Fig5 Audit Runner (EL/CEL)          #
# =============================================== #
# 公平性审计：统一调用 EL/CEL，小样本或数值失败时返回 NA 而非中断
# test：test1/greater、test2/less、test3/or（区间原假设）；CEL 支持三种 alternative
source("Functions/EL_CI.R")  # ELCI 暂不计算
source("Functions/CEL_CI.R")

normalize_audit_test <- function(test) {
  key <- tolower(trimws(as.character(test)))
  switch(
    key,
    "test1" = "test1",
    "1" = "test1",
    "greater" = "test1",
    "大于" = "test1",
    "test2" = "test2",
    "2" = "test2",
    "less" = "test2",
    "小于" = "test2",
    "test3" = "test3",
    "3" = "test3",
    "or" = "test3",
    "two.sided" = "test3",
    "双边" = "test3",
    stop(
      "test 须为 test1 / test2 / test3（greater / less / or），当前为: ",
      test
    )
  )
}

audit_test_alternative <- function(test) {
  switch(
    normalize_audit_test(test),
    test1 = "greater",
    test2 = "less",
    test3 = "two.sided"
  )
}

ci_al_name <- function(ci_method) {
  switch(ci_method,
    "ELCI" = "ELAL",
    "CELCI" = "CELAL",
    stop("未知 ci_method: ", ci_method)
  )
}

# 读 CSV：Inf/-Inf 字符串转 numeric（与 Fig6-Data.R fmt_audit_csv 对应）
parse_audit_csv_num <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  x <- trimws(as.character(x))
  out <- rep(NA_real_, length(x))
  blank <- is.na(x) | x == ""
  out[!blank & x == "Inf"] <- Inf
  out[!blank & x == "-Inf"] <- -Inf
  rest <- !blank & x != "Inf" & x != "-Inf"
  if (any(rest)) {
    out[rest] <- suppressWarnings(as.numeric(x[rest]))
  }
  out
}

read_audit_csv <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    dplyr::mutate(
      lb = parse_audit_csv_num(lb),
      ub = parse_audit_csv_num(ub),
      DeltaG = parse_audit_csv_num(DeltaG),
      ELAL = parse_audit_csv_num(ELAL)
    )
}

ci_na_result <- function(ci_method, n, DeltaG = NA_real_) {
  al_name <- ci_al_name(ci_method)
  out <- list(n = n, DeltaG = DeltaG, lb = NA_real_, ub = NA_real_)
  out[[al_name]] <- NA_real_
  out
}

.warn_skip_ci <- function(ci_method, group, msg) {
  grp <- trimws(as.character(group))
  if (!is.na(grp) && nzchar(grp)) {
    warning(sprintf("%s [%s]: %s", ci_method, grp, msg), call. = FALSE)
  } else {
    warning(sprintf("%s: %s", ci_method, msg), call. = FALSE)
  }
}

run_ci <- function(group_specific_data, ci_method, target = theta_P, alpha = 0.95,
                   test = "test1",
                   group = NULL, plot = TRUE, plot_label = NULL) {
  y <- as.numeric(group_specific_data$Y)
  n <- length(y)
  Delta_hat <- if (n > 0L && is.finite(target) && target != 0) {
    mean(y) / target
  } else {
    NA_real_
  }
  min_n <- 2L

  if (n < min_n) {
    .warn_skip_ci(
      ci_method, group,
      sprintf("交叉组 n=%d < %d，跳过置信区间。", n, min_n)
    )
    return(ci_na_result(ci_method, n, Delta_hat))
  }

  x <- y - target
  if (!is.finite(sd(x)) || sd(x) == 0) {
    y_vals <- paste(sort(unique(y)), collapse = ",")
    .warn_skip_ci(
      ci_method, group,
      sprintf(
        "交叉组 n=%d 时 disparity 无变异（Y=%s），跳过。",
        n, y_vals
      )
    )
    return(ci_na_result(ci_method, n, Delta_hat))
  }

  alt <- audit_test_alternative(test)

  plot_title <- if (!is.null(plot_label) && nzchar(trimws(plot_label))) {
    trimws(plot_label)
  } else {
    group
  }

  tryCatch(
    switch(ci_method,
      "ELCI" = EL_CI(group_specific_data, target, alpha,
                     plot = plot, plot_label = plot_title),
      "CELCI" = CEL_CI(
        group_specific_data, target, alpha,
        alternative = alt,
        plot = plot,
        plot_label = plot_title
      ),
      stop("未知 ci_method: ", ci_method)
    ),
    error = function(e) {
      .warn_skip_ci(
        ci_method, group,
        sprintf("交叉组 n=%d 计算失败（%s），记为 NA。", n, conditionMessage(e))
      )
      ci_na_result(ci_method, n, Delta_hat)
    }
  )
}
