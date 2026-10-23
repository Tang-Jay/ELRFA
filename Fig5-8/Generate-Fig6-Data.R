# =============================================== #
#             Compas Fig7 Audit Data (CEL)        #
# =============================================== #
# test1 / greater: H0: Delta_G <= 1.25        vs H1: Delta_G > 1.25
# test2 / less:    H0: Delta_G >= 1.25        vs H1: Delta_G < 1.25
# test3 / or:      H0: 0.8 <= Delta_G <= 1.25 vs H1: Delta_G < 0.8 or Delta_G > 1.25
# 当前代码根据 test 生成对应方向的 CEL 置信区间，不直接输出检验拒绝结论。

# ====================初始化与依赖================== #
rm(list = ls())
source("Functions/audit_run_ci.R")
library(dplyr)
library(readr)

# =====================参数设置=================== #
test <- "test2"
ci_methods <- c("ELCI", "CELCI")
alpha <- 0.9
dataname <- "Female-Race-Age"
save_data <- TRUE  # FALSE 时不写 CSV
save_data_dir <- "data"
test <- normalize_audit_test(test)

# ===================数据读取与处理================= #
compas_data <- read_csv("data/compas_data.csv")
table(compas_data$sex)
table(compas_data$age_cat)
table(compas_data$race)


# 一、删除重复列
identical(compas_data$decile_score...12, compas_data$decile_score...40)
identical(compas_data$priors_count...15, compas_data$priors_count...49)
compas_data <- compas_data %>%
  rename(decile_score = decile_score...12) %>%
  select(-decile_score...40)  # 删除重复列
compas_data <- compas_data %>%
  rename(priors_count = priors_count...15) %>%
  select(-priors_count...49)  # 删除重复列

# 二、数据筛选与变量构造
clean_data <- compas_data %>%
  select(sex, age_cat, race, decile_score, two_year_recid) %>%
  filter(
    !is.na(decile_score),
    !is.na(two_year_recid)
  ) %>%
  rename(Y = two_year_recid) %>%
  mutate(Y_h = ifelse(decile_score >= 5, 1, 0)) %>%
  select(-decile_score)

# 生成 L 列（Y == Y_h 时 L=0，否则 L=1）
clean_data <- clean_data %>%
  mutate(L = ifelse(Y == Y_h, 0, 1))
table(clean_data$L)

# 选出目标组：女性且预测为阳性
postive_target_data <- clean_data %>%
  filter(sex == "Female", Y_h == 1)
table(postive_target_data$L)

# =====================基准组==================== #
# 全体阳性预测样本的平均值
postive_data <- clean_data %>%
  filter(Y_h == 1)
table(postive_data$Y)
theta_P <- mean(postive_data$Y)
cat('theta_P =', theta_P, '\n')

# ==================分组置信区间计算================ #
audit_rows <- list(
  list("F", postive_target_data),
  list("< 25", filter(postive_target_data, age_cat == "Less than 25")),
  list("25 - 45", filter(postive_target_data, age_cat == "25 - 45")),
  list("> 45", filter(postive_target_data, age_cat == "Greater than 45")),
  # list("Caucasian", filter(postive_target_data, race == "Caucasian")),
  # list("African-American", filter(postive_target_data, race == "African-American")),
  list("< 25, Caucasian", filter(postive_target_data, age_cat == "Less than 25", race == "Caucasian")),
  list("< 25, African-American", filter(postive_target_data, age_cat == "Less than 25", race == "African-American")),
  list("25 - 45, Caucasian", filter(postive_target_data, age_cat == "25 - 45", race == "Caucasian")),
  list("25 - 45, African-American", filter(postive_target_data, age_cat == "25 - 45", race == "African-American")),
  list("> 45, Caucasian", filter(postive_target_data, age_cat == "Greater than 45", race == "Caucasian")),
  list("> 45, African-American", filter(postive_target_data, age_cat == "Greater than 45", race == "African-American"))
)

for (ci_method in ci_methods) {
  cat("\n==========", ci_method, "==========\n")

  for (row in audit_rows) {
    group_label <- row[[1]]
    group_specific_data <- row[[2]]
    cat("\n-------------------------------\n")
    cat("[", group_label, "] 分组样本数:", nrow(group_specific_data), "\n", sep = "")
    run_ci(group_specific_data, ci_method, theta_P, alpha,
           test = test,
           group = group_label)
  }
}

# ===================结果整理与保存================= #
# save_data = TRUE 时写入 save_data_dir/*.csv
if (isTRUE(save_data)) {
  if (!dir.exists(save_data_dir) &&
      !dir.create(save_data_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("无法创建数据保存目录: ", save_data_dir)
  }
  for (ci_method in ci_methods) {
    al_name <- ci_al_name(ci_method)
    audit_df <- bind_rows(lapply(audit_rows, function(row) {
      r <- run_ci(row[[2]], ci_method, theta_P, alpha,
                  test = test,
                  group = row[[1]], plot = FALSE)
      tibble(
        group = row[[1]],
        lb = round(r$lb, 4),
        ub = round(r$ub, 4),
        DeltaG = round(r$DeltaG, 4),
        ELAL = round(r[[al_name]], 4),
        n = r$n
      )
    }))
    out_name <- if (ci_method == "ELCI") {
      sprintf("%s-Fairness-Audit-%.2f-%s.csv", dataname, alpha, ci_method)
    } else {
      sprintf(
        "%s-Fairness-Audit-%.2f-%s-%s.csv",
        dataname, alpha, test, ci_method
      )
    }
    out_path <- file.path(save_data_dir, out_name)
    write_csv(audit_df, out_path)
    cat("已保存:", out_path, "\n")
  }
} else {
  cat("\n未保存 CSV（save_data = FALSE）\n")
}
