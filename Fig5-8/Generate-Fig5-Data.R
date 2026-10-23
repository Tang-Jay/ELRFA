# =============================================== #
#             Compas Fig5 Audit Data (CEL)        #
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

# =====================辅助函数=================== #
age_cat_label <- function(x) {
  switch(x,
    "Less than 25" = "< 25",
    "Greater than 45" = "> 45",
    x
  )
}

sex_label <- function(x) {
  switch(x, "Male" = "M", "Female" = "F", x)
}

# =====================参数设置=================== #
test <- "test1"                        
ci_methods <- c("ELCI", "CELCI")
alpha <- 0.9
dataname <- "African-American"
save_data <- T  # FALSE 时不写 CSV
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
  select(sex,age_cat,race,decile_score,two_year_recid
  ) %>%
  filter(
    !is.na(decile_score), 
    !is.na(two_year_recid)
  ) %>%
  rename(Y = two_year_recid)  %>%
  mutate(Y_h = ifelse(decile_score >= 5, 1, 0)) %>% 
  select(-decile_score)

# 生成 L 列（Y == Y_h 时 L=0，否则 L=1）
clean_data <- clean_data %>%
  mutate(L = ifelse(Y == Y_h, 0, 1))
table(clean_data$L)

# 选出目标族-非裔族
african_american_data <- clean_data %>%
  filter(race == "African-American")
# 选出阳预测非裔族
postive_african_american_data <- african_american_data %>%
  filter(Y_h == 1)
table(postive_african_american_data$L)

# =====================基准组==================== #
# 选出白人族
caucasian_data <- clean_data %>%
  filter(race == "Caucasian")
# 选出阳预测非裔族
postive_caucasian_data <- caucasian_data %>%
  filter(Y_h == 1)
# 查看数据
table(postive_caucasian_data$Y)
table(postive_caucasian_data$L)
# 基准线
theta_P = mean(postive_caucasian_data$Y)
cat('theta_P =', theta_P,'\n')

# ==================分组置信区间计算================ #
sexs <- unique(postive_african_american_data$sex)
age_cats <- unique(postive_african_american_data$age_cat)

for (ci_method in ci_methods) {
  cat("\n==========", ci_method, "==========\n")

  # ALL
  group_specific_data <- postive_african_american_data
  cat("\n--------------------------------\n")
  cat("ALL,")
  cat("分组样本数:", nrow(group_specific_data), "\n")
  run_ci(group_specific_data, ci_method, theta_P, alpha,
         test = test,
         group = "All")

  # Sex
  for (sex_name in sexs) {
    group_specific_data <- postive_african_american_data %>%
      filter(sex == sex_name)
    cat("\n--------------------------------\n")
    cat(sex_name, ",")
    cat("分组样本数:", nrow(group_specific_data), "\n")
    run_ci(group_specific_data, ci_method, theta_P, alpha,
           test = test,
           group = sex_label(sex_name))
  }

  # age_cat
  for (age_cat_name in age_cats) {
    group_specific_data <- postive_african_american_data %>%
      filter(age_cat == age_cat_name)
    cat("\n--------------------------------\n")
    cat(age_cat_name, ",")
    cat("分组样本数:", nrow(group_specific_data), "\n")
    run_ci(group_specific_data, ci_method, theta_P, alpha,
           test = test,
           group = age_cat_label(age_cat_name))
  }

  # sex * age_cat
  for (sex_name in sexs) {
    for (age_cat_name in age_cats) {
      group_specific_data <- postive_african_american_data %>%
        filter(sex == sex_name, age_cat == age_cat_name)
      cat("\n--------------------------------\n")
      cat(sex_name, ",", age_cat_name, ",")
      cat("分组样本数:", nrow(group_specific_data), "\n")
      run_ci(group_specific_data, ci_method, theta_P, alpha,
             test = test,
             group = paste0(sex_label(sex_name), " (", age_cat_label(age_cat_name), ")"))
    }
  }
}

# ===================结果整理与保存================= #
# save_data = TRUE 时写入 save_data_dir/*.csv
audit_rows <- list(
  list("All", postive_african_american_data),
  list("< 25", filter(postive_african_american_data, age_cat == "Less than 25")),
  list("25 - 45", filter(postive_african_american_data, age_cat == "25 - 45")),
  list("> 45", filter(postive_african_american_data, age_cat == "Greater than 45")),
  list("M", filter(postive_african_american_data, sex == "Male")),
  list("F", filter(postive_african_american_data, sex == "Female")),
  list("M (< 25)", filter(postive_african_american_data, sex == "Male", age_cat == "Less than 25")),
  list("F (< 25)", filter(postive_african_american_data, sex == "Female", age_cat == "Less than 25")),
  list("M (25 - 45)", filter(postive_african_american_data, sex == "Male", age_cat == "25 - 45")),
  list("F (25 - 45)", filter(postive_african_american_data, sex == "Female", age_cat == "25 - 45")),
  list("M (> 45)", filter(postive_african_american_data, sex == "Male", age_cat == "Greater than 45")),
  list("F (> 45)", filter(postive_african_american_data, sex == "Female", age_cat == "Greater than 45"))
)
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
