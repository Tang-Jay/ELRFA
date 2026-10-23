#===============================================
#   p14  african_american VS caucasian  分组   #
#===============================================
rm(list = ls())
source('Functions/Compas-p14-values.R')
source('Functions/GlambdaChen.R')
# BH程序
ELBH_procedure <- function(p_values, alpha = 0.05) {
  
  m <- length(p_values)
  
  sorted_p <- sort(p_values)
  k <- max(which(sorted_p <= (1:m) * alpha / m), 0)
  
  if (k == 0) {
    rejected <- rep(FALSE, m)
  } else {
    threshold <- sorted_p[k]
    rejected <- p_values <= threshold
  }
  
  # 获取被拒绝的p值的索引
  rejected_indices <- which(rejected)
  
  # 输出结果
  # cat('BH procedure results (alpha =', alpha, '):\n')
  # cat('Threshold p-value:', if(k == 0) "None" else sorted_p[k], '\n')
  # cat('Number of rejections:', length(rejected_indices), '\n')
  cat('Rejected p-values at indices:', rejected_indices, '\n')
  # cat('Corresponding p-values:', p_values[rejected_indices], '\n')
  cat('-------------------------------\n')
  
  # 返回拒绝的索引（可选）
  return(rejected_indices)
}

# ---------------------获得数据------------------- #
library(dplyr)
library(readr)
# 读取本地 COMPAS 数据集
compas_data <- read_csv("data/compas_data.csv")
table(compas_data$sex)
table(compas_data$age_cat)
table(compas_data$race)

# 一、数据处理
identical(compas_data$decile_score...12, compas_data$decile_score...40)
identical(compas_data$priors_count...15, compas_data$priors_count...49)
compas_data <- compas_data %>%
  rename(decile_score = decile_score...12) %>%
  select(-decile_score...40)  # 删除重复列
compas_data <- compas_data %>%
  rename(priors_count = priors_count...15) %>%
  select(-priors_count...49)  # 删除重复列 

# 二、数据筛选
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

# 选出目标-女性
Female_data <- clean_data %>%
  filter(sex == "Female")

# ---------------------基准线-------------------- #
# 全种族平均值基准线
postive_data <- clean_data %>%
  filter(Y_h == 1)
table(postive_data$Y)
theta_P = mean(postive_data$Y)  #(0.6135062)

# ---------------------调整参数------------------- #
cat('\ntheta_P =', theta_P)
epsilon0 <- 0.93
alpha <- 0.05
target_data <- Female_data
postive_target_data <- target_data %>%
  filter(Y_h == 1)
sexs <- unique(postive_target_data$sex)          
age_cats <- unique(postive_target_data$age_cat)       
races <- unique(postive_target_data$race) 
# ---------------------开始计算------------------- #
m = 1
p_values <- rep(0,m)
# ALL
# 输出
group_specific_data <- postive_target_data
cat('\n-------------------------------\n')
cat(m,'Female,')
cat("分组样本数:", nrow(group_specific_data), "\n")
# 计算置信区间
p_values[m] <- P14_value(group_specific_data,theta_P,epsilon0)
ELBH_procedure(p_values,alpha)


# Female * age_cat
m = 3
p_values <- rep(0,m)
i = 0
# for(sex_name in sexs){
  for(age_cat_name in age_cats){
    i = i + 1
    # 选择分组数据
    group_specific_data <- postive_target_data %>%
      filter(
        sex == "Female",
        age_cat == age_cat_name
      )
    # 输出
    # cat('\n-------------------------------\n')
    cat(i,age_cat_name,',')
    cat("分组样本数:", nrow(group_specific_data), "\n")
    # 计算置信区间
    p_values[i] <- P14_value(group_specific_data,theta_P,epsilon0)
  }
# }
ELBH_procedure(p_values,alpha)


# Less than 25  * race
races = c("Caucasian","African-American","Hispanic","Other" )
m = length(races)
p_values <- rep(0,m)
i = 0
for(race_name in races){
  i = i + 1
  # 选择分组数据
  group_specific_data <- postive_target_data %>%
    filter(
      sex == "Female",
      age_cat == "Less than 25",
      race == race_name
    )
  p14 = P14_value(group_specific_data,theta_P,epsilon0)
  # 输出
  # cat('\n-------------------------------\n')
  cat(i,"Less than 25",',',race_name,',')
  cat(nrow(group_specific_data),'\n')
  # cat(mean(group_specific_data$Y)-theta_P,',',p14,"\n")
  p_values[i] <- p14
}
ELBH_procedure(p_values,alpha)


# 25 - 45  * race
races = c("Caucasian","African-American","Hispanic","Native American" )
m = length(races)
p_values <- rep(0,m)
i = 0
for(race_name in races){
  i = i + 1
  # 选择分组数据
  group_specific_data <- postive_target_data %>%
    filter(
      sex == "Female",
      age_cat == "25 - 45",
      race == race_name
    )
  p14 = P14_value(group_specific_data,theta_P,epsilon0)
  # 输出
  # cat('\n-------------------------------\n')
  cat(i,"25 - 45",',',race_name,',')
  cat(nrow(group_specific_data),'\n')
  # cat(mean(group_specific_data$Y)-theta_P,',',p14,"\n")
  p_values[i] <- p14
}
ELBH_procedure(p_values,alpha)


# Greater than 45  * race
races = c("Caucasian","African-American","Hispanic","Native American" )
m = length(races)
p_values <- rep(0,m)
i = 0
for(race_name in races){
  i = i + 1
  # 选择分组数据
  group_specific_data <- postive_target_data %>%
    filter(
      sex == "Female",
      age_cat == "Greater than 45",
      race == race_name
    )
  p14 = P14_value(group_specific_data,theta_P,epsilon0)
  # 输出
  # cat('\n-------------------------------\n')
  cat(i,"Greater than 45",',',race_name,',')
  cat(nrow(group_specific_data),'\n')
  # cat(mean(group_specific_data$Y)-theta_P,',',p14,"\n")
  p_values[i] <- p14
}
ELBH_procedure(p_values,alpha)

