#===============================================
#                 Compas P14-value
#        H0: ε_G ≥ ε_0  vs  H1: ε_G < ε_0
#===============================================
P14_value<-function(group_specific_data,target=theta_P,Delta0=0.05){
  # 超参数赋值
  m <- 1
  target <- theta_P
  sex_name <- unique(group_specific_data$sex)
  age_cat_name <- unique(group_specific_data$age_cat)
  
  # print(sex_name)
  if(length(sex_name)!=1){ sex_name = ''}
  if(length(age_cat_name)!=1){ age_cat_name = ''}
  
  # 计算Loss(已知Y_h=1,选出Y=1,则L=0,即L取Y)
  Li <- group_specific_data$Y
  n <- length(Li)
  L <- matrix(Li, n, m)
  Theta <- matrix(theta_P, n, m)
  
  # 估计方程
  # 检验 H0: Delta <= Delta0
   E <-  matrix(Delta0, n, m, byrow = TRUE)
  L1 <- L/Theta
   g <- L  - Theta * E
  
  # 计算EL值
  z = g
  lam = c(lambdaChen(z))
  T04 = 2*sum(log(1+t(lam)%*%t(z)))
  if(mean(L1) < Delta0){T14 <- T04} else {T14 <- 0}
  
  # 调包计算
  # library(emplik)
  # if (mean(L1) <= Delta0) {
  #   # 如果估计值小于等于Delta0，T14 = 0
  #   T14 <- 0
  # } else {
  #   # 检验 H0: Delta = Delta0
  #   el_result <- el.test(L- Theta, mu = Delta0)
  #   T04 <- el_result$`-2LLR`
  #   T14 <- T04
  # }
  # print(T14)
  
  # 计算p值
  if (T14 == 0) {
    p_value <- 1
  } else {
    p_value <- 0.5 * (1 - pchisq(T14, df = 1))
  }
  
  return(p_value)
}

