#==============================================#
#           Compare EL and T-test Power        #  
#==============================================#
# T04：等式检验；
# T14：单边检验（Lower-备择）
# T24：单边检验（Upper-备择）
# T34：区间检验
#==============================================#
rm(list = ls()) 
source('GlambdaChen.R')
generate_one_hot_matrix <- function(X, m) {
  breaks <- seq(0, 1, length.out = m + 1)
  group_indices <- findInterval(X, breaks, rightmost.closed = TRUE)
  one_hot_matrix <- matrix(0, nrow = length(X), ncol = m)
  one_hot_matrix[cbind(1:length(X), group_indices)] <- 1
  return(one_hot_matrix)
}

# 四类检验分别使用各自的 tau 取值
default_tau_grids <- function() {
  list(
    T04 = seq(0.50, 1.50, 0.10),
    T14 = seq(0.30, 1.30, 0.10),
    T24 = seq(0.75, 1.75, 0.10),
    T34 = seq(0.50, 1.50, 0.05)
  )
}

# 区间假设检验函数（T-test）
interval_test <- function(sample_data, mu1, mu2, alpha = 0.05){
  test_upper <- t.test(sample_data, mu = mu2, alternative = "greater", conf.level = 1 - alpha)
  test_lower <- t.test(sample_data, mu = mu1, alternative = "less", conf.level = 1 - alpha)
  reject_H0 <- (test_upper$p.value < alpha) || (test_lower$p.value < alpha)
  
  return(list(
    reject = reject_H0,
    p_upper = test_upper$p.value,
    p_lower = test_lower$p.value,
    test_upper = test_upper,
    test_lower = test_lower
  ))
}

#==============================================#
#           Power 计算函数                     #
#==============================================#
calculate_power_EL_Ttest <- function(n = 100,
                                     errorType = "error1",
                                     Delta0 = 0.05,
                                     alpha = 0.05,
                                     nsim = 1000,
                                     m = 1,
                                     tau_grids = default_tau_grids(),
                                     save_csv = FALSE){

  test_types <- c("T04", "T14", "T24", "T34")
  if(!identical(sort(names(tau_grids)), sort(test_types))){
    stop("tau_grids must contain T04, T14, T24, and T34")
  }

  # 临界值
  cut1 <- qchisq(1 - alpha, 1)      # 用于双侧检验
  cut2 <- qchisq(1 - 2*alpha, 1)    # 用于单侧检验

  cat("\n 开始计算 power...\n")
  cat("参数: n =", n, ", nsim =", nsim, ", errorType =", errorType, ", m =", m, "\n\n")

  results <- lapply(test_types, function(testType) {
    taus <- tau_grids[[testType]]
    power_el <- numeric(length(taus))
    power_ttest <- numeric(length(taus))

    for(j in seq_along(taus)){
      tau <- taus[j]
      beta0 <- 2
      beta1 <- beta0 - 2*tau
      reject_el <- logical(nsim)
      reject_ttest <- logical(nsim)

      for(i in seq_len(nsim)){
        Xi <- runif(n, min = 0, max = 1)
        if(errorType == "error1"){
          v <- rnorm(n, mean = 0, sd = 1)
        } else if(errorType == "error2"){
          v <- rexp(n, rate = 1) - 1
        } else {
          v <- runif(n, min = -1, max = 1 + tau)
        }

        Yi <- beta0 * Xi + v
        X <- matrix(Xi, ncol = 1)
        Y <- matrix(Yi, ncol = 1)
        Y_h <- beta1 * X
        Li <- Y - Y_h
        S <- generate_one_hot_matrix(X, m)
        L <- matrix(Li, n, m)
        E0 <- matrix(Delta0, n, m, byrow = TRUE)
        E1 <- matrix(Delta1, n, m, byrow = TRUE)
        E2 <- matrix(Delta2, n, m, byrow = TRUE)

        if(testType == "T04"){
          z <- (L - E0) * S
          lam <- lambdaChen(z)
          statistic_el <- 2*sum(log(1 + t(lam) %*% t(z)))
          reject_el[i] <- statistic_el > cut1
          reject_ttest[i] <- t.test(Li, mu = Delta0, alternative = "two.sided")$p.value < alpha
        } else if(testType == "T14"){
          statistic_el <- if(mean(L) <= Delta1){
            z <- (L - E1) * S
            lam <- lambdaChen(z)
            2*sum(log(1 + t(lam) %*% t(z)))
          } else {
            0
          }
          reject_el[i] <- statistic_el > cut2
          reject_ttest[i] <- t.test(Li, mu = Delta1, alternative = "less")$p.value < alpha
        } else if(testType == "T24"){
          statistic_el <- if(mean(L) >= Delta2){
            z <- (L - E2) * S
            lam <- lambdaChen(z)
            2*sum(log(1 + t(lam) %*% t(z)))
          } else {
            0
          }
          reject_el[i] <- statistic_el > cut2
          reject_ttest[i] <- t.test(Li, mu = Delta2, alternative = "greater")$p.value < alpha
        } else {
          if(mean(L) <= Delta1){
            z <- (L - E1) * S
            lam <- lambdaChen(z)
            statistic_el <- 2*sum(log(1 + t(lam) %*% t(z)))
          } else if(mean(L) >= Delta2){
            z <- (L - E2) * S
            lam <- lambdaChen(z)
            statistic_el <- 2*sum(log(1 + t(lam) %*% t(z)))
          } else {
            statistic_el <- 0
          }
          reject_el[i] <- statistic_el >= cut2
          reject_ttest[i] <- interval_test(L, Delta1, Delta2, alpha = alpha)$reject
        }
      }

      power_el[j] <- mean(reject_el)
      power_ttest[j] <- mean(reject_ttest)
      if(j %% 5 == 0 || j == length(taus)){
        cat(sprintf("%s, tau = %.2f: EL = %.4f, T-test = %.4f\n",
                    testType, tau, power_el[j], power_ttest[j]))
      }
    }

    list(taus = taus, EL = power_el, Ttest = power_ttest)
  })
  names(results) <- test_types

  # 保存CSV数据（循环保存所有检验类型）
  if(save_csv){
    for(testType in test_types){
      data_df <- data.frame(
        tau = results[[testType]]$taus,
        EL = results[[testType]]$EL,
        Ttest = results[[testType]]$Ttest
      )
      
      filename <- paste0("PowerCompare-", testType, "-", errorType, "-", n, ".csv")
      write.csv(data_df, file = filename, row.names = FALSE)
      cat("已保存:", filename, "\n")
    }
  }


  return(list(
    taus = lapply(results, `[[`, "taus"),
    EL = lapply(results, `[[`, "EL"),
    Ttest = lapply(results, `[[`, "Ttest")
  ))
}

#==============================================#
#           多样本量对比绘图函数               #
# 与 Figs2-3.R 同款图，数据来自内存  #
#==============================================#
plot_multi_n_compare <- function(power_results,
                                 ns,
                                 errorType,
                                 testTypes = c("T04", "T14", "T24", "T34"),
                                 save_eps = FALSE){
  # power_results: list, names = as.character(n)，每个元素为 calculate_power_EL_Ttest 的返回值
  for(testType in testTypes){
    if(save_eps){
      setEPS()
      postscript(paste0('Compare-Multi-', testType, '-', errorType, '.eps'),
                 width = 5.5, height = 5)
    }
    taus <- power_results[[as.character(ns[1])]]$taus[[testType]]
    tau_l <- min(taus)
    tau_u <- max(taus)
    par(oma = c(0, 0, 0, 0), mar = c(2, 3.2, 2, 1))
    plot(taus, taus,
         type = "l",
         xaxs = "i",
         yaxs = "i",
         ylim = c(0, 1),
         xlab = "",
         ylab = "",
         col = "white",
         xaxt = "n",
         yaxt = "n")
    axis(1, seq(tau_l, tau_u, 0.1), mgp = c(4, 0, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
    axis(2, seq(0, 1, 0.1), mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
    axis(3, seq(tau_l, tau_u, 0.1), tck = 0.01, labels = F, tcl = 0, lwd = 2)
    axis(4, seq(0, 1, 0.1), labels = F, tcl = 0, lwd = 2)
    hypothesis_map <- list(
      "T04" = expression(paste(H[0], ": ", Delta == Delta[0], " vs ", H[1], ": ", Delta != Delta[0])),
      "T14" = expression(paste(H[0], ": ", Delta >= Delta[1], " vs ", H[1], ": ", Delta < Delta[1])),
      "T24" = expression(paste(H[0], ": ", Delta <= Delta[2], " vs ", H[1], ": ", Delta > Delta[2])),
      "T34" = expression(paste(H[0], ": ", Delta[1] <= Delta, "" <= Delta[2], " vs ", H[1], ": ", Delta, " < ",Delta[1], " or ", Delta, " > ",Delta[2]))
    )
    hypothesis_text <- if(testType %in% names(hypothesis_map)) hypothesis_map[[testType]] else paste0("H? vs H", testType)
    title(main = hypothesis_text, cex.main = 1.7)
    title(xlab = "true disparity", line = 1, col.lab = 1, font.lab = 2, cex.lab = 1)
    title(ylab = "power", line = 2, col.lab = 1, font.lab = 2, cex.lab = 1)
    abline(h = 0.05, lwd = 1.5, col = "gray", lty = 2)
    abline(v = 0, lwd = 1.5, col = "gray", lty = 2)
    color_palette_el <- c("coral", "red1", "red3", "red2")
    color_palette_ttest <- c("dodgerblue", "blue1", "blue3", "blue4")
    colors_el <- setNames(color_palette_el[1:length(ns)], as.character(ns))
    colors_ttest <- setNames(color_palette_ttest[1:length(ns)], as.character(ns))
    lty_el <- setNames(rep(1, length(ns)), as.character(ns))
    lty_ttest <- setNames(rep(2, length(ns)), as.character(ns))
    legend_labels <- c()
    legend_cols <- c()
    legend_lty <- c()
    for(n in ns){
      n_str <- as.character(n)
      if(!is.null(power_results[[n_str]]$EL[[testType]])){
        lines(taus, power_results[[n_str]]$EL[[testType]],
              col = colors_el[n_str], lty = lty_el[n_str], lwd = 2)
        legend_labels <- c(legend_labels, paste0("EL(n=", n, ")"))
        legend_cols <- c(legend_cols, colors_el[n_str])
        legend_lty <- c(legend_lty, lty_el[n_str])
      }
    }
    for(n in ns){
      n_str <- as.character(n)
      if(!is.null(power_results[[n_str]]$Ttest[[testType]])){
        lines(taus, power_results[[n_str]]$Ttest[[testType]],
              col = colors_ttest[n_str], lty = lty_ttest[n_str], lwd = 2)
        legend_labels <- c(legend_labels, paste0("T-test(n=", n, ")"))
        legend_cols <- c(legend_cols, colors_ttest[n_str])
        legend_lty <- c(legend_lty, lty_ttest[n_str])
      }
    }
    legend(x = tau_l + 0.03, y = 0.3, inset = 0.02, cex = 0.6,
           legend = legend_labels,
           lty = legend_lty,
           col = legend_cols,
           lwd = 2)
    if(save_eps) dev.off()
  }
}

#==============================================#
#           使用示例                           #
#==============================================#
# 全局参数设置（多样本量对比图）
set.seed(123)
nsim <- 2000  # 我们模拟采用 2000，数据结果储存在脚本同路径中
Delta0 <- 1
Delta1 <- 0.8
Delta2 <- 1.25
ns <- c(20, 50, 500)

# 与 Figs2-3 同款：多样本量对比图（数据由本脚本模拟生成，不读 CSV）
for(errorType in c("error1", "error2")){
  power_results <- list()
  for(n in ns){
    power_results[[as.character(n)]] <- calculate_power_EL_Ttest(
      n = n,
      errorType = errorType,
      Delta0 = Delta0,
      nsim = nsim,
      save_csv = T
    )
  }
  plot_multi_n_compare(power_results,
                      ns = ns,
                      errorType = errorType,
                      testTypes = c("T04", "T14", "T24", "T34"),
                      save_eps = FALSE)
}
