#==============================================#
#           Compare EL and T-test Power        #  
#==============================================#
rm(list = ls()) 
# 读取函数
read_power_from_csv <- function(n = 100,
                                errorType,
                                testTypes = c("T04", "T14", "T24", "T34")){
  
  # 初始化数据结构
  taus <- list()
  EL_list <- list()
  Ttest_list <- list()
  
  # 循环读取每个检验类型的CSV文件
  for(testType in testTypes){
    filename <- file.path("data", paste0("PowerCompare-", testType, "-", errorType, "-", n, ".csv"))
    # filename <- paste0("PowerCompare-", testType, "-", errorType, "-", n, ".csv")
    
    if(file.exists(filename)){
      # 读取CSV数据（移除row.names参数，read.csv不接受FALSE值）
      power_data <- read.csv(file = filename)
      
      required_columns <- c("tau", "EL", "Ttest")
      if(!all(required_columns %in% names(power_data))){
        stop("数据文件缺少 tau、EL 或 Ttest 列: ", filename)
      }

      taus[[testType]] <- power_data$tau
      EL_list[[testType]] <- power_data$EL
      Ttest_list[[testType]] <- power_data$Ttest
      
      cat("已读取:", filename, "\n")
    } else {
      stop("文件不存在: ", filename)
    }
  }
  
  # 返回与 calculate_power_EL_Ttest 相同的格式
  return(list(
    taus = taus,
    EL = EL_list,
    Ttest = Ttest_list
  ))
}

# 绘图函数
plot_multi_n_compare <- function(ns = c(50, 100),
                                 errorType,
                                 testTypes = c("T04", "T14", "T24", "T34"),
                                 save_eps = FALSE){
  
  # 调取函数读取所有样本量的数据
  power_results <- list()
  for(n in ns){
    power_results[[as.character(n)]] <- read_power_from_csv(
      n = n,
      errorType = errorType,
      testTypes = testTypes
    )
  }
  
  # 对每个检验类型分别绘图
  for(testType in testTypes){
    if(save_eps){
      setEPS()
      postscript(paste0('Compare-Multi-', testType, '-', errorType, '.eps'), 
                 width = 5.5, height = 5)
    }
    
    # 每类检验使用各自的 tau 网格
    taus <- power_results[[as.character(ns[1])]]$taus[[testType]]
    for(n in ns){
      current_taus <- power_results[[as.character(n)]]$taus[[testType]]
      if(!isTRUE(all.equal(taus, current_taus))){
        stop(testType, " 在不同样本量下的 tau 取值不一致")
      }
    }
    tau_l <- min(taus)
    tau_u <- max(taus)
    
    par(oma = c(0, 0, 0, 0), mar = c(2, 3.2, 2, 1))
    
    # 创建空白图
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
    
    # 设置坐标轴
    axis(1, seq(tau_l, tau_u, 0.1), mgp = c(4, 0, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
    axis(2, seq(0, 1, 0.1), mgp = c(4, 0.2, 0), tcl = 0.5, font = 2, lwd = 2, las = 1)
    axis(3, seq(tau_l, tau_u, 0.1), tck = 0.01, labels = FALSE, tcl = 0, lwd = 2)
    axis(4, seq(0, 1, 0.1), labels = FALSE, tcl = 0, lwd = 2)
    
    # 添加标题和标签（根据testType显示对应的假设）
    # 使用 expression 以便 !=, <=, >= 正确渲染为 ≠, ≤, ≥
    hypothesis_map <- list(
      "T04" = expression(paste(H[0], ": ", Delta == Delta[0], " vs ", H[1], ": ", Delta != Delta[0])),
      "T14" = expression(paste(H[0], ": ", Delta >= Delta[1], " vs ", H[1], ": ", Delta < Delta[1])),
      "T24" = expression(paste(H[0], ": ", Delta <= Delta[2], " vs ", H[1], ": ", Delta > Delta[2])),
      "T34" = expression(paste(H[0], ": ", Delta[1] <= Delta, "" <= Delta[2], " vs ", H[1], ": ", Delta, " < ",Delta[1], " or ", Delta, " > ",Delta[2]))
    )
    null_boundaries <- list(T04 = 1, T14 = 0.8, T24 = 1.25, T34 = c(0.8, 1.25))
    hypothesis_text <- if(testType %in% names(hypothesis_map)) hypothesis_map[[testType]] else paste0("H? vs H", testType)
    title(main = hypothesis_text, cex.main = 1.7)
    title(xlab = "true disparity", line = 1, col.lab = 1, font.lab = 2, cex.lab = 1)
    title(ylab = "power", line = 2, col.lab = 1, font.lab = 2, cex.lab = 1)
    
    # 添加参考线
    abline(h = 0.05, lwd = 1.5, col = "gray", lty = 2)
    abline(v = null_boundaries[[testType]], lwd = 1.5, col = "gray", lty = 2)
    
    # 绘制不同样本量的power曲线
    # 颜色配置（根据样本量动态生成）
    color_palette_el <- c("coral","red1", "red3", "red2" )
    color_palette_ttest <- c("dodgerblue", "blue1", "blue3", "blue4")
    
    # 根据样本量数量动态分配颜色
    colors_el <- setNames(color_palette_el[1:length(ns)], as.character(ns))
    colors_ttest <- setNames(color_palette_ttest[1:length(ns)], as.character(ns))
    lty_el <- setNames(rep(1, length(ns)), as.character(ns))      # EL统一实线
    lty_ttest <- setNames(rep(2, length(ns)), as.character(ns))    # T-test统一虚线
    
    # 图例标签
    legend_labels <- c()
    legend_cols <- c()
    legend_lty <- c()
    
    # 绘制EL方法（所有样本量）
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
    
    # 绘制T-test方法（所有样本量）
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
    
    # 添加图例
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
# 全局参数设置
# 多样本量对比绘图（在同一图中显示EL和T-test在不同样本量下的表现）
ns <- c(20, 50, 500)
errorTypes <- c("error1", "error2")
for(errorType in errorTypes){
  plot_multi_n_compare(ns = ns,
                       errorType = errorType,
                       testTypes = c("T04", "T14", "T24", "T34"),
                       save_eps = T)
}

