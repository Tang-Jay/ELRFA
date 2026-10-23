# =============================================== #
#                 Compas CI Figure                #
#             （Fig5-7-Data.R 输出）              #
# =============================================== #
rm(list = ls())
library(ggplot2)
library(dplyr)
library(readr)
library(Cairo)

source("Functions/audit_run_ci.R")

alpha <- 0.9     # 0.9 / 0.95
test <- "test1"  # test1 / test2
datanames <- switch(test,
  test1 = c("African-American", "Sex-Age"), # test1 = c("African-American", "Sex-Age"),
  test2 = "Female-Race-Age",
  stop("test 须为 test1 / test2，当前为: ", test)
)
ci_methods <- c("ELCI", "CELCI")
save_plot <- T  # FALSE 时仅 print，不写 eps

# ELCI 无 test 参数，输出文件名不含 test；CELCI 含 test
audit_csv_path <- function(dataname, alpha, test = "test1", ci_method = "") {
  if (!nzchar(ci_method)) {
    stop("须指定 ci_method（ELCI / CELCI）。")
  }
  if (ci_method == "ELCI") {
    sprintf("data/%s-Fairness-Audit-%.2f-%s.csv", dataname, alpha, ci_method)
  } else {
    sprintf(
      "data/%s-Fairness-Audit-%.2f-%s-%s.csv",
      dataname, alpha, test, ci_method
    )
  }
}

audit_plot_path <- function(dataname, alpha, test, ci_method) {
  base <- paste0("PPV-", dataname, "-Disparity-", alpha)
  if (ci_method == "ELCI") {
    paste0(base, "-", ci_method, ".eps")
  } else {
    paste0(base, "-", test, "-", ci_method, ".eps")
  }
}

# Female-Race-Age 按年龄组固定配色；其余情形沿用 Fig4 配色
audit_line_colors <- function(groups, dataname) {
  groups <- as.character(groups)
  if (dataname == "Female-Race-Age") {
    return(case_when(
      grepl("^< 25", groups) ~ "blue",
      grepl("^25 - 45", groups) ~ "green3",
      grepl("^> 45", groups) ~ "red",
      TRUE ~ "black"
    ))
  }

  base4 <- c("black", "blue", "green3", "red")
  fig4_colors <- c(base4, rep(base4, each = 2L))
  rep(fig4_colors, length.out = length(groups))
}

audit_line_types <- function(n) {
  ifelse(seq_len(n) <= 4L, "solid",
         ifelse(seq_len(n) %% 2L == 1L, "dotted", "dashed"))
}

elci_spread_groups <- function(dataname) {
  switch(dataname,
    "Sex-Age" = c("All", "M"),
    "African-American" = c("All", "M"),
    "Female-Race-Age" = "F",
    character(0)
  )
}

plot_height_in <- function(n) {
  if (n <= 12L) 6 else 6 + 0.35 * (n - 12L)
}

prep_audit_plot_data <- function(audit_data, ci_method, dataname, x_lim, test = "test1",
                                 group_levels = NULL) {
  normalize_audit_test(test)
  if (is.null(group_levels)) {
    group_levels <- unique(audit_data$group)
  }
  audit_data <- tibble::tibble(group = group_levels) %>%
    dplyr::left_join(audit_data, by = "group")
  n_grp <- length(group_levels)
  spread_grp <- elci_spread_groups(dataname)
  inf_pad <- max(0.008, 0.015 * diff(x_lim))
  lb_inf_plot <- x_lim[1] + inf_pad
  ub_inf_plot <- x_lim[2] - inf_pad
  out <- audit_data %>%
    mutate(
      line_color = audit_line_colors(group, dataname),
      line_type = audit_line_types(n_grp),
      lb_label = case_when(
        is.na(lb) ~ NA_character_,
        is.infinite(lb) & lb < 0 ~ "-Inf",
        is.finite(lb) ~ sprintf("%.4f", lb),
        TRUE ~ NA_character_
      ),
      ub_label = case_when(
        is.na(ub) ~ NA_character_,
        is.infinite(ub) & ub > 0 ~ "Inf",
        is.finite(ub) ~ sprintf("%.4f", ub),
        TRUE ~ NA_character_
      ),
      lb_plot = case_when(
        is.infinite(lb) & lb < 0 ~ lb_inf_plot,
        is.finite(lb) ~ pmin(pmax(lb, lb_inf_plot), ub_inf_plot),
        TRUE ~ NA_real_
      ),
      ub_plot = case_when(
        is.infinite(ub) & ub > 0 ~ ub_inf_plot,
        is.finite(ub) ~ pmin(pmax(ub, lb_inf_plot), ub_inf_plot),
        TRUE ~ NA_real_
      ),
      seg_lb = lb_plot,
      seg_ub = ub_plot,
      elci_spread_label = ci_method == "ELCI" & group %in% spread_grp &
        is.finite(lb) & is.finite(ub),
      lb_label_x = if_else(elci_spread_label, lb_plot - 0.00005, lb_plot),
      ub_label_x = if_else(elci_spread_label, ub_plot + 0.00005, ub_plot),
      lb_label_hjust = if_else(elci_spread_label, 1, 0.5),
      ub_label_hjust = if_else(elci_spread_label, 0, 0.5)
    ) %>%
    mutate(
      lb_label = if_else(is.finite(lb_plot), lb_label, NA_character_),
      ub_label = if_else(is.finite(ub_plot), ub_label, NA_character_)
    )
  out$group <- factor(out$group, levels = group_levels)
  out
}

# x 轴以 DeltaG = 1 为中心，固定为 0.4–1.6
audit_x_axis_from_elci <- function(elci_data) {
  x_lim <- c(0.4, 1.6)
  list(limits = x_lim, breaks = seq(0.4, 1.6, by = 0.2))
}

# 循环 ci_methods × datanames，为每个绘制一张图
for (dataname in datanames) {
  elci_path <- audit_csv_path(dataname, alpha, test, "ELCI")
  elci_ref <- read_audit_csv(elci_path)
  group_levels <- elci_ref$group
  x_axis <- audit_x_axis_from_elci(elci_ref)
  x_lim <- x_axis$limits
  x_breaks <- x_axis$breaks

  for (ci_method in ci_methods) {
    audit_path <- audit_csv_path(dataname, alpha, test, ci_method)
    audit_data <- prep_audit_plot_data(
      read_audit_csv(audit_path), ci_method, dataname, x_lim, test = test,
      group_levels = group_levels
    )
    n_grp <- nrow(audit_data)

    lb_label_data <- dplyr::filter(audit_data, !is.na(lb_label))
    ub_label_data <- dplyr::filter(audit_data, !is.na(ub_label))
    ci_segment_data <- dplyr::filter(audit_data, is.finite(seg_lb), is.finite(seg_ub))
    lb_point_data <- dplyr::filter(audit_data, is.finite(lb_plot))
    ub_point_data <- dplyr::filter(audit_data, is.finite(ub_plot))

    PPV_plot <- ggplot(audit_data, aes(y = group)) +
      geom_segment(
        data = ci_segment_data,
        aes(x = seg_lb, xend = seg_ub, y = group, yend = group,
            color = line_color, linetype = line_type),
        lwd = 0.8, alpha = 0.7
      ) +
      geom_point(
        data = lb_point_data,
        aes(x = lb_plot, color = line_color), shape = 19, size = 2.5, alpha = 0.9
      ) +
      geom_text(
        data = lb_label_data,
        aes(x = lb_label_x, label = lb_label, color = line_color, hjust = lb_label_hjust),
        vjust = 1.8, size = 2.8, show.legend = FALSE
      ) +
      geom_point(
        data = ub_point_data,
        aes(x = ub_plot, color = line_color), shape = 19, size = 2.5, alpha = 0.9
      ) +
      geom_text(
        data = ub_label_data,
        aes(x = ub_label_x, label = ub_label, color = line_color, hjust = ub_label_hjust),
        vjust = 1.8, size = 2.8, show.legend = FALSE
      ) +
      geom_point(aes(x = DeltaG, color = line_color),
                 shape = 18, size = 3) +
      scale_color_identity() +
      scale_linetype_identity() +
      geom_vline(xintercept = 1, linetype = "solid",
                 color = "gray60", lwd = 0.6) +
      scale_x_continuous(
        limits = x_lim,
        breaks = x_breaks,
        labels = scales::number_format(accuracy = 0.01),
        expand = c(0, 0)
      ) +
      labs(x = "Postive Predictive Value Disparity", y = NULL,
           title = paste0(
             dataname, " Subgroups with ", alpha * 100, "% Confidence Intervals",
             if (ci_method == "CELCI") " (one-side)" else ""
           )) +
      theme_minimal() +
      theme(
        panel.border = element_rect(color = "black", fill = NA),
        axis.line.x = element_line(color = "black"),
        panel.grid.major.x = element_line(color = "gray90"),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 11, color = "black", margin = margin(r = 10)),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.margin = unit(c(1, 1, 1, 1), "cm")
      ) +
      scale_y_discrete(limits = rev(group_levels)) +
      coord_cartesian(clip = "off")

    print(PPV_plot)

    if (isTRUE(save_plot)) {
      ggsave(
        filename = audit_plot_path(dataname, alpha, test, ci_method),
        plot = PPV_plot,
        device = cairo_ps,
        width = 8.5,
        height = plot_height_in(n_grp),
        units = "in",
        dpi = 300,
        bg = "transparent"
      )
    }
  }
}
