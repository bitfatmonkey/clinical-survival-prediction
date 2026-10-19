portfolio_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
}

save_plot <- function(plot, name, width = 8, height = 5) {
  ggplot2::ggsave(file.path("figures", paste0(name, ".png")), plot,
                  width = width, height = height, dpi = 160, bg = "white", type = "cairo")
  invisible(plot)
}

plot_eda <- function(df) {
  labels <- c(age = "Age (years)", tsize = "Tumour size (mm)",
              pnodes = "Positive nodes (count)", progrec = "Progesterone receptor (fmol)",
              estrec = "Oestrogen receptor (fmol)", time = "Observed follow-up (days)")
  continuous <- tidyr::pivot_longer(df[, names(labels)], cols = dplyr::everything())
  continuous$name <- factor(continuous$name, levels = names(labels), labels = labels)
  save_plot(ggplot2::ggplot(continuous, ggplot2::aes(value)) +
    ggplot2::geom_histogram(bins = 30, fill = "#197a8c", colour = "white") +
    ggplot2::facet_wrap(~name, scales = "free", ncol = 2) + portfolio_theme() +
    ggplot2::labs(title = "Clinical measurements and observed follow-up", x = NULL, y = "Patients"),
    "clinical-distributions", height = 9)
  categorical <- data.frame(
    variable = rep(c("Hormonal therapy", "Menopausal status", "Tumour grade", "Outcome"), each = nrow(df)),
    value = c(as.character(df$horTh), as.character(df$menostat), as.character(df$tgrade),
              ifelse(df$cens == 1, "Event", "Censored")))
  save_plot(ggplot2::ggplot(categorical, ggplot2::aes(value)) +
    ggplot2::geom_bar(fill = "#555e89") + ggplot2::facet_wrap(~variable, scales = "free_x") +
    portfolio_theme() + ggplot2::labs(title = "Patient groups and observed event status", x = NULL, y = "Patients"),
    "clinical-categories", height = 6)
}

plot_km <- function(fit, title, name) {
  d <- broom::tidy(fit, conf.int = TRUE)
  if (!"strata" %in% names(d)) d$strata <- "Overall"
  zeros <- d[!duplicated(d$strata), ]
  zeros$time <- 0; zeros$estimate <- 1; zeros$conf.low <- 1; zeros$conf.high <- 1
  d <- rbind(zeros, d)
  d <- d[order(d$strata, d$time), ]
  save_plot(ggplot2::ggplot(d, ggplot2::aes(time / 365, estimate, colour = strata, fill = strata)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high), alpha = .12, colour = NA) +
    ggplot2::geom_step(linewidth = .65) + portfolio_theme() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = title, subtitle = "Kaplan–Meier estimates with 95% confidence intervals",
                  x = "Time (years; 365 days)", y = "Recurrence-free survival probability", colour = NULL, fill = NULL), name)
}

plot_performance <- function(concordance, auc, brier, calibration) {
  save_plot(ggplot2::ggplot(concordance, ggplot2::aes(model, estimate, colour = model)) +
    ggplot2::geom_hline(yintercept = .5, linetype = 2, colour = "grey60") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper)) +
    ggplot2::facet_wrap(~metric) + ggplot2::coord_cartesian(ylim = c(.4, 1)) + portfolio_theme() +
    ggplot2::labs(title = "Held-out discrimination", subtitle = "95% influence-function confidence intervals; higher is better",
                  x = NULL, y = "Concordance"), "concordance-comparison")
  for (item in c("AUC", "Brier")) {
    d <- if (item == "AUC") auc else brier
    save_plot(ggplot2::ggplot(d, ggplot2::aes(times / 365, .data[[item]], colour = model)) +
      ggplot2::geom_line() + ggplot2::geom_point() + portfolio_theme() +
      ggplot2::labs(title = paste("Held-out", if (item == "AUC") "time-dependent AUC" else "IPCW Brier score"),
                    subtitle = if (item == "AUC") "Cumulative/dynamic cases and controls; higher is better" else "Includes training Kaplan–Meier benchmark; lower is better",
                    x = "Time (years)", y = item), paste0(tolower(item), "-over-time"))
  }
  save_plot(ggplot2::ggplot(calibration, ggplot2::aes(predicted_risk, observed_risk, colour = model)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper)) +
    ggplot2::geom_line() + portfolio_theme() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = "Five-year calibration in test-set risk tertiles",
                  subtitle = "Observed risks: Kaplan–Meier with 95% confidence intervals; descriptive and imprecise",
                  x = "Mean predicted event probability", y = "Observed event probability"), "calibration")
}
