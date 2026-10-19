write_table <- function(x, name) {
  utils::write.csv(as.data.frame(x), file.path("reports", paste0(name, ".csv")), row.names = FALSE)
}

characteristics <- function(cohorts) {
  rows <- list()
  for (cohort in names(cohorts)) {
    df <- cohorts[[cohort]]
    rows[[length(rows) + 1L]] <- data.frame(cohort, variable = "Patients", level = "N", value = as.character(nrow(df)))
    for (var in c("age", "tsize", "pnodes", "progrec", "estrec", "time")) {
      q <- stats::quantile(df[[var]], c(.25, .5, .75))
      rows[[length(rows) + 1L]] <- data.frame(cohort, variable = var, level = "Median [IQR]",
        value = sprintf("%.1f [%.1f, %.1f]", q[2], q[1], q[3]))
    }
    for (var in c("horTh", "menostat", "tgrade", "cens")) {
      tab <- table(df[[var]])
      rows[[length(rows) + 1L]] <- data.frame(cohort, variable = var, level = names(tab),
        value = sprintf("%d (%.1f%%)", as.integer(tab), as.numeric(tab) / nrow(df) * 100))
    }
  }
  do.call(rbind, rows)
}

km_analysis <- function(df, horizons) {
  df$size_category <- cut(df$tsize, breaks = c(-Inf, 20, 50, Inf),
                           labels = c("≤20 mm", ">20–50 mm", ">50 mm"))
  fits <- list(
    overall = survival::survfit(Surv(time, cens) ~ 1, data = df),
    hormone = survival::survfit(Surv(time, cens) ~ horTh, data = df),
    menopause = survival::survfit(Surv(time, cens) ~ menostat, data = df),
    size = survival::survfit(Surv(time, cens) ~ size_category, data = df)
  )
  titles <- c(overall = "Overall recurrence-free survival", hormone = "Recurrence-free survival by hormonal therapy",
              menopause = "Recurrence-free survival by menopausal status", size = "Recurrence-free survival by tumour size")
  risks <- list(); probabilities <- list(); medians <- list()
  for (name in names(fits)) {
    fit <- fits[[name]]
    plot_km(fit, titles[[name]], paste0("km-", name))
    s <- summary(fit, times = c(0, horizons), extend = FALSE)
    risk <- data.frame(analysis = name, stratum = if (is.null(s$strata)) "Overall" else as.character(s$strata),
                       days = s$time, n_risk = s$n.risk)
    risks[[name]] <- risk
    probabilities[[name]] <- data.frame(risk, survival = s$surv, lower = s$lower, upper = s$upper)
    tab <- as.data.frame(if (is.null(dim(summary(fit)$table))) t(summary(fit)$table) else summary(fit)$table)
    tab$stratum <- if (is.null(fit$strata)) "Overall" else names(fit$strata)
    tab$analysis <- name
    medians[[name]] <- tab
  }
  risk <- do.call(rbind, risks)
  save_plot(ggplot2::ggplot(risk, ggplot2::aes(days / 365, stratum, label = n_risk)) +
    ggplot2::geom_text() + ggplot2::facet_grid(analysis ~ ., scales = "free_y", space = "free_y") +
    portfolio_theme() + ggplot2::scale_x_continuous(breaks = c(0, horizons / 365)) +
    ggplot2::labs(title = "Numbers at risk", x = "Time (years)", y = NULL), "km-risk-table", height = 8)
  lr <- lapply(c("horTh", "menostat", "size_category"), function(v) {
    f <- stats::as.formula(paste("Surv(time, cens) ~", v))
    z <- survival::survdiff(f, data = df)
    data.frame(variable = v, chisq = z$chisq, df = length(z$n) - 1,
                p_value = stats::pchisq(z$chisq, length(z$n) - 1, lower.tail = FALSE))
  })
  write_table(do.call(rbind, lr), "logrank-tests")
  write_table(risk, "numbers-at-risk")
  write_table(do.call(rbind, probabilities), "km-probabilities")
  write_table(do.call(rbind, medians), "km-medians")
  invisible(fits)
}

cox_diagnostics <- function(train, selected, linear, spline) {
  terms <- c("age10", "tsize10", "log_nodes", "log_progrec", "log_estrec", "horTh", "menostat", "tgrade")
  uni <- do.call(rbind, lapply(terms, function(v) {
    fit <- survival::coxph(stats::as.formula(paste("Surv(time, cens) ~", v)), data = train)
    out <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    out$predictor <- v; out
  }))
  write_table(uni, "cox-univariate-training")
  adjusted <- broom::tidy(linear, exponentiate = TRUE, conf.int = TRUE)
  write_table(adjusted, "cox-adjusted-linear-training")
  write_table(broom::tidy(selected, exponentiate = TRUE, conf.int = TRUE), "cox-selected-coefficients")
  zph <- survival::cox.zph(selected)
  ph <- as.data.frame(zph$table); ph$term <- rownames(ph)
  write_table(ph, "cox-ph-test")
  grDevices::png("figures/cox-schoenfeld.png", width = 1600, height = 1800, res = 160, type = "cairo")
  graphics::par(mfrow = c(4, 2), mar = c(4, 4, 3, 1))
  plot(zph)
  grDevices::dev.off()
  lrt <- as.data.frame(stats::anova(linear, spline, test = "LRT"))
  lrt$model <- c("linear", "spline")
  write_table(lrt, "cox-nonlinearity-lrt")
  save_plot(ggplot2::ggplot(adjusted, ggplot2::aes(estimate, reorder(term, estimate))) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = .2) +
    ggplot2::geom_point(colour = "#197a8c") + ggplot2::scale_x_log10() + portfolio_theme() +
    ggplot2::labs(title = "Adjusted hazard ratios: prespecified linear Cox model",
                  subtitle = "Training data only; 95% confidence intervals; transformations defined in methodology",
                  x = "Hazard ratio (log scale)", y = NULL), "cox-forest", height = 6)
  # Compare linear and spline age/size shapes at a fixed training reference profile.
  ref <- train[1, , drop = FALSE]
  for (v in c("age", "tsize", "pnodes", "progrec", "estrec")) ref[[v]] <- median(train[[v]])
  ref$horTh <- factor("no", levels = levels(train$horTh))
  ref$menostat <- factor("Post", levels = levels(train$menostat))
  ref$tgrade <- factor("II", levels = levels(train$tgrade))
  curves <- list()
  for (v in c("age", "tsize")) {
    nd <- ref[rep(1, 80), ]
    nd[[v]] <- seq(quantile(train[[v]], .05), quantile(train[[v]], .95), length.out = 80)
    nd <- prepare_predictors(nd)
    for (name in c("Linear", "Spline")) {
      model <- if (name == "Linear") linear else spline
      p <- stats::predict(model, newdata = nd, type = "lp")
      baseline <- stats::predict(model, newdata = prepare_predictors(ref), type = "lp")
      curves[[length(curves) + 1L]] <- data.frame(variable = v, value = nd[[v]],
        hazard_ratio = exp(p - baseline), model = name)
    }
  }
  curves <- do.call(rbind, curves)
  write_table(curves, 'cox-functional-form-curves')
  stopifnot(all(is.finite(curves$hazard_ratio)))
  save_plot(ggplot2::ggplot(curves, ggplot2::aes(value, hazard_ratio, colour = model)) +
    ggplot2::geom_line() + ggplot2::facet_wrap(~variable, scales = "free_x") +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) + portfolio_theme() +
    ggplot2::labs(title = "Training-only functional-form sensitivity", x = "Age (years) or tumour size (mm)",
                  y = "Hazard ratio versus training median", colour = NULL), "cox-functional-forms")
  invisible(list(ph = ph, adjusted = adjusted))
}
