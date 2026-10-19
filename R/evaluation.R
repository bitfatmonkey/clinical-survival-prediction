score_predictions <- function(models, data, times, ci = TRUE) {
  # IPCW and influence-function standard errors come from riskRegression::Score.
  # Censoring nuisance distribution is estimated in each assessment sample,
  # not used for model fitting or selection outside the training CV.
  riskRegression::Score(
    object = models, formula = Surv(time, cens) ~ 1,
    data = data, times = times, metrics = c("auc", "brier"),
    cens.model = "km", cens.method = "ipcw", null.model = FALSE,
    conf.int = if (ci) 0.95 else FALSE, contrasts = FALSE,
    split.method = "none", verbose = 0
  )
}

integrated_brier <- function(times, brier) {
  ord <- order(times)
  times <- times[ord]; brier <- brier[ord]
  stopifnot(length(times) > 1L, !anyDuplicated(times), all(is.finite(brier)))
  sum(diff(times) * (head(brier, -1) + tail(brier, -1)) / 2) /
    diff(range(times))
}

concordance_metrics <- function(data, risks, tau = 1825) {
  rows <- list()
  for (model in names(risks)) {
    d <- data; d$risk <- risks[[model]]
    for (method in c("Harrell", "Uno (5-year truncated)")) {
      fit <- if (method == "Harrell") {
        survival::concordance(Surv(time, cens) ~ risk, data = d,
                              reverse = TRUE, timewt = "n")
      } else {
        survival::concordance(Surv(time, cens) ~ risk, data = d,
                              reverse = TRUE, timewt = "n/G2", ymax = tau)
      }
      se <- sqrt(fit$var)
      rows[[length(rows) + 1L]] <- data.frame(
        model = model, metric = method, estimate = fit$concordance,
        lower = max(0, fit$concordance - 1.96 * se),
        upper = min(1, fit$concordance + 1.96 * se)
      )
    }
  }
  do.call(rbind, rows)
}

followup_support <- function(data, horizons) {
  g <- survival::survfit(Surv(time, 1 - cens) ~ 1, data = data)
  gs <- summary(g, times = horizons, extend = FALSE)$surv
  stopifnot(length(gs) == length(horizons))
  out <- data.frame(days = horizons,
                    at_risk = vapply(horizons, function(t) sum(data$time >= t), 1L),
                    events_by_time = vapply(horizons, function(t) {
                      sum(data$time <= t & data$cens == 1)
                    }, 1L), censor_survival = gs)
  stopifnot(all(out$at_risk >= 20), all(out$events_by_time >= 5),
            all(out$censor_survival > 0.1))
  out
}

calibration_groups <- function(data, risk, horizon = 1825, model) {
  # Tertiles describe calibration on the assessment sample; no model is updated.
  group <- cut(rank(risk, ties.method = "first"),
               breaks = c(0, nrow(data)/3, 2*nrow(data)/3, nrow(data)), labels = FALSE)
  do.call(rbind, lapply(sort(unique(group)), function(g) {
    d <- data[group == g, , drop = FALSE]
    km <- survival::survfit(Surv(time, cens) ~ 1, data = d)
    s <- summary(km, times = horizon, extend = FALSE)
    stopifnot(length(s$surv) == 1)
    data.frame(model = model, group = g, n = nrow(d), at_risk = s$n.risk,
               predicted_risk = mean(risk[group == g]), observed_risk = 1 - s$surv,
               lower = 1 - s$upper, upper = 1 - s$lower)
  }))
}
