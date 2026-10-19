cox_formula <- function(shape = "linear") {
  if (shape == "linear") {
    Surv(time, cens) ~ age10 + tsize10 + log_nodes + log_progrec +
      log_estrec + horTh + menostat + tgrade
  } else {
    Surv(time, cens) ~ splines::ns(age10, df = 3) +
      splines::ns(tsize10, df = 3) + log_nodes + log_progrec + log_estrec +
      horTh + menostat + tgrade
  }
}

fit_cox <- function(train, shape = "linear") {
  survival::coxph(cox_formula(shape), data = train, x = TRUE, y = TRUE,
                  model = TRUE, ties = "efron")
}

fit_rsf <- function(train, mtry = 2L, nodesize = 15L, ntree = 500L,
                    seed = 20261021L, importance = FALSE) {
  randomForestSRC::rfsrc(
    Surv(time, cens) ~ age10 + tsize10 + log_nodes + log_progrec +
      log_estrec + horTh + menostat + tgrade,
    data = train, ntree = ntree, mtry = mtry, nodesize = nodesize,
    nsplit = 10, splitrule = "logrank", nthread = 2,
    seed = -abs(as.integer(seed)), importance = importance,
    forest = TRUE, na.action = "na.omit"
  )
}

predict_survival <- function(model, newdata, times) {
  risk <- riskRegression::predictRisk(model, newdata = newdata, times = times)
  risk <- matrix(risk, nrow = nrow(newdata), ncol = length(times))
  stopifnot(all(is.finite(risk)), all(risk >= 0 & risk <= 1))
  1 - risk
}

# This function accepts only the training partition, never the test partition.
tune_training <- function(train, fold, times, trees = 200L) {
  candidates <- data.frame(
    model = c(rep("Cox", 2), rep("RSF", 4)),
    shape = c("linear", "spline", rep(NA_character_, 4)),
    mtry = c(NA, NA, 2, 2, 4, 4),
    nodesize = c(NA, NA, 10, 25, 10, 25)
  )
  scores <- list()
  for (i in seq_len(nrow(candidates))) {
    for (f in sort(unique(fold))) {
      tr <- train[fold != f, , drop = FALSE]
      va <- train[fold == f, , drop = FALSE]
      fit <- if (candidates$model[i] == "Cox") {
        fit_cox(tr, candidates$shape[i])
      } else {
        fit_rsf(tr, candidates$mtry[i], candidates$nodesize[i], ntree = trees,
                seed = 20261021L + f)
      }
      score <- score_predictions(list(candidate = fit), va, times, ci = FALSE)
      brier <- as.data.frame(score$Brier$score)
      brier <- brier[brier$model == "candidate", ]
      scores[[length(scores) + 1L]] <- data.frame(
        candidate = i, fold = f, ibs = integrated_brier(brier$times, brier$Brier)
      )
      message("CV candidate ", i, "/", nrow(candidates), ", fold ", f,
              ": IBS = ", signif(tail(scores, 1)[[1]]$ibs, 4))
    }
  }
  scores <- do.call(rbind, scores)
  means <- aggregate(ibs ~ candidate, scores, mean)
  candidates$candidate <- seq_len(nrow(candidates))
  candidates <- merge(candidates, means, by = "candidate", sort = TRUE)
  best <- lapply(split(candidates, candidates$model), function(x) {
    x[order(x$ibs, x$candidate)[1], , drop = FALSE]
  })
  list(candidates = candidates, folds = scores, best = best,
       training_ids = train$id)
}
