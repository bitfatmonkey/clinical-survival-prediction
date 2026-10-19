Surv <- survival::Surv
# Data contract is taken from help("GBSG2", package = "TH.data").
load_gbsg2 <- function() {
  e <- new.env(parent = emptyenv())
  data("GBSG2", package = "TH.data", envir = e)
  df <- e$GBSG2
  validate_gbsg2(df)
  df$id <- seq_len(nrow(df))
  # Unordered grade gives separate II vs I and III vs I hazard ratios.
  df$tgrade <- factor(df$tgrade, levels = c("I", "II", "III"), ordered = FALSE)
  df
}

validate_gbsg2 <- function(df) {
  expected <- c("horTh", "age", "menostat", "tsize", "tgrade", "pnodes",
                "progrec", "estrec", "time", "cens")
  stopifnot(nrow(df) == 686L, identical(names(df), expected), !anyNA(df),
            all(df$time > 0), setequal(unique(df$cens), 0:1),
            all(df$age > 0), all(df$tsize > 0), all(df$pnodes >= 0),
            all(df$progrec >= 0), all(df$estrec >= 0))
  invisible(TRUE)
}

prepare_predictors <- function(df) {
  # Fixed, prespecified transformations: no estimated preprocessing parameters.
  df$age10 <- df$age / 10
  df$tsize10 <- df$tsize / 10
  df$log_nodes <- log1p(df$pnodes)
  df$log_progrec <- log1p(df$progrec)
  df$log_estrec <- log1p(df$estrec)
  df
}

stratified_split <- function(df, fraction = 0.7, seed = 20261019L) {
  set.seed(seed)
  train <- sort(unlist(lapply(split(seq_len(nrow(df)), df$cens), function(i) {
    sample(i, floor(length(i) * fraction))
  }), use.names = FALSE))
  list(train = train, test = setdiff(seq_len(nrow(df)), train))
}

stratified_folds <- function(df, k = 3L, seed = 20261020L) {
  set.seed(seed)
  fold <- integer(nrow(df))
  for (indices in split(seq_len(nrow(df)), df$cens)) {
    fold[indices] <- sample(rep(seq_len(k), length.out = length(indices)))
  }
  fold
}
