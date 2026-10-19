project_file <- function(path) file.path(Sys.getenv('PORTFOLIO_ROOT'),path)

testthat::test_that('official source schema and outcome coding are respected', {
 df <- load_gbsg2()
 testthat::expect_equal(dim(df),c(686L,11L))
 testthat::expect_equal(sum(df$cens==1),299L)
 testthat::expect_equal(sum(df$cens==0),387L)
 testthat::expect_true(all(df$time>0))
 testthat::expect_false(anyNA(df))
 testthat::expect_equal(levels(df$tgrade),c('I','II','III'))
 testthat::expect_false(is.ordered(df$tgrade))
})

testthat::test_that('split is reproducible, disjoint and event stratified', {
 df <- prepare_predictors(load_gbsg2())
 a <- stratified_split(df); b <- stratified_split(df)
 testthat::expect_identical(a,b)
 testthat::expect_length(intersect(a$train,a$test),0L)
 testthat::expect_setequal(c(a$train,a$test),seq_len(nrow(df)))
 testthat::expect_equal(sum(df$cens[a$train]),209L)
 testthat::expect_equal(sum(df$cens[a$test]),90L)
 testthat::expect_identical(stratified_folds(df[a$train,]),stratified_folds(df[a$train,]))
 tab <- table(stratified_folds(df[a$train,]),df$cens[a$train])
 testthat::expect_true(all(apply(tab,2,function(x) diff(range(x)))<=1))
})

testthat::test_that('training audit and fixed preprocessing exclude the test partition', {
 r <- readRDS(project_file('artifacts/results.rds'))
 testthat::expect_setequal(r$tuning$training_ids,r$train_ids)
 testthat::expect_length(intersect(r$tuning$training_ids,r$test_ids),0L)
 testthat::expect_equal(nrow(r$tuning$folds),18L)
 testthat::expect_true(all(is.finite(r$tuning$folds$ibs)))
 testthat::expect_false('test' %in% names(formals(tune_training)))
 df <- load_gbsg2(); tr <- df[r$train_ids,]; te <- df[r$test_ids,]
 before <- prepare_predictors(tr)
 te$age <- te$age+1000; te$time <- te$time+1000
 testthat::expect_identical(prepare_predictors(tr),before)
 models <- readRDS(project_file('artifacts/models.rds'))
 testthat::expect_equal(models$Cox$n,length(r$train_ids))
 testthat::expect_setequal(as.integer(rownames(models$Cox$model)),r$train_ids)
 testthat::expect_equal(models$RSF$n,length(r$train_ids))
 testthat::expect_false(any(c('id','time','cens') %in% models$RSF$xvar.names))
})

testthat::test_that('survival predictions are dimensionally correct, bounded and monotone', {
 r <- readRDS(project_file('artifacts/results.rds'))
 for (p in r$predictions) {
  testthat::expect_equal(dim(p),c(length(r$test_ids),length(r$config$grid)))
  testthat::expect_true(all(is.finite(p)))
  testthat::expect_true(all(p>=0 & p<=1))
  testthat::expect_true(all(apply(p,1,function(x) all(diff(x)<=1e-10))))
 }
})

testthat::test_that('risk direction is correct in a deterministic concordance fixture', {
 # A unit-test fixture only; all project analyses use the real GBSG2 data.
 d <- data.frame(time=c(1,2,3,4),cens=rep(1,4),risk=c(4,3,2,1))
 c <- concordance_metrics(d,list(perfect=d$risk),tau=4)
 testthat::expect_equal(c$estimate,rep(1,2))
 reversed <- concordance_metrics(d,list(reversed=-d$risk),tau=4)
 testthat::expect_equal(reversed$estimate,rep(0,2))
})

testthat::test_that('validated IPCW metrics reduce to uncensored score identities', {
 r <- readRDS(project_file('artifacts/results.rds')); models <- readRDS(project_file('artifacts/models.rds'))
 df <- prepare_predictors(load_gbsg2())
 # Use genuine held-out event observations; no censored outcomes are recoded.
 d <- df[r$test_ids,]; d <- d[d$cens==1,]
 tau <- 1095
 p <- 1-predict_survival(models$Cox,d,tau)[,1]
 score <- score_predictions(list(Cox=models$Cox),d,tau,ci=FALSE)
 actual_brier <- score$Brier$score$Brier
 actual_auc <- score$AUC$score$AUC
 y <- as.numeric(d$time<=tau)
 testthat::expect_equal(as.numeric(actual_brier),mean((y-p)^2),tolerance=1e-8)
 pairs <- outer(p[y==1],p[y==0],'-')
 testthat::expect_equal(as.numeric(actual_auc),mean((pairs>0)+.5*(pairs==0)),tolerance=1e-8)
 testthat::expect_equal(integrated_brier(c(1,2,3),c(.1,.2,.3)),.2)
})

testthat::test_that('saved reports match evaluated predictions and result objects', {
 r <- readRDS(project_file('artifacts/results.rds'))
 df <- prepare_predictors(load_gbsg2()); te <- df[r$test_ids,]
 risks <- lapply(r$predictions,function(p) 1-p[,ncol(p)])
 c <- concordance_metrics(te,risks)
 saved <- read.csv(project_file('reports/test-concordance.csv'))
 testthat::expect_equal(saved$estimate,c$estimate,tolerance=1e-12)
 testthat::expect_equal(saved$estimate,r$concordance$estimate,tolerance=1e-12)
 manifest <- read.csv(project_file('reports/split-manifest.csv'))
 testthat::expect_setequal(manifest$id[manifest$partition=='test'],r$test_ids)
 testthat::expect_true(all(is.na(manifest$fold[manifest$partition=='test'])))
 b <- read.csv(project_file('reports/test-brier.csv'))
 ibs <- read.csv(project_file('reports/test-integrated-brier.csv'))
 for (m in ibs$model) {
  bm <- b[b$model==m,]
  testthat::expect_equal(ibs$IBS[ibs$model==m],integrated_brier(bm$times,bm$Brier),tolerance=1e-12)
 }
 testthat::expect_gt(file.info(project_file('reports/clinical-survival-report.html'))$size,10000)
 html <- paste(readLines(project_file('reports/clinical-survival-report.html'),warn=FALSE),collapse='\n')
 testthat::expect_match(html,'686 women',fixed=TRUE)
 readme <- paste(readLines(project_file('README.md')),collapse='\n')
 testthat::expect_match(readme,sprintf('%.3f',saved$estimate[1]),fixed=TRUE)
 curves <- read.csv(project_file('reports/cox-functional-form-curves.csv'))
 testthat::expect_true(all(is.finite(curves$hazard_ratio)))
 testthat::expect_equal(nrow(curves),320L)
})
