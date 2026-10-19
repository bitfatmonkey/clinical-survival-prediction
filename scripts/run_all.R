#!/usr/bin/env Rscript
# Run from the repository root: Rscript scripts/run_all.R
if (!file.exists('DESCRIPTION') || !dir.exists('R')) stop('Run from the repository root.')
source('renv/activate.R')
options(warn = 1)
Sys.setenv(OMP_NUM_THREADS = '2')
for (file in sort(list.files('R', pattern = '\\.R$', full.names = TRUE))) source(file)
for (dir in c('reports','figures','artifacts')) dir.create(dir, showWarnings = FALSE)
started <- Sys.time()
config <- list(seed = 20261019L, train_fraction = .7, folds = 3L,
               horizons = c(365,1095,1825), grid = seq(365,1825,length.out=13),
               tuning_trees = 200L, final_trees = 500L)
# Read and retain the official installed documentation as provenance evidence.
h <- utils:::.getHelpFile(help('GBSG2', package='TH.data'))
text <- capture.output(tools::Rd2txt(h))
text <- gsub('.\b', '', text)
writeLines(text, 'docs/GBSG2-official-documentation.txt')
df <- prepare_predictors(load_gbsg2())
split <- stratified_split(df, fraction=config$train_fraction, seed=config$seed)
train <- df[split$train, , drop=FALSE]
test <- df[split$test, , drop=FALSE]
fold <- stratified_folds(train,k=config$folds)
stopifnot(!length(intersect(train$id,test$id)))
write_table(data.frame(id=df$id, partition=ifelse(df$id %in% train$id,'training','test'),
                       fold=ifelse(df$id %in% train$id,fold[match(df$id,train$id)],NA)), 'split-manifest')
write_table(data.frame(variable=names(df), type=vapply(df,function(x) paste(class(x),collapse='/'),''),
                        missing=vapply(df,function(x) sum(is.na(x)),1L)), 'data-quality')
write_table(characteristics(list(Overall=df,Training=train,Test=test)), 'patient-characteristics')
write_table(followup_support(test,config$horizons), 'test-followup-support')
plot_eda(df)
km_analysis(df,config$horizons)
# All candidate selection is completed before the test models are scored.
tuning <- tune_training(train,fold,config$grid,trees=config$tuning_trees)
write_table(tuning$candidates,'training-cv-candidates')
write_table(tuning$folds,'training-cv-fold-scores')
cox <- fit_cox(train,tuning$best$Cox$shape)
rsf <- fit_rsf(train,tuning$best$RSF$mtry,tuning$best$RSF$nodesize,
               ntree=config$final_trees,importance='permute')
models <- list(Cox=cox,RSF=rsf)
linear <- fit_cox(train,'linear'); spline <- fit_cox(train,'spline')
cox_diag <- cox_diagnostics(train,cox,linear,spline)
# Final evaluation on the same untouched cases and the same time grid.
km_benchmark <- survival::survfit(Surv(time,cens) ~ 1, data=train)
score <- score_predictions(c(models,list(KM=km_benchmark)),test,config$grid,ci=TRUE)
auc <- as.data.frame(score$AUC$score)
brier <- as.data.frame(score$Brier$score)
write_table(auc,'test-time-dependent-auc')
write_table(brier,'test-brier')
ibs <- do.call(rbind,lapply(split(brier,brier$model),function(d) {
 data.frame(model=as.character(d$model[1]),from=365,to=1825,IBS=integrated_brier(d$times,d$Brier))
}))
write_table(ibs,'test-integrated-brier')
predictions <- lapply(models,predict_survival,newdata=test,times=config$grid)
risks <- lapply(predictions,function(p) 1-p[,ncol(p)])
concordance <- concordance_metrics(test,risks)
write_table(concordance,'test-concordance')
# Paired influence-function comparison accounts for using the same test patients.
pair_data <- test; pair_data$cox_risk <- risks$Cox; pair_data$rsf_risk <- risks$RSF
paired <- survival::concordance(Surv(time,cens) ~ cox_risk + rsf_risk,
                               data=pair_data, reverse=TRUE)
contrast <- c(1,-1)
difference <- sum(contrast*paired$concordance)
se <- sqrt(as.numeric(t(contrast) %*% paired$var %*% contrast))
write_table(data.frame(contrast='Cox minus RSF Harrell C',estimate=difference,
                       lower=difference-1.96*se,upper=difference+1.96*se),'test-paired-concordance')
calibration <- do.call(rbind,lapply(names(risks),function(name) calibration_groups(test,risks[[name]],model=name)))
write_table(calibration,'test-calibration')
plot_performance(concordance,auc,brier,calibration)
importance <- data.frame(variable=names(rsf$importance),importance=as.numeric(rsf$importance))
write_table(importance,'rsf-training-importance')
save_plot(ggplot2::ggplot(importance,ggplot2::aes(importance,reorder(variable,importance))) +
 ggplot2::geom_col(fill='#197a8c') + portfolio_theme() +
 ggplot2::labs(title='Training out-of-bag permutation importance',
               subtitle='Increase in forest prediction error; correlated predictors can share importance',
               x='Permutation importance',y=NULL),'rsf-importance')
write_table(data.frame(trees=seq_along(rsf$err.rate),oob_error=as.numeric(rsf$err.rate)), 'rsf-oob-error')
# Select genuine training patients representing 20th/50th/80th percentiles of
# training Cox risk; these are illustrations, not recommendations.
train_risk <- 1-predict_survival(cox,train,1825)[,1]
ids <- vapply(quantile(train_risk,c(.2,.5,.8)),function(q) which.min(abs(train_risk-q)),1L)
profiles <- train[ids, ,drop=FALSE]
profiles$profile <- c('Lower training risk','Middle training risk','Higher training risk')
write_table(profiles[,c('id','profile','age','tsize','tgrade','pnodes','progrec','estrec','horTh','menostat')],'illustrative-training-profiles')
curves <- do.call(rbind,lapply(names(models),function(name) {
 times <- seq(0,1825,by=25)
 p <- predict_survival(models[[name]],profiles,times)
 data.frame(model=name,profile=rep(profiles$profile,each=length(times)),
             time=rep(times,nrow(profiles)),survival=as.vector(t(p)))
}))
write_table(curves,'illustrative-survival-curves')
save_plot(ggplot2::ggplot(curves,ggplot2::aes(time/365,survival,colour=profile)) +
 ggplot2::geom_step() + ggplot2::facet_wrap(~model) + portfolio_theme() +
 ggplot2::scale_y_continuous(limits=c(0,1)) +
 ggplot2::labs(title='Predicted survival for three genuine training profiles',
               subtitle='Selected from training Cox risk percentiles; not clinical advice',
               x='Time (years)',y='Predicted recurrence-free survival',colour=NULL),'individual-survival')
result <- list(config=config,split=split,train_ids=train$id,test_ids=test$id,
               tuning=tuning,selected_cox=tuning$best$Cox$shape,
               test_events=sum(test$cens),train_events=sum(train$cens),
               coefficients=length(stats::coef(cox)),concordance=concordance,
               ibs=ibs,auc=auc,brier=brier,calibration=calibration,predictions=predictions,
               ph=cox_diag$ph,elapsed=as.numeric(difftime(Sys.time(),started,units='secs')))
saveRDS(result,'artifacts/results.rds')
saveRDS(models,'artifacts/models.rds')
saveRDS(score,'artifacts/test-score.rds')
writeLines(capture.output(sessionInfo()),'reports/session-info.txt')
packages <- c('TH.data','survival','dplyr','tidyr','ggplot2','broom','randomForestSRC','riskRegression','pec','testthat','renv','knitr','rmarkdown')
write_table(data.frame(package=packages,version=vapply(packages,function(p) as.character(packageVersion(p)),'')), 'package-versions')
source('scripts/render_report.R')
source('scripts/test.R')
cat('\nANALYSIS, REPORT AND TESTS COMPLETED\n')
print(concordance); print(ibs)
