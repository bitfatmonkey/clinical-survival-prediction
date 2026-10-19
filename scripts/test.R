if (!file.exists('artifacts/results.rds')) stop('Execute the analysis first.')
for (file in sort(list.files('R', pattern='\\.R$', full.names=TRUE))) source(file)
Sys.setenv(PORTFOLIO_ROOT=normalizePath('.'))
test_results <- testthat::test_dir('tests/testthat', reporter='summary', stop_on_failure=TRUE)
test_summary <- as.data.frame(test_results)
utils::write.csv(test_summary[,c('file','test','nb','failed','skipped','error','warning','passed')],
                 'reports/test-results.csv',row.names=FALSE)
