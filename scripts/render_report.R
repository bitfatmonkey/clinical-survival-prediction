if (!file.exists('artifacts/results.rds')) stop('Run scripts/run_all.R before rendering.')
rmarkdown::render('notebooks/analysis.Rmd', output_file='clinical-survival-report.html',
                  output_dir='reports', knit_root_dir=normalizePath('.'),
                  envir=new.env(parent=globalenv()), quiet=TRUE)
source('scripts/write_readme.R')
