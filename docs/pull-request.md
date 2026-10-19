## Summary

Adds an executable R-only portfolio comparing Cox regression and Random Survival Forests for recurrence-free survival on the real TH.data GBSG2 dataset. The pipeline validates source definitions, creates a fixed event-stratified holdout, tunes exclusively inside training folds, and generates actual clinical EDA, Kaplan–Meier analyses, Cox inference/diagnostics, forest interpretation and censoring-aware performance reports.

Includes reusable R functions, a single analysis entry point, renv.lock, a detailed R Markdown HTML report, 16 generated figures, CSV evidence, methodology/limitations documentation and an R 4.5.3 GitHub Actions workflow. Fitted model objects remain ignored.

## Observed results

- 686 patients / 299 events; training 479 / 209 events; test 207 / 90 events.
- Training-only three-fold CV selected linear Cox and RSF mtry=4 / nodesize=25; final forest has 500 trees.
- Held-out Harrell C: Cox 0.6809 (95% CI 0.6303–0.7314); RSF 0.6866 (0.6338–0.7394).
- Uno C truncated at five years: Cox 0.6763; RSF 0.6812.
- IPCW integrated Brier score over days 365–1825: Cox 0.1911; RSF 0.1894; training-fitted KM benchmark 0.2158.
- Paired Cox-minus-RSF C difference -0.0057 (95% CI -0.0367 to 0.0252); no convincing evidence of general model superiority.

## Validation and reproduction

Executed `Rscript scripts/run_all.R` using R 4.5.3 through the cloud runtime wrapper: exit 0, analysis and HTML rendering completed, 7 test groups / 49 passing assertions, no failed or skipped tests. Tested schema/event coding, split separation/reproducibility, training-only audit, prediction dimensions/bounds/monotonicity, risk direction, no-censoring scoring identities and report consistency. Repeated execution reproduced split and key metrics.

From the root, restore packages with `renv::restore()` and run `Rscript scripts/run_all.R`. See docs/reproducibility.md for native dependencies and cloud/macOS instructions, and reports/validation.md for detailed execution evidence and non-failing warnings. CI is configured; remote Actions execution is not claimed as passed.

## Limitations

Educational retrospective analysis, not a clinically validated system. Historical node-positive cohort; no external/NHS validation. Cox global PH test p=0.0201; late follow-up and calibration are sparse (35 held-out patients at risk at five years). Marginal independent-censoring assumptions, single-split variability and conditional uncertainty intervals limit inference. Native Apple Silicon restoration has not been executed. Pandoc deprecation and renv helper shutdown warnings are recorded in validation evidence.
