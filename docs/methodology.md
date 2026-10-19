# Analysis protocol and methodological choices

## Question and source

The objective is prognostic prediction of recurrence-free survival in GBSG2, not mortality prediction or causal treatment-effect estimation. TH.data's installed GBSG2 documentation is extracted by the executed pipeline and committed with its source references. The dataset contains 686 women, 299 recorded events and 387 right-censored observations. The extract does not provide event-type detail sufficient to separate recurrence and death components. No endpoint component or receptor unit denominator is inferred beyond the documentation.

## Data quality and preprocessing

The ten source columns, row count, missingness, positive follow-up, binary outcome coding and nonnegative clinical measurements are validated. There are no missing observations, so no imputation or exclusions are applied. Receptor values of zero are retained. Extreme observations are retained rather than winsorised. Grade is converted from ordered I/II/III to an unordered factor with I as reference. Row IDs are deterministic row positions and are excluded from models.

All preprocessing is prespecified and pointwise: age/10, size/10, and log(1+x) for positive nodes and both receptors. These transforms have no parameters learned from the full cohort. The Cox spline candidate learns age/size knots only within each training fold. Forests receive the same underlying eight predictors, including these monotone transforms; the transformation does not create extra clinical information.

## Holdout and training-only selection

A seed-20261019 event-stratified split assigns floor(70% of each event stratum) to training: 479 observations and 209 events. The untouched test set has 207 observations and 90 events. Each training subject is assigned to one of three event-stratified folds with seed 20261020. The split manifest records row IDs, partition and fold; test rows have no training-fold assignment.

Candidate selection receives only the training partition and its fold IDs. Each candidate is fitted three times, without its corresponding validation fold. No variable selection or test-set tuning occurs. The criterion is mean cross-validated IPCW integrated Brier score over days 365–1825, approximated by a 13-point trapezoid rule. This restricted interval deliberately avoids unsupported late follow-up; it is not zero-to-five-year IBS. The small grid is fixed before final test evaluation.

Cox candidates:

1. Linear age and size; log(1+x) nodes/receptors; hormonal therapy, menopausal status and nominal grade (9 coefficients).
2. Natural cubic splines with 3 df each for age and size; other terms unchanged (13 coefficients).

Forest candidates cross mtry {2,4} with node size {10,25}, using 200 trees, log-rank splitting, 10 random split points and fixed fold-specific seeds. The selected forest is refitted with 500 trees and training OOB permutation importance. The same seed is used for the final fit on repeated runs; two CPU/OpenMP threads are requested. Floating-point and parallel reductions can vary slightly across architectures.

With 209 training events, this is a deliberately modest Cox parameter count; events per coefficient is a rough diagnostic, not a formal overfitting guarantee. Penalised Cox was not added because it would add another model family and tuning exercise to a small comparison. A shrinkage analysis with nested resampling is a future extension, not a missing unreported result.

## Descriptive and inferential analysis

Full-cohort EDA and Kaplan–Meier are descriptive and not used to select predictors. Size categories use <=20 mm, >20–50 mm and >50 mm conventional cutpoints, not complete TNM stages. Overall and grouped survival, pointwise 95% confidence intervals, risk tables, medians and exploratory log-rank tests are exported. An unestimable median or upper limit is recorded as missing, not zero. Risk sets fall over time; tail interpretation requires caution.

Univariate and adjusted Cox association tables use training data. The adjusted HR forest plot always depicts the linear specification so coefficients have interpretable units. If splines win prediction selection, their basis coefficients are exported but not interpreted as unit clinical effects. Continuous shape plots compare the two candidates at training reference values, over the central 90% of observed training age/size. A training-only likelihood-ratio test compares functional forms as sensitivity analysis. It is not another test-set selection rule.

Age and size HRs are per 10 units. A log(1+x) HR refers to multiplying (1+x) by e. The HR for doubling (1+x) is exp(beta*log(2)). Reference factor levels are no hormonal therapy, premenopausal status and grade I. Contrasting crude and adjusted HRs can show sensitivity to adjustment, but changes can also reflect non-collapsibility; no causal effect is estimated.

Schoenfeld residuals and cox.zph assess PH on training data for the selected specification. PH violations are reported without post-hoc test-driven model revision. HRs and absolute predictions remain model-dependent. P-values and intervals are exploratory, conditional on specification and unadjusted for multiple comparisons or selection.

## Censoring-aware held-out evaluation

Both final models are assessed once on the same 207 held-out subjects, event coding and grid. No fitted model is altered after evaluation. The marginal censoring model used to assess a fold or the test set is estimated within that assessment sample. This is standard IPCW performance estimation, not preprocessing or model selection on test outcomes. It assumes independent censoring; conditional or informative censoring can invalidate the weights.

- **Harrell C:** survival::concordance with reverse=TRUE and timewt='n', ranking both models by five-year predicted event probability. Comparable pairs cover observed follow-up.
- **Uno C:** the same ranking with survival's documented timewt='n/G2', truncated using ymax=1825. This is explicitly a five-year truncated statistic.
- **AUC:** riskRegression::Score cumulative/dynamic time-dependent AUC, with observed events by t as cases and event-free beyond t as controls, IPCW-adjusted for censoring. Pointwise influence-function intervals are exported.
- **Brier:** riskRegression::Score IPCW Brier score at each grid time, not a binary event-status MSE that ignores censoring. A Kaplan–Meier benchmark is fitted in training, not in test. AUC for that constant-risk benchmark is uninformative (0.5).
- **IBS:** time-average of Brier over days 365–1825 with trapezoidal integration on 13 points; no IBS interval is claimed.
- **Calibration:** within-model test risk tertiles at five years; mean predicted event risk versus 1-KM with 95% confidence limits and group risk-set sizes. This is coarse descriptive calibration, not recalibration or a validated clinical threshold rule.

Cox and forest risk predictions both use riskRegression::predictRisk. The forest method is a censoring-aware survival prediction interface; mortality scores are not probabilities. All curves are checked for correct dimensions, finiteness, probability bounds and non-increasing survival over time.

At one, three and five years the test at-risk counts are 184, 97 and 35; events by those times are 18, 74 and 87. Marginal censoring survival is approximately .976, .757 and .336. The pipeline requires >=20 at risk, >=5 events and censoring survival >.1 at the predeclared reporting horizons; these are practical support checks, not proof of stable estimation. No extrapolation beyond five years is evaluated.

Concordance intervals use survival's infinitesimal-jackknife variance. A paired Cox-minus-RSF Harrell comparison uses the joint covariance to respect paired subjects. Pointwise AUC and Brier intervals come from Score. These intervals condition on the fitted models and do not incorporate all uncertainty from fitting/tuning, splitting or selection. No simultaneous coverage, multiplicity control or definitive superiority claim is implied.

## Quality assurance

Tests check source schema and event count, split reproducibility and disjointness, training audit, fold balance, absence of outcome/ID predictors, correct prediction shape and monotonicity, risk direction, and report consistency. On genuine uncensored test observations only, the validated library's Brier and AUC are cross-checked against their no-censoring MSE and pairwise-ranking identities. A deterministic four-observation concordance fixture is used only as a unit test; project analyses use no synthetic data.

Actual figures, CSV tables, session information, package versions and the R Markdown HTML report are committed. Model objects remain ignored. Re-executing the pipeline regenerates the evidence rather than reading fabricated performance constants. Remote CI status and independent external validation must be reported separately.
