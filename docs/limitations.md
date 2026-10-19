# Limitations and clinical scope

- This is a retrospective educational portfolio, not a clinically validated tool or medical advice. It is not deployed for patient care.
- GBSG2 is a historical cohort of women with node-positive breast cancer. Contemporary UK/NHS transportability, different stages, sexes, ethnicities, calendar periods and treatment pathways are not established.
- The documentation identifies recurrence-free survival and 0/1 event coding, but the supplied extract does not provide separate event components. Death and recurrence must not be inferred from this binary column. Competing risks cannot be disentangled here.
- Missing modern biological/treatment predictors, selection into the original study and historical care limit prediction. The data do not establish hospital mortality risk.
- Prognostic HRs, treatment-group KM curves and variable importance do not estimate causal treatment effects. Univariate-to-adjusted changes can reflect confounding, correlation or non-collapsibility.
- Marginal censoring weights and Kaplan–Meier require independent censoring. A covariate-dependent censoring mechanism would require a justified conditional model and sensitivity analysis.
- A small single train/test split gives uncertain comparisons. Cross-validation was restricted to training but is not a substitute for external validation or a repeated nested validation study.
- PH findings and continuous functional forms constrain Cox interpretation. The model was not revised after test performance was observed. A PH departure does not prove a forest is better calibrated or clinically useful.
- Forest importance can be unstable or shared among correlated predictors. Illustrative training profiles are not patient recommendations. Out-of-bag performance is a training diagnostic.
- At five years only 35 test patients remain at risk. Calibration groups are coarse, their uncertainty can be large, and late follow-up weights can be unstable. Horizons beyond five years were not evaluated.
- Concordance and AUC assess ranking, not absolute risk accuracy. Brier and grouped calibration add complementary information but no decision thresholds, utility or net clinical benefit were validated.
- CIs condition on fitted models and do not incorporate the full uncertainty in splitting, tuning or training. Exploratory log-rank/Cox tests are not multiplicity-adjusted. IBS has no claimed CI.
- No penalised Cox model, time-varying-effect refit, external cohort, subgroup fairness validation, competing-risk analysis or deployment evaluation is included. These are explicitly future work, not hidden results.
- CPU reproducibility is seed-controlled; compiler, BLAS and parallel arithmetic can cause small platform differences. Exact package versions are locked, but R itself and system dependencies must also be installed compatibly.
