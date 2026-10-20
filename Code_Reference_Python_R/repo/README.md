# Causal covariate selection: three-stage framework and OAENet

Python and R implementations of covariate-selection methods for causal inference, with code to
replicate the synthetic experiments of two papers:

* **Enhanced three-stage framework** (Yang & Noor-E-Alam, manuscript): exposure model (logistic
  regression or linear SVM) → penalty smoothing (tanh or sigmoid) → weighted elastic-net outcome
  model → adaptive elastic-net outcome model.
* **Outcome-adaptive elastic net, OAENet** (Islam, Shikalgar & Noor-E-Alam, 2025): outcome-model
  coefficients define adaptive weights for a penalised logistic exposure model.

> **Which code produced the published results?**
> The simulation results of the three-stage paper were produced with `FS_CI_Final.R`.
> `mode = "paper"` in `causal_fs.py` / `causal_fs.R` reproduces that code; the R version calls the
> same packages (glmnet, msaenet, e1071, MatchIt) and is the reference implementation. The Python
> version re-implements the same estimators and agrees with R in distribution (see *Validation*).
> `mode = "general"` (the default) is a general-purpose variant and does **not** reproduce the papers.

## Repository structure

```
Code_Reference_Python_R/
├── README.md
├── OAENet_Paper_Replication_Python_STABLE.ipynb   replication of the OAENet (IISE) simulations, Python
├── OAENet_paper_replication.R                     replication of the OAENet (IISE) simulations, R
└── repo/
    ├── R/
    │   ├── causal_fs.R                        three_stage(), oaenet(), att_matching(), detect_outcome_type()
    │   ├── simulate_demo.R                    three-stage paper simulations (Figures 3 and 4), paper mode
    │   ├── Scenario3_OAENet_vs_ThreeStage.R   stand-alone check: Scenario 3, coefficients x1 vs x2,
    │   │                                      code blocks copied verbatim from FS_CI_Final.R
    │   ├── simulation_results_paper_R.csv     output of simulate_demo.R (x2, N = 200, 30 seeds)
    │   └── scenario3_x1_vs_x2.csv             output of Scenario3_OAENet_vs_ThreeStage.R
    └── python/
        ├── causal_fs.py                       three_stage(), oaenet(), att_matching(), detect_outcome_type()
        ├── simulate_demo.py                   three-stage paper simulations, paper mode (default)
        ├── simulate_demo_general.py           earlier general-mode demonstration (continuous and binary Y)
        ├── ThreeStage_vs_OAENet_paper_mode.ipynb  notebook version of simulate_demo.py with figures
        │                                          and a cell-by-cell comparison with the R results
        ├── simulation_results_paper_py_mult2.csv  output of the notebook / simulate_demo.py (x2)
        └── requirements.txt
```

## Implementation modes

| | `mode = "paper"` | `mode = "iise"` (OAENet only) | `mode = "general"` (default) |
|---|---|---|---|
| Purpose | reproduce `FS_CI_Final.R` | the published OAENet | general-purpose use |
| Covariates | standardised with `scale()` | standardised | standardised internally |
| Three-stage Stage 1 | `glm` logit or linear SVM (cost 1) | – | same |
| Three-stage Stage 2 | `cv.glmnet(alpha = .5, relax = TRUE, gamma = c(0,.25,.5,.75,1))`, `lambda.1se` | – | elastic net, α = 0.5, 1-SE rule |
| Three-stage Stage 3 | `msaenet::aenet`, α grid 0.05–0.95, `lambda.1se` | – | adaptive elastic net, α = 0.5 |
| Three-stage selection | \|Θ\| ≥ 0.1 | – | Θ ≠ 0 |
| Outcome models | always Gaussian | Gaussian (OLS) | Gaussian or logistic (auto-detected) |
| OAENet | **ONet**, the OAENet implementation of `FS_CI_Final.R` (see below) | OLS Y ~ X, w = \|θ\|^(−3), logistic adaptive elastic net, 5-fold CV, 1-SE | w = \|θ\|^(−γ), (λ, γ) by wAMD |
| ATT | MatchIt nearest 1:1, caliper 0.25, `discard = "both"`; coefficient of A in `lm(Y ~ A + X)` | same as paper | binary Y → matched risk difference |

### OAENet implementations

The OAENet used in the three-stage paper (**ONet**, `mode = "paper"`) is an implementation of the
published OAENet (Islam et al., 2025) with the following modifications: the outcome model adjusts for
the treatment, the adaptive weights are transformed and truncated, and (λ, γ) are tuned jointly by
wAMD instead of cross-validation with a fixed γ. `mode = "iise"` follows the published description.

| | Published description (`mode = "iise"`) | ONet, as used in the three-stage paper (`mode = "paper"`) | `mode = "general"` |
|---|---|---|---|
| Outcome model | OLS: Y ~ X | lm: Y ~ A + X (adjusts for treatment) | Y ~ A + X (logistic for binary Y) |
| Penalty weights | \|θ\|^(−γ) | log, truncate at 0, normalise, sigmoid; values > 0.9 set to 1000 | \|θ\|^(−γ) |
| γ | fixed at 3 | tuned jointly with λ | tuned jointly with λ |
| Tuning | 5-fold CV, 1-SE | wAMD over λ = n^e | wAMD over λ = n^e |
| Propensity score in wAMD | – | linear predictor (`predict()` default); `ps_type = "response"` for probabilities | probabilities |
| Selection | coefficient ≠ 0 | normalised coefficient ≥ 0.01 | coefficient ≠ 0 |

The published description does not specify the elastic-net mixing parameter; α = 0.5 is used in
`mode = "iise"`.

## Usage

### Inputs

| Argument | Requirement |
|---|---|
| `X` | numeric matrix, n × p, no missing values; encode categorical covariates (e.g. one-hot) beforehand |
| `Y` | numeric outcome, continuous or binary 0/1 |
| `T` | treatment coded 0/1 |

Use `mode = "general"` (default) for applied work and `mode = "paper"` to reproduce the three-stage
paper. `selected` is returned as column indices of `X`: **0-based in Python, 1-based in R**.

### Task 1: covariate selection

**Python**
```python
from causal_fs import three_stage, oaenet

res = three_stage(X, Y, T, prior="SVM", smoother="sigmoid")   # Enh-ESVMS
# res = three_stage(X, Y, T, prior="LR", smoother="tanh")     # Enh-ELRT
sel = res["selected"]                  # indices of the selected columns
coef = res["Theta"]                    # Stage-3 coefficients

oae = oaenet(X, Y, T)                  # OAENet
sel_oae = oae["selected"]

selected_names = [colnames[j] for j in sel]   # colnames: list of column names of X
```

**R**
```r
source("causal_fs.R")

res <- three_stage(X, Y, T, prior = "SVM", smoother = "sigmoid")   # Enh-ESVMS
# res <- three_stage(X, Y, T, prior = "LR", smoother = "tanh")     # Enh-ELRT
sel  <- res$selected                   # indices of the selected columns
coef <- res$Theta                      # Stage-3 coefficients

oae <- oaenet(X, Y, T)                 # OAENet
sel_oae <- oae$selected

colnames(X)[sel]
```

Add `mode = "paper"` to either function to use the implementation of `FS_CI_Final.R`, or
`mode = "iise"` in `oaenet()` for the published description.

### Task 2: ATT estimation

`att_matching()` estimates the average treatment effect on the treated by **propensity score
matching (PSM)**, the only estimator implemented:

1. propensity score from a logistic regression of `T` on the covariates passed to the function;
2. common support: units outside the overlap of the treated and control propensity scores are dropped;
3. 1:1 nearest-neighbour matching without replacement, treated units in decreasing order of the
   propensity score, caliper 0.25 standard deviations of the propensity score; treated units without
   a control within the caliper are left unmatched;
4. effect on the matched sample:

| Outcome | `mode = "general"` | `mode = "paper"` |
|---|---|---|
| continuous | coefficient of `T` in `lm(Y ~ T + X)` | coefficient of `T` in `lm(Y ~ T + X)` |
| binary | difference in mean outcome between matched treated and controls (risk difference) | coefficient of `T` in `lm(Y ~ T + X)` |

Pass only the covariates to adjust for, normally the output of Task 1.

**Python** — pass the selected columns:
```python
from causal_fs import three_stage, att_matching

sel = three_stage(X, Y, T, prior="SVM", smoother="sigmoid")["selected"]
att = att_matching(X[:, sel], Y, T)              # ATT after covariate selection
att_all = att_matching(X, Y, T)                  # ATT adjusting for all covariates
```

**R** — pass the full matrix and the selected indices:
```r
sel <- three_stage(X, Y, T, prior = "SVM", smoother = "sigmoid")$selected
att <- att_matching(X, Y, T, selected = sel)     # ATT after covariate selection
att_all <- att_matching(X, Y, T)                 # ATT adjusting for all covariates
```

`att_matching()` returns a point estimate only. For a confidence interval, resample the data and
repeat both steps; note that the plain bootstrap is not valid for nearest-neighbour matching without
modification (Abadie & Imbens, 2008).

### Running the simulations

**R** (run from `repo/R/`)
```bash
Rscript simulate_demo.R
```
Settings at the top of `simulate_demo.R`: `COEF_MULT`, `SCENARIOS`, `RHO_LIST`, `ITR_MAX`, `N`, `P`,
`OUT_FILE`, `EXPORT_DIR`.

**Python** (run from `repo/python/`)
```bash
python simulate_demo.py --mode paper --coef-mult 2 --n 200 --reps 30   # three-stage paper simulations
python simulate_demo.py --mode paper --from-r ../R/r_export            # same-data comparison with R
python simulate_demo.py --mode general --reps 10 --n 500               # earlier general demonstration
```
or open `ThreeStage_vs_OAENet_paper_mode.ipynb`.

Packages — R: `glmnet`, `msaenet`, `e1071`, `MatchIt`, `MASS`. Python: see `requirements.txt`
(`numpy`, `pandas`, `scipy`, `scikit-learn`; `matplotlib` for the notebook).

## Simulation design (three-stage paper, paper mode)

* X ~ N(0, Σ), Σ_jk = ρ, p = 20, ρ ∈ {0, 0.25, 0.5, 0.75}; X is standardised after generation.
* A = 1{U < expit(Xα)}, U ~ Uniform(0, 1); Y = Xβ + ε + 0.5A, ε ~ N(0, 1).
* X1–X2 confounders, X3–X4 pure outcome predictors, X5–X6 pure treatment predictors, X7–X20 noise.

| Scenario | outcome coefficients β (X1–X4) | treatment coefficients α (X1, X2, X5, X6) |
|---|---|---|
| 1 | 0.6, 0.6, 0.6, 0.6 | 1, 1, 1, 1 |
| 2 | 0.6, 0.6, 0.6, 0.6 | 0.4, 0.4, 1, 1 |
| 3 | 0.2, 0.2, 0.6, 0.6 | 1, 1, 1, 1 |
| 4 | 0.6, 0.6, 0.6, 0.6 | 1, 1, 1.8, 1.8 |

**Coefficient multiplier.** `FS_CI_Final.R` multiplies all coefficients by 2
(`beta_v = 2*c(...)`, `alpha_v = 2*c(...)`); the published figures correspond to `COEF_MULT = 2`.
`COEF_MULT = 1` gives the values above as printed in Section 6.1 of the paper.

**Metric.** Percentage selection bias = |ATT(selected) − ATT(target)| / 0.5 × 100, where the target
model matches on X1–X4; reported with 90% confidence intervals of the mean, pooled over ρ
(Figure 3), together with the probability that each covariate is selected (Figure 4).

## Validation

**R paper mode vs. the published figures** (`simulation_results_paper_R.csv`, x2, N = 200, 30 seeds,
pooled over ρ):

| Scenario | Enh-ELRT | Enh-ESVMS | ONet |
|---|---|---|---|
| 1 | 0.8 [0.0, 2.1] | 1.2 [0.0, 2.7] | 8.6 [5.6, 11.7] |
| 2 | 2.4 [0.6, 4.2] | 1.6 [0.1, 3.1] | 7.5 [5.0, 10.0] |
| 3 | 13.1 [9.0, 17.3] | 10.9 [7.1, 14.7] | 17.4 [13.3, 21.6] |
| 4 | 1.4 [0.0, 2.9] | 1.4 [0.0, 2.9] | 7.1 [3.8, 10.3] |

The three-stage results agree with Figures 3 and 4 in all scenarios; ONet agrees in Scenarios 1–3.

**Python vs. R paper mode.** In 47 of 48 (scenario × ρ × method) cells the distributions of the
selection bias do not differ significantly (Mann–Whitney, p ≥ 0.05). ONet agrees closely. The Python
three-stage implementation selects pure treatment predictors and noise variables slightly less often
than R, which gives somewhat lower bias in Scenario 3; for reporting, use the R results.

**Coefficient multiplier** (`scenario3_x1_vs_x2.csv`, Scenario 3): with `COEF_MULT = 1`, confounder
selection by the three-stage framework drops to about 0.4–0.5 and its bias rises to 24–28%; it
remains lower than that of ONet (34%).

## Reproducibility notes

* Cross-validation folds are random in both languages, and R and Python use different random-number
  generators, so results agree in distribution, not data set by data set. For a same-data comparison,
  set `EXPORT_DIR` in `simulate_demo.R` and run `simulate_demo.py --from-r`.
* `numpy.random.Generator.multivariate_normal` can give different draws for correlated designs
  (ρ > 0) on different platforms or numpy versions; Python results for ρ > 0 may therefore not be
  identical across machines.
* Python paper mode re-implements glmnet behaviour: penalty-factor rescaling, 1/n standardisation,
  Gaussian y standardisation, λ sequence and early stopping, interpolation in `s`, grouped CV and the
  1-SE rule. Three details follow the logic of the R packages rather than a line-by-line port: the
  relaxed `gamma.1se` rule (largest λ, then largest γ, within 1 SE), msaenet's choice of α by
  min √CV-error, and the partial matching of `penalty.factor` to `penalty.factor.init` in `aenet()`.
* Paper mode reproduces two behaviours of `FS_CI_Final.R`: the ATT is `NA` when exactly one covariate
  is selected (the column name is lost in `xx[, idx]`), and ONet uses the linear predictor as the
  propensity score in wAMD. Neither changes the conclusions in the simulations above.

## OAENet (IISE) replication

`OAENet_Paper_Replication_Python_STABLE.ipynb` and `OAENet_paper_replication.R` re-run the
simulations of Islam et al. (2025) following the description in that paper (OLS outcome model,
γ = 3, 5-fold CV with the 1-SE rule). In Scenario 2 the selection of confounders and outcome
predictors is reproduced. In Scenario 1A the pure outcome predictors X11–X20 are selected in about
8% of runs, against about 90% reported in the paper. The data-generating process as printed has no
outcome noise, so the variance of the ATT is not comparable with the reported values.

## General mode

`mode = "general"` (the default) is intended for applied work:

* Outcome type is detected automatically (two distinct values → binary, logistic outcome models;
  otherwise Gaussian). Override with `outcome_type`. The treatment must be coded 0/1.
* glmnet-style fits with a data-driven λ path and the 1-SE rule; α = 0.5 in Stages 2 and 3;
  γ1 = 1 (sigmoid) or 0.5 (tanh); γ2 = 1; Θ ≠ 0 selects a covariate.
* OAENet with w = |θ|^(−γ) and (λ, γ) chosen by minimum wAMD; ridge term λ₂ = λ.

## References

* Islam MS, Shikalgar S, Noor-E-Alam M. A two-stage feature selection approach for robust evaluation
  of treatment effects in high-dimensional observational data. *IISE Transactions on Healthcare
  Systems Engineering* 15(2):117–136, 2025. https://doi.org/10.1080/24725579.2024.2447715
* Shortreed SM, Ertefaie A. Outcome-adaptive lasso: variable selection for causal inference.
  *Biometrics* 73(4):1111–1122, 2017.
* Yang T, Noor-E-Alam M. Optimizing feature selection in causal inference: a three-stage
  computational framework for unbiased estimation (manuscript).
* Abadie A, Imbens GW. On the failure of the bootstrap for matching estimators. *Econometrica*
  76(6):1537–1557, 2008.
* Zou H, Zhang HH. On the adaptive elastic-net with a diverging number of parameters.
  *Annals of Statistics* 37(4):1733–1751, 2009.
