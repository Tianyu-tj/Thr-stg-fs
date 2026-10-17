# Causal covariate selection: three-stage framework and OAENet

Python and R implementations of two covariate-selection methods for causal inference, with a
synthetic-data demonstration.

* **Enhanced three-stage framework** (Yang & Noor-E-Alam): exposure model (logistic regression or
  linear SVM) → penalty smoothing (sigmoid or tanh) → weighted elastic-net outcome model →
  adaptive elastic-net outcome model. Variables with non-zero Stage-3 coefficients are selected.
* **Outcome-adaptive elastic net, OAENet** (Islam, Shikalgar & Noor-E-Alam, 2024): outcome-model
  coefficients define adaptive weights for a penalised logistic exposure model; (λ, γ) are chosen by
  the weighted absolute mean difference, wAMD (Shortreed & Ertefaie, 2017).

## Outcome type is detected automatically

| Outcome | Detected when | Outcome models used |
|---|---|---|
| binary | exactly two distinct values (recoded to 0/1 if needed) | logistic / binomial elastic net |
| continuous | more than two distinct values | linear / Gaussian elastic net |

Override with `outcome_type="binary"` or `outcome_type="continuous"`. A warning is issued when an
outcome with five or fewer distinct values is treated as continuous. The treatment must be coded 0/1.

## Files

```
python/
  causal_fs.py        three_stage(), oaenet(), att_matching(), detect_outcome_type()
  simulate_demo.py    synthetic experiments (continuous and binary outcomes)
  requirements.txt
R/
  causal_fs.R         three_stage(), oaenet(), att_matching(), detect_outcome_type()
  simulate_demo.R     the same experiments in R  (packages: glmnet, e1071, MatchIt, MASS)
```

## Usage

**Python**
```python
from causal_fs import three_stage, oaenet, att_matching
res = three_stage(X, Y, T, prior="SVM", smoother="sigmoid")   # Enh-ESVMS; ("LR", "tanh") = Enh-ELRT
sel = res["selected"]                                         # column indices of X
res_oae = oaenet(X, Y, T)
att = att_matching(X[:, sel], Y, T)
```
```bash
cd python && python simulate_demo.py --reps 10 --n 500
```

**R**
```r
source("causal_fs.R")
res <- three_stage(X, Y, T, prior = "SVM", smoother = "sigmoid")
res_oae <- oaenet(X, Y, T)
att <- att_matching(X[, res$selected, drop = FALSE], Y, T)
```
```bash
cd R && Rscript simulate_demo.R
```

## Simulation design

`X ~ N(0, Σ)` with `Σ_jk = ρ`, `p = 20`; `T ~ Bernoulli(expit(Xβ))`.
Continuous outcome: `Y = 0.5T + Xθ + ε`. Binary outcome: `Y ~ Bernoulli(expit(−1 + 0.5T + Xθ))`.
X1–X2 are confounders, X3–X4 pure outcome predictors, X5–X6 pure treatment predictors, X7–X20 noise;
the desired selection is {X1, X2, X3, X4}. Scenarios 1–4 follow Section 6.1 of the three-stage paper.
Reported: recovery of {X1–X4}, selection of X5/X6, number of noise variables, and absolute ATT bias
relative to the target model (matching on X1–X4) and to the true ATT.

## Implementation notes

* Penalised fits follow glmnet conventions: covariates standardised internally, data-driven λ path,
  penalty `λ Σ pf_j [α|b_j| + (1−α)/2 b_j²]` with penalty factors rescaled to sum to the number of
  variables, λ chosen by the 1-SE rule of 5-fold cross-validation (MSE for Gaussian, binomial deviance
  for binomial). Variables with zero Stage-2 coefficient receive an infinite Stage-3 penalty (excluded).
* Defaults: α = 0.5 in Stages 2 and 3; γ1 = 1 (sigmoid) or 0.5 (tanh); γ2 = 1.
* OAENet grid: λ = n^e for e ∈ {−10, −5, −2, −1, −0.75, −0.5, −0.25, 0.25, 0.49} (summed-loss scale),
  γ = 2(2 − e + 1), ridge term λ₂ = λ; the pair with minimum wAMD is retained.
* Linear SVM uses the hinge loss with cost 1 (libsvm: `sklearn.svm.SVC(kernel="linear")`,
  `e1071::svm(kernel = "linear")`).
* `att_matching()`: 1:1 nearest-neighbour matching on the logistic propensity score, caliper 0.25 SD,
  common support enforced; continuous outcome → coefficient of T in `lm(Y ~ T + X)` on the matched
  sample; binary outcome → matched-sample risk difference.
* The Python and R versions implement the same estimators; numerical results can differ slightly
  because of solver tolerances, λ paths and random fold assignment.

## Changes relative to earlier versions of this code

1. **Outcome-type detection.** The papers write the outcome stages as least-squares models. For a
   binary outcome, logistic (binomial) outcome models are now used in Stages 2–3 of the three-stage
   framework and for the OAENet weights and wAMD.
2. **OAENet tuning** by minimum wAMD over (λ, γ), instead of cross-validated exposure deviance with a
   fixed γ.
3. **Linear SVM** with the hinge loss (libsvm), matching Eq. (6), instead of the squared hinge loss.
4. **Data-driven λ path** with the 1-SE rule, instead of a fixed λ grid.
5. **Selection rule**: non-zero Stage-3 coefficients; no fixed coefficient-size threshold (such a
   threshold depends on the outcome scale).
6. **Stage 3** is fitted directly with `glmnet` (R) using penalty factors and exclusions.
7. **Simulation**: treatment assigned as `Bernoulli(expit(Xβ))`.

## References

* Islam MS, Shikalgar S, Noor-E-Alam M. A two-stage feature selection approach for robust evaluation of
  treatment effects in high-dimensional observational data. *IISE Transactions on Healthcare Systems
  Engineering*, 2024.
* Shortreed SM, Ertefaie A. Outcome-adaptive lasso: variable selection for causal inference.
  *Biometrics* 73(4):1111–1122, 2017.
* Yang T, Noor-E-Alam M. Optimizing feature selection in causal inference: a three-stage computational
  framework for unbiased estimation (manuscript).
* Zou H, Zhang HH. On the adaptive elastic-net with a diverging number of parameters.
  *Annals of Statistics* 37(4):1733–1751, 2009.
