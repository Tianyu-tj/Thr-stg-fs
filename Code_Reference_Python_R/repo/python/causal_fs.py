"""
causal_fs.py
============
Covariate selection for causal inference:

* Enhanced three-stage framework (Yang & Noor-E-Alam): exposure model -> penalty
  smoothing (sigmoid / tanh) -> weighted elastic-net outcome model -> adaptive
  elastic-net outcome model.
* Outcome-adaptive elastic net, OAENet (Islam, Shikalgar & Noor-E-Alam, 2024),
  tuned by the weighted absolute mean difference (wAMD; Shortreed & Ertefaie, 2017).

The outcome type is detected automatically:
    binary outcome     -> logistic (binomial) outcome models
    continuous outcome -> linear (Gaussian) outcome models
It can also be set explicitly with ``outcome_type="binary" | "continuous"``.

Penalised fits follow glmnet conventions:
    objective = loss/n + lam * sum_j pf_j * [alpha*|b_j| + (1-alpha)/2 * b_j^2]
with covariates standardised internally, a data-driven lambda path and
lambda selected by the 1-SE rule of K-fold cross-validation
(mean squared error for Gaussian, binomial deviance for binomial).

Dependencies: numpy, scipy, scikit-learn.
"""
import warnings
import numpy as np
from scipy.special import expit
from sklearn.linear_model import Lasso, LogisticRegression
from sklearn.svm import SVC
from sklearn.model_selection import KFold, StratifiedKFold

__all__ = ["detect_outcome_type", "three_stage", "oaenet", "att_matching", "CONFIG"]

CONFIG = dict(
    nlambda=60,
    lambda_min_ratio=1e-4,
    nfolds=5,
    alpha_stage2=0.5,     # elastic-net mixing, Stage 2 (Eq. 12)
    alpha_stage3=0.5,     # elastic-net mixing, Stage 3 (Eq. 13)
    gamma1_sigmoid=1.0,   # Online Supplement K.2
    gamma1_tanh=0.5,
    gamma2=1.0,           # Online Supplement K.3
    rule="1se",           # "1se" or "min"
    huge=1e10,            # penalty factors above this are treated as exclusions
)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def detect_outcome_type(y, outcome_type="auto", verbose=True):
    """Return ("binary" | "continuous", y as float array).

    binary: exactly two distinct non-missing values (recoded to 0/1 if needed).
    continuous: more than two distinct values.
    """
    y = np.asarray(y, dtype=float)
    if np.isnan(y).any():
        raise ValueError("Outcome contains missing values; handle them before selection.")
    vals = np.unique(y)
    if outcome_type == "auto":
        if len(vals) < 2:
            raise ValueError("Outcome is constant.")
        outcome_type = "binary" if len(vals) == 2 else "continuous"
        if outcome_type == "continuous" and len(vals) <= 5 and verbose:
            warnings.warn(f"Outcome has only {len(vals)} distinct values; treated as continuous. "
                          "Set outcome_type explicitly if this is not intended.")
    if outcome_type == "binary":
        if len(vals) != 2:
            raise ValueError(f"outcome_type='binary' but the outcome has {len(vals)} distinct values.")
        if not np.array_equal(vals, [0.0, 1.0]):
            if verbose:
                warnings.warn(f"Binary outcome recoded: {vals[0]} -> 0, {vals[1]} -> 1.")
            y = (y == vals[1]).astype(float)
    elif outcome_type != "continuous":
        raise ValueError("outcome_type must be 'auto', 'binary' or 'continuous'.")
    if verbose:
        print(f"[causal_fs] outcome type: {outcome_type} -> "
              f"{'logistic' if outcome_type == 'binary' else 'linear'} outcome models")
    return outcome_type, y


def _check_treatment(t):
    t = np.asarray(t)
    vals = np.unique(t)
    if not np.array_equal(vals, [0, 1]):
        raise ValueError("Treatment must be coded 0/1.")
    return t.astype(int)


def _unpenalised_logit():
    """Unpenalised logistic regression across scikit-learn versions."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        try:
            return LogisticRegression(penalty=None, max_iter=10000)
        except Exception:                       # very old scikit-learn
            return LogisticRegression(C=1e10, max_iter=10000)


def _fit_logit(X, t):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        return _unpenalised_logit().fit(X, t)


def _standardize(X):
    X = np.asarray(X, dtype=float)
    mu, sd = X.mean(0), X.std(0)
    sd = np.where(sd > 1e-12, sd, 1.0)
    return (X - mu) / sd, mu, sd


# ---------------------------------------------------------------------------
# Penalised regression (glmnet conventions)
# ---------------------------------------------------------------------------
def _gauss_fit(Xc, yc, pf, lam, alpha, warm=None):
    """(1/2n)||yc - Xc b||^2 + lam*sum pf_j [alpha|b_j| + (1-alpha)/2 b_j^2]; Xc, yc centred.
    Solved exactly by the Lasso on an augmented, column-rescaled design."""
    n, p = Xc.shape
    l2 = lam * (1 - alpha) * pf
    Xa = np.vstack([Xc, np.diag(np.sqrt(n * l2))]) if np.any(l2 > 0) else Xc
    ya = np.concatenate([yc, np.zeros(p)]) if np.any(l2 > 0) else yc
    na = Xa.shape[0]
    m = Lasso(alpha=max(lam * alpha * n / na, 1e-12), fit_intercept=False, max_iter=20000, tol=1e-7,
              warm_start=warm is not None)
    if warm is not None:
        m.coef_ = warm * pf
    m.fit(Xa / pf, ya)
    return m.coef_ / pf


def _binom_fit(Xs, y, pf, lam, alpha, b=None, b0=None, iters=25, tol=1e-6):
    """-loglik/n + penalty, by IRLS; each step is a weighted Gaussian elastic net."""
    n, p = Xs.shape
    b = np.zeros(p) if b is None else b.copy()
    b0 = np.log(y.mean() / (1 - y.mean())) if b0 is None else b0
    for _ in range(iters):
        eta = b0 + Xs @ b
        pr = np.clip(expit(eta), 1e-6, 1 - 1e-6)
        W = pr * (1 - pr)
        z = eta + (y - pr) / W
        sw = W / W.mean()
        xm = (sw[:, None] * Xs).mean(0); zm = (sw * z).mean()
        r = np.sqrt(sw)
        bn = _gauss_fit((Xs - xm) * r[:, None], (z - zm) * r, pf, lam / W.mean(), alpha, warm=b)
        b0n = zm - xm @ bn
        done = np.max(np.abs(bn - b)) < tol and abs(b0n - b0) < tol
        b, b0 = bn, b0n
        if done:
            break
    return b0, b


def _path(X, y, pf, alpha, family, lambdas=None):
    """Coefficient path on the original X scale. pf entries >= CONFIG['huge'] or inf are excluded."""
    n, p = X.shape
    Xs, mu, sd = _standardize(X)
    pf = np.asarray(pf, float)
    act = np.isfinite(pf) & (pf < CONFIG["huge"])
    pfa = np.maximum(pf[act], 1e-8)
    pfa = pfa * act.sum() / pfa.sum()            # glmnet: penalty factors rescaled to sum to #variables
    Xa = Xs[:, act]
    ybar = y.mean()
    if lambdas is None:
        g = np.abs(Xa.T @ (y - ybar)) / n
        lmax = np.max(g / (pfa * max(alpha, 1e-3))) if act.any() else 1.0
        lambdas = lmax * np.geomspace(1, CONFIG["lambda_min_ratio"], CONFIG["nlambda"])
    coefs = np.zeros((len(lambdas), p)); b0s = np.zeros(len(lambdas))
    b, b0 = None, None
    for k, lam in enumerate(lambdas):
        full = np.zeros(p)
        if act.any():
            if family == "gaussian":
                b = _gauss_fit(Xa - Xa.mean(0), y - ybar, pfa, lam, alpha, warm=b)
                b0 = ybar
            else:
                b0, b = _binom_fit(Xa, y, pfa, lam, alpha, b, b0)
            full[act] = b
        else:
            b0 = ybar if family == "gaussian" else np.log(ybar / (1 - ybar))
        beta = full / sd
        coefs[k] = beta
        b0s[k] = b0 - mu @ beta
    return np.asarray(lambdas), b0s, coefs


def cv_enet(X, y, pf, alpha, family, seed=1, rule=None):
    """cv.glmnet analogue; returns coefficients at lambda.1se (or lambda.min)."""
    rule = rule or CONFIG["rule"]
    lambdas, b0f, cf = _path(X, y, pf, alpha, family)
    K = CONFIG["nfolds"]
    splitter = (StratifiedKFold(K, shuffle=True, random_state=seed).split(X, y) if family == "binomial"
                else KFold(K, shuffle=True, random_state=seed).split(X))
    err = np.zeros((K, len(lambdas)))
    for f, (tr, va) in enumerate(splitter):
        _, b0, c = _path(X[tr], y[tr], pf, alpha, family, lambdas)
        for l in range(len(lambdas)):
            eta = b0[l] + X[va] @ c[l]
            if family == "gaussian":
                err[f, l] = np.mean((y[va] - eta) ** 2)
            else:
                pr = np.clip(expit(eta), 1e-10, 1 - 1e-10)
                err[f, l] = -2 * np.mean(y[va] * np.log(pr) + (1 - y[va]) * np.log(1 - pr))
    cvm = err.mean(0); cvsd = err.std(0, ddof=1) / np.sqrt(K)
    i = int(np.argmin(cvm))
    l = int(np.flatnonzero(cvm <= cvm[i] + cvsd[i])[0]) if rule == "1se" else i
    return dict(coef=cf[l], intercept=b0f[l], lam=lambdas[l], cvm=cvm[l])


# ---------------------------------------------------------------------------
# Enhanced three-stage framework (Algorithm 1)
# ---------------------------------------------------------------------------
def _exposure_coef(Xs, t, prior):
    if prior == "LR":      # glm(T ~ X, family = binomial)
        b = _fit_logit(Xs, t).coef_.ravel()
    elif prior == "SVM":   # e1071::svm(kernel = "linear", cost = 1): libsvm, hinge loss
        b = SVC(kernel="linear", C=1.0).fit(Xs, t).coef_.ravel()
    else:
        raise ValueError("prior must be 'LR' or 'SVM'")
    return np.where(np.isfinite(b), b, 0.0)


def _smooth(beta, smoother, gamma1):
    a = np.abs(beta)
    if smoother == "sigmoid":
        return (1.0 / (1.0 + np.exp(-a))) ** gamma1
    if smoother == "tanh":
        return np.maximum(np.tanh(a), 1e-8) ** gamma1
    raise ValueError("smoother must be 'sigmoid' or 'tanh'")


def three_stage(X, Y, T, prior="SVM", smoother="sigmoid", gamma1=None, gamma2=None,
                outcome_type="auto", seed=1, verbose=True):
    """Enhanced three-stage feature selection.

    prior/smoother: ("SVM","sigmoid") = Enh-ESVMS, ("LR","tanh") = Enh-ELRT, etc.
    Returns dict with 'selected' (column indices), 'Theta' (Stage-3 coefficients on the
    standardised scale) and intermediate quantities.
    """
    otype, Y = detect_outcome_type(Y, outcome_type, verbose)
    T = _check_treatment(T)
    family = "binomial" if otype == "binary" else "gaussian"
    if gamma1 is None:
        gamma1 = CONFIG["gamma1_sigmoid"] if smoother == "sigmoid" else CONFIG["gamma1_tanh"]
    gamma2 = CONFIG["gamma2"] if gamma2 is None else gamma2
    Xs, _, _ = _standardize(X)
    beta = _exposure_coef(Xs, T, prior)                                  # Stage 1
    w = _smooth(beta, smoother, gamma1)                                  # penalty smoothing
    s2 = cv_enet(Xs, Y, w, CONFIG["alpha_stage2"], family, seed)         # Stage 2, Eq. (12)
    theta = s2["coef"]
    with np.errstate(divide="ignore"):
        psi = np.where(np.abs(theta) > 0, (1.0 / np.abs(theta)) ** gamma2, np.inf)
    s3 = cv_enet(Xs, Y, psi, CONFIG["alpha_stage3"], family, seed + 1)  # Stage 3, Eq. (13)
    Theta = s3["coef"]
    return dict(selected=np.flatnonzero(np.abs(Theta) > 1e-10), Theta=Theta, theta_stage2=theta,
                beta_exposure=beta, w=w, psi=psi, outcome_type=otype)


# ---------------------------------------------------------------------------
# OAENet with wAMD tuning
# ---------------------------------------------------------------------------
def _logistic_enet_fista(X, t, l1w, l2, max_iter=5000, tol=1e-8):
    """min mean logistic loss + sum_j l1w_j|b_j| + l2/2 ||b||^2 (FISTA; intercept unpenalised)."""
    n, p = X.shape
    b0 = np.log(t.mean() / (1 - t.mean())); b = np.zeros(p)
    L = np.linalg.norm(X, 2) ** 2 / (4 * n) + l2 + 1e-8
    z0, z, tk = b0, b.copy(), 1.0
    for _ in range(max_iter):
        r = expit(z0 + X @ z) - t
        g0 = r.mean(); g = X.T @ r / n + l2 * z
        nb0 = z0 - g0 / L
        u = z - g / L
        nb = np.sign(u) * np.maximum(np.abs(u) - l1w / L, 0)
        d = max(abs(nb0 - b0), np.max(np.abs(nb - b)) if p else 0.0)
        tn = (1 + np.sqrt(1 + 4 * tk * tk)) / 2
        z0 = nb0 + (tk - 1) / tn * (nb0 - b0); z = nb + (tk - 1) / tn * (nb - b)
        b0, b, tk = nb0, nb, tn
        if d < tol:
            break
    return b0, b


def _wamd(Xs, t, ps, theta):
    w = np.where(t == 1, 1 / ps, 1 / (1 - ps))
    m1 = (w[t == 1, None] * Xs[t == 1]).sum(0) / w[t == 1].sum()
    m0 = (w[t == 0, None] * Xs[t == 0]).sum(0) / w[t == 0].sum()
    return float(np.sum(np.abs(theta) * np.abs(m1 - m0)))


def oaenet(X, Y, T, outcome_type="auto", lam_exps=(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
           conv=2.0, ridge_ratio=1.0, verbose=True):
    """OAENet: outcome model Y ~ T + X (linear or logistic) gives theta; exposure model is a
    logistic elastic net with L1 weights |theta|^-gamma; (lambda, gamma) chosen by minimum wAMD.
    lambda = n^lam_exp on the summed-loss scale, gamma = 2*(conv - lam_exp + 1); L2 = ridge_ratio*lambda."""
    otype, Y = detect_outcome_type(Y, outcome_type, verbose)
    T = _check_treatment(T)
    n, p = np.asarray(X).shape
    Xs, _, _ = _standardize(X)
    D = np.column_stack([T, Xs])
    if otype == "binary":
        theta = _fit_logit(D, Y).coef_[0][1:]
    else:
        theta = np.linalg.lstsq(np.column_stack([np.ones(n), D]), Y, rcond=None)[0][2:]
    best = None
    for le in lam_exps:
        lam_sum = float(n) ** le
        gam = 2 * (conv - le + 1)
        logw = np.log(lam_sum / n) - gam * np.log(np.maximum(np.abs(theta), 1e-300))
        keep = logw < np.log(0.6)   # exact screen: |gradient| <= 0.5 for standardised X, so larger penalties give 0
        b = np.zeros(p)
        if keep.any():
            _, bk = _logistic_enet_fista(Xs[:, keep], T, np.exp(logw[keep]), ridge_ratio * lam_sum / n)
            b[keep] = bk
        eta = Xs @ b
        b0 = np.log(T.mean() / (1 - T.mean()))
        for _ in range(50):
            pr = expit(b0 + eta)
            b0 -= np.mean(pr - T) / (np.mean(pr * (1 - pr)) + 1e-12)
        ps = np.clip(expit(b0 + eta), 1e-6, 1 - 1e-6)
        val = _wamd(Xs, T, ps, theta)
        if best is None or val < best["wAMD"]:
            best = dict(wAMD=val, lam_exp=le, gamma=gam, beta=b.copy())
    best.update(selected=np.flatnonzero(np.abs(best["beta"]) > 1e-10), theta_outcome=theta, outcome_type=otype)
    return best


# ---------------------------------------------------------------------------
# ATT by propensity-score matching on the selected covariates
# ---------------------------------------------------------------------------
def att_matching(X, Y, T, outcome_type="auto", caliper=0.25, verbose=False):
    """1:1 nearest-neighbour matching on the logistic propensity score (largest PS first, no replacement,
    caliper in SD of the PS, common support enforced). Continuous outcome: coefficient of T in
    lm(Y ~ T + X) on the matched sample. Binary outcome: matched-pair risk difference."""
    otype, Y = detect_outcome_type(Y, outcome_type, verbose)
    T = _check_treatment(T)
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    if X.shape[1] == 0:
        return float(Y[T == 1].mean() - Y[T == 0].mean())
    Xs, _, _ = _standardize(X)
    ps = _fit_logit(Xs, T).predict_proba(Xs)[:, 1]
    lo, hi = max(ps[T == 1].min(), ps[T == 0].min()), min(ps[T == 1].max(), ps[T == 0].max())
    keep = (ps >= lo) & (ps <= hi)
    width = caliper * ps.std(ddof=1)
    tr = np.flatnonzero((T == 1) & keep); tr = tr[np.argsort(-ps[tr])]
    ctrl = np.flatnonzero((T == 0) & keep); used = np.zeros(len(T), bool); pairs = []
    order = ctrl[np.argsort(ps[ctrl])]; cps = ps[order]
    for i in tr:
        pos = np.searchsorted(cps, ps[i]); best, bd = -1, np.inf
        for direction in (-1, 1):                       # nearest unused control on each side
            j = pos if direction == 1 else pos - 1
            while 0 <= j < len(order) and used[order[j]]:
                j += direction
            if 0 <= j < len(order) and abs(cps[j] - ps[i]) < bd:
                best, bd = order[j], abs(cps[j] - ps[i])
        if best >= 0 and bd <= width:
            pairs.append((i, best)); used[best] = True
    if not pairs:
        return float("nan")
    pairs = np.array(pairs)
    if otype == "binary":
        return float(np.mean(Y[pairs[:, 0]] - Y[pairs[:, 1]]))
    idx = pairs.ravel()
    D = np.column_stack([np.ones(len(idx)), T[idx], X[idx]])
    return float(np.linalg.lstsq(D, Y[idx], rcond=None)[0][1])
