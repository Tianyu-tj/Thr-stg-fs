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


def _three_stage_general(X, Y, T, prior="SVM", smoother="sigmoid", gamma1=None, gamma2=None,
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


def _oaenet_general(X, Y, T, outcome_type="auto", lam_exps=(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
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
def _att_matching_general(X, Y, T, outcome_type="auto", caliper=0.25, verbose=False):
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


# ===========================================================================
# PAPER MODE: replication of FS_CI_Final.R
# ===========================================================================
# The functions below reproduce the R code used for the three-stage paper:
#   Stage 2  cv.glmnet(alpha = .5, nfolds = 5, relax = TRUE, gamma = c(0,.25,.5,.75,1)),
#            coef(fit, s = lambda.1se) with gamma = gamma.1se
#   Stage 3  msaenet::aenet(alphas = seq(.05,.95,.05), nfolds = 5, rule = "lambda.1se",
#            penalty.factor.init = Psi)   [R partial matching of `penalty.factor`]
#   select   |Theta| >= 0.1;  outcome models are always Gaussian (as in the R code)
#   ONet     the modified OAENet benchmark of FS_CI_Final.R (sigmoid-transformed log weights,
#            glmnet logistic path evaluated at s = n^lambda, wAMD, |normalised coef| >= 0.01)
#   ATT      MatchIt nearest 1:1, caliper .25, discard = "both"; coefficient of A in lm(Y ~ A + X)
# glmnet details reproduced: penalty factors rescaled to sum to nvars (Inf -> excluded, counted as 1),
# 1/n standardisation, Gaussian y standardisation, lambda.min.ratio, path early stopping
# (fdev = 1e-5, devmax = 0.999, mnlam = 5) on automatic paths only, linear interpolation in s,
# grouped CV with fold-size weights and the 1-SE rule; relaxed fits are unpenalised refits on the
# active set; gamma.1se = largest lambda, then largest gamma, within 1 SE of the minimum.
# Not reproducible across languages: R's random fold assignment (results agree in distribution).

GLMNET_CTRL = dict(nlambda=100, fdev=1e-5, devmax=0.999, mnlam=5, thresh=1e-7, maxit=10000)
_EPS = np.finfo(float).eps


def _std_r(X):
    """R's scale(): centre and divide by the (n-1) standard deviation."""
    X = np.asarray(X, float)
    sd = X.std(0, ddof=1)
    return (X - X.mean(0)) / np.where(sd > 0, sd, 1.0)


def _vp_glmnet(pf):
    pf = np.asarray(pf, float).copy()
    excl = ~np.isfinite(pf)
    pf[excl] = 1.0                                   # glmnet: excluded variables count as 1
    pf = np.maximum(pf, 0.0)
    return pf * len(pf) / pf.sum(), ~excl


def _cd(G, c, v, alpha, lam, b, thresh, maxit):
    """Covariance coordinate descent for (1/2) b'Gb - c'b + lam*sum v_j[alpha|b_j| + (1-alpha)/2 b_j^2]."""
    k = len(c)
    if k == 0:
        return b
    grad = c - G @ b
    Gd = np.diag(G)
    l1, l2 = lam * alpha * v, lam * (1 - alpha) * v
    cols = [G[:, j] for j in range(k)]
    for _ in range(maxit):
        dmax = 0.0
        for j in range(k):
            bj = b[j]; r = grad[j] + Gd[j] * bj
            a = abs(r) - l1[j]
            nb = (a / (Gd[j] + l2[j]) if a > 0 else 0.0) * (1.0 if r > 0 else -1.0)
            d = nb - bj
            if d != 0.0:
                b[j] = nb; grad -= cols[j] * d
                dd = Gd[j] * d * d
                if dd > dmax: dmax = dd
        if dmax < thresh:
            break
    return b


def glmnet_path(X, y, pf, alpha, family="gaussian", lambdas=None):
    """glmnet-compatible elastic-net path. Returns dict(lambda, a0, beta (nlam x p), dev_ratio)."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n, p = X.shape
    vp, keep = _vp_glmnet(pf)
    xm = X.mean(0); xs = X.std(0)
    keep &= xs > 0
    idx = np.flatnonzero(keep)
    Z = (X[:, idx] - xm[idx]) / xs[idx]
    v = vp[idx]
    user = lambdas is not None
    ctrl = GLMNET_CTRL
    if family == "gaussian":
        ym = y.mean(); ys = y.std(); ys = ys if ys > 0 else 1.0
        yc = (y - ym) / ys
        G = Z.T @ Z / n; c = Z.T @ yc / n
        grad0 = np.abs(c)
    else:
        ybar = y.mean()
        grad0 = np.abs(Z.T @ (y - ybar)) / n
        nulldev = -2 * n * (ybar * np.log(ybar) + (1 - ybar) * np.log(1 - ybar))
    pen = v > 0
    if user:
        lam = np.asarray(lambdas, float) / (ys if family == "gaussian" else 1.0)
    else:
        lmax = np.max(grad0[pen] / v[pen]) / max(alpha, 1e-3) if pen.any() else 1.0
        ratio = 1e-4 if n > p else 1e-2
        lam = lmax * ratio ** (np.arange(ctrl["nlambda"]) / (ctrl["nlambda"] - 1))
    k = len(idx); b = np.zeros(k)
    a0s, B, devs, lams = [], [], [], []
    b0 = np.log(ybar / (1 - ybar)) if family == "binomial" else 0.0
    rsq_prev = 0.0
    for m, lm in enumerate(lam):
        if family == "gaussian":
            b = _cd(G, c, v, alpha, lm, b, ctrl["thresh"], ctrl["maxit"])
            rsq = 2 * c @ b - b @ G @ b
            a0 = ym - ys * (b / xs[idx]) @ xm[idx]
            beta_k = ys * b / xs[idx]
        else:
            for _ in range(100):                                    # IRLS outer loop
                eta = b0 + Z @ b
                pr = np.clip(expit(eta), 1e-9, 1 - 1e-9)
                w = pr * (1 - pr)
                zr = eta + (y - pr) / w
                sw = w.sum()
                xw = (w @ Z) / sw; zw = (w @ zr) / sw
                Zc = Z - xw
                Gw = (Zc * w[:, None]).T @ Zc / n
                cw = Zc.T @ (w * (zr - zw)) / n
                bold = b.copy(); b0old = b0
                b = _cd(Gw, cw, v, alpha, lm, b, ctrl["thresh"], ctrl["maxit"])
                b0 = zw - xw @ b
                if np.max(np.abs(np.r_[b - bold, b0 - b0old])) < 1e-6:
                    break
            pr = np.clip(expit(b0 + Z @ b), 1e-15, 1 - 1e-15)
            dev = -2 * np.sum(y * np.log(pr) + (1 - y) * np.log(1 - pr))
            rsq = 1 - dev / nulldev
            beta_k = b / xs[idx]
            a0 = b0 - beta_k @ xm[idx]
        full = np.zeros(p); full[idx] = beta_k
        a0s.append(a0); B.append(full); devs.append(rsq); lams.append(lm)
        if not user and m + 1 >= ctrl["mnlam"]:
            if rsq - rsq_prev < ctrl["fdev"] * rsq or rsq > ctrl["devmax"]:
                break
        rsq_prev = rsq
    lams = np.array(lams) * (ys if family == "gaussian" else 1.0)
    return dict(lambda_=lams, a0=np.array(a0s), beta=np.array(B), dev_ratio=np.array(devs))


def _interp(path, s):
    """predict.glmnet / coef.glmnet at s: linear interpolation in lambda, clamped to the path."""
    lam = path["lambda_"]
    if s >= lam[0]:
        return path["a0"][0], path["beta"][0]
    if s <= lam[-1]:
        return path["a0"][-1], path["beta"][-1]
    k = np.flatnonzero(lam >= s)[-1]
    f = (s - lam[k + 1]) / (lam[k] - lam[k + 1])
    return f * path["a0"][k] + (1 - f) * path["a0"][k + 1], f * path["beta"][k] + (1 - f) * path["beta"][k + 1]


def _relaxed(X, y, beta_row):
    """Unpenalised Gaussian refit on the active set (glmnet relax)."""
    a = np.flatnonzero(beta_row != 0)
    out = np.zeros(X.shape[1])
    D = np.column_stack([np.ones(len(y)), X[:, a]])
    coef = np.linalg.lstsq(D, y, rcond=None)[0]
    out[a] = coef[1:]
    return coef[0], out


def cv_glmnet_r(X, y, pf, alpha, nfolds=5, seed=1, relax=False, gammas=(0, .25, .5, .75, 1)):
    """Gaussian cv.glmnet (type.measure = 'mse'), optionally relax = TRUE. Returns full-data fit,
    lambda.min / lambda.1se and, if relax, gamma.1se and the coef(s = lambda.1se) vector."""
    X = np.asarray(X, float); y = np.asarray(y, float); n = len(y)
    full = glmnet_path(X, y, pf, alpha, "gaussian")
    lam = full["lambda_"]
    rng = np.random.default_rng(seed)
    foldid = rng.permutation(np.resize(np.arange(nfolds), n))
    gs = np.asarray(gammas, float) if relax else np.array([1.0])
    cvraw = np.zeros((len(gs), nfolds, len(lam))); wts = np.zeros(nfolds)
    for f in range(nfolds):
        tr, va = foldid != f, foldid == f
        fp = glmnet_path(X[tr], y[tr], pf, alpha, "gaussian", lambdas=lam)
        pred_pen = fp["a0"][None, :] + X[va] @ fp["beta"].T
        if relax:
            rel = [_relaxed(X[tr], y[tr], fp["beta"][k]) for k in range(len(lam))]
            pred_rel = np.array([r[0] for r in rel])[None, :] + X[va] @ np.array([r[1] for r in rel]).T
        for gi, g in enumerate(gs):
            pred = pred_pen if g == 1 else g * pred_pen + (1 - g) * pred_rel
            cvraw[gi, f] = np.mean((y[va][:, None] - pred) ** 2, axis=0)
        wts[f] = va.sum()
    cvm = np.einsum("gfl,f->gl", cvraw, wts) / wts.sum()
    cvsd = np.sqrt(np.einsum("gfl,f->gl", (cvraw - cvm[:, None, :]) ** 2, wts) / wts.sum() / (nfolds - 1))
    g1 = len(gs) - 1                                             # gamma = 1 (unrelaxed)
    i_min = np.flatnonzero(cvm[g1] <= cvm[g1].min())[0]
    i_1se = np.flatnonzero(cvm[g1] <= cvm[g1, i_min] + cvsd[g1, i_min])[0]
    out = dict(path=full, lambda_min=lam[i_min], lambda_1se=lam[i_1se], i_min=i_min, i_1se=i_1se,
               cvm=cvm[g1], cvsd=cvsd[g1])
    if relax:                                                    # getOptcv.relaxed
        flat_cvm = cvm.ravel(); flat_sd = cvsd.ravel()
        lam_rep = np.tile(lam, len(gs)); gam_rep = np.repeat(gs, len(lam))
        idmin = flat_cvm <= flat_cvm.min()
        order = np.lexsort((-gam_rep[idmin], -lam_rep[idmin]))[0]
        imin = np.flatnonzero(idmin)[order]
        ok = flat_cvm <= flat_cvm[imin] + flat_sd[imin]
        order = np.lexsort((-gam_rep[ok], -lam_rep[ok]))[0]
        g1se = gam_rep[np.flatnonzero(ok)[order]]
        a_pen, b_pen = full["a0"][i_1se], full["beta"][i_1se]
        a_rel, b_rel = _relaxed(X, y, b_pen)
        out.update(gamma_1se=g1se, coef=g1se * b_pen + (1 - g1se) * b_rel,
                   intercept=g1se * a_pen + (1 - g1se) * a_rel)
    else:
        out.update(coef=full["beta"][i_1se], intercept=full["a0"][i_1se])
    return out


def aenet_r(X, y, penalty_factor_init, alphas=np.round(np.arange(.05, .951, .05), 2), nfolds=5,
            rule="lambda.1se", seed=1001):
    """msaenet::aenet (tune = 'cv', init = 'enet', scale = 1), Gaussian family."""
    def tune(pf):
        best = None
        for a in alphas:
            r = cv_glmnet_r(X, y, pf, a, nfolds, seed)
            err = np.sqrt(r["cvm"]).min()                 # msaenet: min(sqrt(cvm)) selects alpha
            if best is None or err < best[0]:
                best = (err, a, r)
        _, a, r = best
        lam = r["lambda_1se"] if rule == "lambda.1se" else r["lambda_min"]
        return a, lam
    a1, l1 = tune(penalty_factor_init)
    bhat = glmnet_path(X, y, penalty_factor_init, a1, "gaussian", lambdas=[l1])["beta"][0]
    if np.all(bhat == 0):
        bhat = np.full(len(bhat), 2 * _EPS)
    adpen = np.maximum(np.abs(bhat), _EPS) ** (-1.0)
    a2, l2 = tune(adpen)
    beta = glmnet_path(X, y, adpen, a2, "gaussian", lambdas=[l2])["beta"][0]
    return dict(beta=beta, alpha_enet=a1, lambda_enet=l1, alpha_aenet=a2, lambda_aenet=l2)


def _exposure_coef_r(Xs, t, prior):
    if prior == "LR":                                     # glm(A ~ ., binomial)
        b = _fit_logit(Xs, t).coef_.ravel()
    elif prior == "SVM":                                  # e1071::svm(kernel = 'linear', cost = 1, scale = TRUE)
        b = SVC(kernel="linear", C=1.0, tol=1e-3).fit(Xs, t).coef_.ravel()
    else:
        raise ValueError("prior must be 'LR' or 'SVM'")
    return b


def three_stage_paper(X, Y, T, prior="SVM", smoother="sigmoid", threshold=0.1, seed=1, verbose=False):
    """Three-stage selection exactly as in FS_CI_Final.R (Gaussian outcome models)."""
    T = _check_treatment(T); Y = np.asarray(Y, float)
    Xs = _std_r(X)
    beta = _exposure_coef_r(Xs, T, prior)
    pen = np.abs(beta); pen[~np.isfinite(pen)] = 1e10
    w = np.tanh(pen) ** 0.5 if smoother == "tanh" else (1 / (1 + np.exp(-pen))) ** 1.0
    s2 = cv_glmnet_r(Xs, Y, w, 0.5, 5, seed, relax=True)
    theta = s2["coef"]
    with np.errstate(divide="ignore"):
        psi = (1.0 / np.abs(theta)) ** 1.0               # zero coefficient -> Inf -> excluded
    s3 = aenet_r(Xs, Y, psi)                        # msaenet default seed = 1001
    Theta = s3["beta"]
    return dict(selected=np.flatnonzero(np.abs(Theta) >= threshold), Theta=Theta, theta_stage2=theta,
                beta_exposure=beta, w=w, psi=psi, gamma_1se=s2["gamma_1se"],
                alpha_stage3=(s3["alpha_enet"], s3["alpha_aenet"]), outcome_type="continuous (paper mode)")


def oaenet_paper(X, Y, T, ps_type="link", lam_exps=(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
                 conv=2.0, verbose=False):
    """ONet benchmark exactly as in FS_CI_Final.R. ps_type='link' reproduces predict() without
    type='response' (the R code); 'response' uses probabilities."""
    T = _check_treatment(T); Y = np.asarray(Y, float)
    Xs = _std_r(X); n, p = Xs.shape
    betaXY = np.linalg.lstsq(np.column_stack([np.ones(n), T, Xs]), Y, rcond=None)[0][2:]
    best_val, best_beta = np.nan, None
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for il in lam_exps:
            ig = 2 * (conv - il + 1)
            pre = np.log(np.abs(betaXY) ** (-ig))
            pre = np.where(pre < 0, 0.0, pre)
            nrm = np.sqrt(np.sum(pre ** 2))
            pre = 1 / (1 + np.exp(-(pre / nrm * p)))
            pre = np.where(pre > 0.9, 1e3, pre)
            pen = np.where(np.isfinite(pre), pre, 1e10)
            try:
                path = glmnet_path(Xs, T.astype(float), pen, 0.5, "binomial")
            except Exception:
                continue
            a0, b = _interp(path, n ** il)
            eta = a0 + Xs @ b
            fp = eta if ps_type == "link" else expit(eta)
            wt = np.where(T == 1, 1 / fp, 1 / (1 - fp))
            m1 = (Xs[T == 1] * wt[T == 1, None]).sum(0) / wt[T == 1].sum()
            m0 = (Xs[T == 0] * wt[T == 0, None]).sum(0) / wt[T == 0].sum()
            val = np.sum(np.abs(m1 - m0) * np.abs(betaXY))
            if np.isfinite(val) and (best_beta is None or val < best_val):   # which.min ignores NaN/NA
                best_val, best_beta, best_il = val, b, il
    if best_beta is None:
        return dict(selected=np.array([], int), beta=None, wAMD=np.nan)
    with np.errstate(invalid="ignore"):
        nb = best_beta / np.sqrt(np.sum(best_beta ** 2))
    sel = np.flatnonzero(np.abs(np.nan_to_num(nb, nan=0.0)) >= 0.01)
    return dict(selected=sel, beta=best_beta, wAMD=best_val, lam_exp=best_il, theta_outcome=betaXY)


def att_matching_paper(X, Y, T, caliper=0.25, r_single_variable_failure=True):
    """MatchIt nearest 1:1 (caliper .25 SD of the PS, discard = 'both', largest PS first) and the
    coefficient of A in lm(Y ~ A + X) on the matched data, for any outcome type.
    r_single_variable_failure=True reproduces FS_CI_Final.R, where xx[, idx] with one index drops
    the column name and the formula fails (ATT = NA)."""
    X = np.asarray(X, float)
    if X.ndim == 1:
        X = X[:, None]
    if X.shape[1] == 0 or (r_single_variable_failure and X.shape[1] == 1):
        return float("nan")
    return _att_matching_general(X, np.asarray(Y, float), T, outcome_type="continuous", caliper=caliper)


# ---------------------------------------------------------------------------
# Public functions: mode = "general" (default, earlier behaviour) or "paper"
# ---------------------------------------------------------------------------
def three_stage(X, Y, T, prior="SVM", smoother="sigmoid", gamma1=None, gamma2=None,
                outcome_type="auto", seed=1, verbose=True, mode="general", threshold=0.1):
    """mode='general': glmnet-style implementation with outcome-type detection (earlier behaviour).
    mode='paper'  : exact replication of FS_CI_Final.R (relaxed Stage 2, msaenet::aenet Stage 3,
                    |Theta| >= threshold, Gaussian outcome models)."""
    if mode == "paper":
        return three_stage_paper(X, Y, T, prior, smoother, threshold, seed, verbose)
    return _three_stage_general(X, Y, T, prior, smoother, gamma1, gamma2, outcome_type, seed, verbose)


def oaenet(X, Y, T, outcome_type="auto", lam_exps=(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
           conv=2.0, ridge_ratio=1.0, verbose=True, mode="general", ps_type="link"):
    """mode='general': OAENet with |theta|^-gamma weights and wAMD tuning (earlier behaviour).
    mode='paper'  : the ONet benchmark of FS_CI_Final.R."""
    if mode == "paper":
        return oaenet_paper(X, Y, T, ps_type, lam_exps, conv, verbose)
    return _oaenet_general(X, Y, T, outcome_type, lam_exps, conv, ridge_ratio, verbose)


def att_matching(X, Y, T, outcome_type="auto", caliper=0.25, verbose=False, mode="general"):
    """mode='general': continuous -> lm coefficient, binary -> matched risk difference.
    mode='paper'  : lm coefficient for every outcome type, R single-variable failure reproduced."""
    if mode == "paper":
        return att_matching_paper(X, Y, T, caliper)
    return _att_matching_general(X, Y, T, outcome_type, caliper, verbose)


__all__ += ["three_stage_paper", "oaenet_paper", "att_matching_paper", "glmnet_path", "cv_glmnet_r", "aenet_r"]
