# causal_fs.R
# ============
# Covariate selection for causal inference
#   * Enhanced three-stage framework (Yang & Noor-E-Alam)
#   * Outcome-adaptive elastic net, OAENet (Islam, Shikalgar & Noor-E-Alam, 2024), tuned by wAMD
#     (Shortreed & Ertefaie, 2017)
#
# Outcome type is detected automatically:
#   binary outcome     -> logistic (binomial) outcome models
#   continuous outcome -> linear (Gaussian) outcome models
# or set explicitly with outcome_type = "binary" / "continuous".
#
# Packages: glmnet, e1071, MatchIt (MatchIt only for att_matching()).

suppressPackageStartupMessages({
  library(glmnet)
  library(e1071)
})

CFG <- list(nfolds = 5, alpha_stage2 = 0.5, alpha_stage3 = 0.5,
            gamma1_sigmoid = 1.0, gamma1_tanh = 0.5, gamma2 = 1.0, rule = "lambda.1se")

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
detect_outcome_type <- function(y, outcome_type = "auto", verbose = TRUE) {
  if (anyNA(y)) stop("Outcome contains missing values; handle them before selection.")
  y <- as.numeric(y)
  vals <- sort(unique(y))
  if (outcome_type == "auto") {
    if (length(vals) < 2) stop("Outcome is constant.")
    outcome_type <- if (length(vals) == 2) "binary" else "continuous"
    if (outcome_type == "continuous" && length(vals) <= 5 && verbose)
      warning(sprintf("Outcome has only %d distinct values; treated as continuous. Set outcome_type explicitly if not intended.",
                      length(vals)))
  }
  if (outcome_type == "binary") {
    if (length(vals) != 2) stop(sprintf("outcome_type = 'binary' but the outcome has %d distinct values.", length(vals)))
    if (!identical(vals, c(0, 1))) {
      if (verbose) warning(sprintf("Binary outcome recoded: %s -> 0, %s -> 1.", vals[1], vals[2]))
      y <- as.numeric(y == vals[2])
    }
  } else if (outcome_type != "continuous") {
    stop("outcome_type must be 'auto', 'binary' or 'continuous'.")
  }
  if (verbose) message(sprintf("[causal_fs] outcome type: %s -> %s outcome models", outcome_type,
                               if (outcome_type == "binary") "logistic" else "linear"))
  list(type = outcome_type, y = y)
}

check_treatment <- function(t) {
  t <- as.numeric(t)
  if (!identical(sort(unique(t)), c(0, 1))) stop("Treatment must be coded 0/1.")
  t
}

standardize <- function(X) {
  X <- as.matrix(X); sd <- apply(X, 2, sd); sd[sd < 1e-12] <- 1
  scale(X, center = TRUE, scale = sd)
}

# cv.glmnet with penalty factors; Inf (or zero-coefficient) entries are excluded
cv_enet <- function(X, y, pf, alpha, family, seed = 1) {
  excl <- which(!is.finite(pf) | pf >= 1e10)
  if (length(excl) == ncol(X)) return(rep(0, ncol(X)))
  pf2 <- pf; pf2[excl] <- 1
  set.seed(seed)
  fit <- cv.glmnet(X, y, family = family, alpha = alpha, penalty.factor = pf2,
                   exclude = if (length(excl)) excl else NULL, nfolds = CFG$nfolds,
                   type.measure = if (family == "binomial") "deviance" else "mse", standardize = TRUE)
  as.numeric(coef(fit, s = CFG$rule))[-1]
}

# ---------------------------------------------------------------------------
# Enhanced three-stage framework (Algorithm 1)
# ---------------------------------------------------------------------------
exposure_coef <- function(Xs, t, prior) {
  if (prior == "LR") {
    b <- coef(glm(t ~ Xs, family = binomial()))[-1]
  } else if (prior == "SVM") {   # linear SVM, cost = 1 (libsvm hinge loss); X already standardised
    fit <- svm(x = Xs, y = factor(t), kernel = "linear", cost = 1, scale = FALSE, type = "C-classification")
    b <- as.numeric(t(fit$coefs) %*% fit$SV)
  } else stop("prior must be 'LR' or 'SVM'")
  b[!is.finite(b)] <- 0
  b
}

smooth_penalty <- function(beta, smoother, gamma1) {
  a <- abs(beta)
  if (smoother == "sigmoid") return((1 / (1 + exp(-a)))^gamma1)
  if (smoother == "tanh") return(pmax(tanh(a), 1e-8)^gamma1)
  stop("smoother must be 'sigmoid' or 'tanh'")
}

three_stage <- function(X, Y, T, prior = "SVM", smoother = "sigmoid", gamma1 = NULL, gamma2 = NULL,
                        outcome_type = "auto", seed = 1, verbose = TRUE) {
  ot <- detect_outcome_type(Y, outcome_type, verbose); Y <- ot$y
  T <- check_treatment(T)
  family <- if (ot$type == "binary") "binomial" else "gaussian"
  if (is.null(gamma1)) gamma1 <- if (smoother == "sigmoid") CFG$gamma1_sigmoid else CFG$gamma1_tanh
  if (is.null(gamma2)) gamma2 <- CFG$gamma2
  Xs <- standardize(X)
  beta  <- exposure_coef(Xs, T, prior)                                   # Stage 1
  w     <- smooth_penalty(beta, smoother, gamma1)                        # penalty smoothing
  theta <- cv_enet(Xs, Y, w, CFG$alpha_stage2, family, seed)             # Stage 2, Eq. (12)
  psi   <- ifelse(abs(theta) > 0, (1 / abs(theta))^gamma2, Inf)
  Theta <- cv_enet(Xs, Y, psi, CFG$alpha_stage3, family, seed + 1)       # Stage 3, Eq. (13)
  list(selected = which(abs(Theta) > 1e-10), Theta = Theta, theta_stage2 = theta,
       beta_exposure = beta, w = w, psi = psi, outcome_type = ot$type)
}

# ---------------------------------------------------------------------------
# OAENet with wAMD tuning
# ---------------------------------------------------------------------------
# min mean logistic loss + sum_j l1w_j |b_j| + l2/2 ||b||^2   (FISTA; intercept unpenalised)
logistic_enet_fista <- function(X, t, l1w, l2, max_iter = 5000, tol = 1e-8) {
  n <- nrow(X); p <- ncol(X)
  b0 <- log(mean(t) / (1 - mean(t))); b <- rep(0, p)
  L <- max(svd(X, nu = 0, nv = 0)$d)^2 / (4 * n) + l2 + 1e-8
  z0 <- b0; z <- b; tk <- 1
  for (it in seq_len(max_iter)) {
    r <- plogis(z0 + X %*% z)[, 1] - t
    g0 <- mean(r); g <- as.numeric(crossprod(X, r)) / n + l2 * z
    nb0 <- z0 - g0 / L
    u <- z - g / L
    nb <- sign(u) * pmax(abs(u) - l1w / L, 0)
    d <- max(abs(nb0 - b0), if (p) max(abs(nb - b)) else 0)
    tn <- (1 + sqrt(1 + 4 * tk^2)) / 2
    z0 <- nb0 + (tk - 1) / tn * (nb0 - b0); z <- nb + (tk - 1) / tn * (nb - b)
    b0 <- nb0; b <- nb; tk <- tn
    if (d < tol) break
  }
  list(b0 = b0, b = b)
}

wamd <- function(Xs, t, ps, theta) {
  w <- ifelse(t == 1, 1 / ps, 1 / (1 - ps))
  m1 <- colSums(w[t == 1] * Xs[t == 1, , drop = FALSE]) / sum(w[t == 1])
  m0 <- colSums(w[t == 0] * Xs[t == 0, , drop = FALSE]) / sum(w[t == 0])
  sum(abs(theta) * abs(m1 - m0))
}

oaenet <- function(X, Y, T, outcome_type = "auto",
                   lam_exps = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
                   conv = 2, ridge_ratio = 1, verbose = TRUE) {
  ot <- detect_outcome_type(Y, outcome_type, verbose); Y <- ot$y
  T <- check_treatment(T)
  Xs <- standardize(X); n <- nrow(Xs); p <- ncol(Xs)
  theta <- if (ot$type == "binary") coef(glm(Y ~ T + Xs, family = binomial()))[-(1:2)]
           else coef(lm(Y ~ T + Xs))[-(1:2)]
  theta[!is.finite(theta)] <- 0
  best <- NULL
  for (le in lam_exps) {
    lam_sum <- n^le
    gam <- 2 * (conv - le + 1)
    logw <- log(lam_sum / n) - gam * log(pmax(abs(theta), 1e-300))
    keep <- logw < log(0.6)   # exact screen: |gradient| <= 0.5 for standardised X
    b <- rep(0, p)
    if (any(keep)) b[keep] <- logistic_enet_fista(Xs[, keep, drop = FALSE], T, exp(logw[keep]),
                                                  ridge_ratio * lam_sum / n)$b
    eta <- as.numeric(Xs %*% b); b0 <- log(mean(T) / (1 - mean(T)))
    for (k in 1:50) { pr <- plogis(b0 + eta); b0 <- b0 - mean(pr - T) / (mean(pr * (1 - pr)) + 1e-12) }
    ps <- pmin(pmax(plogis(b0 + eta), 1e-6), 1 - 1e-6)
    val <- wamd(Xs, T, ps, theta)
    if (is.null(best) || val < best$wAMD) best <- list(wAMD = val, lam_exp = le, gamma = gam, beta = b)
  }
  best$selected <- which(abs(best$beta) > 1e-10)
  best$theta_outcome <- theta
  best$outcome_type <- ot$type
  best
}

# ---------------------------------------------------------------------------
# ATT by propensity-score matching on the selected covariates
# ---------------------------------------------------------------------------
att_matching <- function(X, Y, T, outcome_type = "auto", caliper = 0.25) {
  suppressPackageStartupMessages(library(MatchIt))
  ot <- detect_outcome_type(Y, outcome_type, verbose = FALSE); Y <- ot$y
  X <- as.matrix(X)
  if (ncol(X) == 0) return(mean(Y[T == 1]) - mean(Y[T == 0]))
  colnames(X) <- paste0("X", seq_len(ncol(X)))
  d <- data.frame(Y = Y, T = T, X)
  form <- as.formula(paste("T ~", paste(colnames(X), collapse = " + ")))
  m <- matchit(form, data = d, method = "nearest", caliper = caliper, discard = "both", estimand = "ATT")
  md <- match.data(m)
  if (ot$type == "binary") {            # matched-sample risk difference
    return(weighted.mean(md$Y[md$T == 1], md$weights[md$T == 1]) - weighted.mean(md$Y[md$T == 0], md$weights[md$T == 0]))
  }
  fit <- lm(as.formula(paste("Y ~ T +", paste(colnames(X), collapse = " + "))), data = md, weights = weights)
  unname(coef(fit)["T"])
}
