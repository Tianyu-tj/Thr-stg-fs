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
# Packages: glmnet, e1071, MatchIt (MatchIt only for att_matching()); msaenet for mode = "paper".
# mode = "paper" reproduces FS_CI_Final.R exactly (see the PAPER MODE section at the end).

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

three_stage_general <- function(X, Y, T, prior = "SVM", smoother = "sigmoid", gamma1 = NULL, gamma2 = NULL,
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

oaenet_general <- function(X, Y, T, outcome_type = "auto",
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
att_matching_general <- function(X, Y, T, outcome_type = "auto", caliper = 0.25) {
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


# ===========================================================================
# PAPER MODE: code of FS_CI_Final.R, wrapped as functions
# ===========================================================================
# Stage 2: cv.glmnet(relax = TRUE), coef at lambda.1se; Stage 3: msaenet::aenet with an alpha grid;
# selection |Theta| >= 0.1; outcome models always Gaussian; ONet = the benchmark of FS_CI_Final.R;
# ATT = MatchIt nearest 1:1 (caliper .25, discard = "both") + lm(Y ~ A + X).
# X is standardised with scale() first, as the simulation code does before selection.

paper_data <- function(X, Y, T) {
  xx <- as.data.frame(scale(as.matrix(X)))
  if (is.null(colnames(X))) colnames(xx) <- paste0("V", seq_len(ncol(xx))) else colnames(xx) <- make.names(colnames(X))
  xa <- cbind(xx, A = as.numeric(T))
  list(xx = xx, xa = xa, yy = as.numeric(Y), aa = as.numeric(T),
       Data = cbind(xa, Y = as.numeric(Y)), p = ncol(xx), var.list = colnames(xx))
}

three_stage_paper <- function(X, Y, T, prior = "SVM", smoother = "sigmoid", threshold = 0.1) {
  suppressPackageStartupMessages(library(msaenet))
  d <- paper_data(X, Y, T); xx <- d$xx; xa <- d$xa; yy <- d$yy; p <- d$p
  if (prior == "LR") {
    fit1 <- glm(A ~ ., family = binomial(link = 'logit'), data = xa)
  } else {
    fit1 <- svm(A ~ ., data = xa, type = 'C-classification', kernel = 'linear', scale = TRUE)
  }
  weight_logit <- coef(fit1)[-c(1)]
  penalty <- abs(weight_logit)
  penalty[is.na(penalty)] <- 10000000000
  penalty <- if (smoother == "tanh") tanh(penalty)^0.5 else (1/(1+exp(-penalty)))^1
  fit2 <- suppressWarnings(cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                   nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                   relax = TRUE, parallel = TRUE, standardize = TRUE,
                   type.measure = 'mse', penalty.factor = penalty))
  beta_el <- as.matrix(coef(fit2, s = fit2$lambda.1se))
  beta_el_allvar <- as.matrix(beta_el[2:(p+1)])
  penalty_adel <- (1/abs(beta_el_allvar))^1
  penalty_adel[is.na(penalty_adel)] <- 10000000000
  fit3 <- suppressWarnings(aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                nfolds = 5L, rule = "lambda.1se",
                penalty.factor = penalty_adel, parallel = TRUE))
  beta_aden_raw <- as.matrix(coef(fit3, s = fit3$lambda.1se))
  Theta <- as.numeric(beta_aden_raw[1:p])
  list(selected = which(abs(Theta) >= threshold), Theta = Theta, theta_stage2 = as.numeric(beta_el_allvar),
       beta_exposure = as.numeric(weight_logit), w = as.numeric(penalty), psi = as.numeric(penalty_adel),
       outcome_type = "continuous (paper mode)")
}

oaenet_paper <- function(X, Y, T, ps_type = "link",
                         lambda_vec = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49), conv = 2) {
  d <- paper_data(X, Y, T); xx <- d$xx; aa <- d$aa; Data <- d$Data; p <- d$p; var.list <- d$var.list
  n <- nrow(xx)
  names(lambda_vec) <- as.character(lambda_vec)
  gamma_vals <- 2 * (conv - lambda_vec + 1); names(gamma_vals) <- names(lambda_vec)
  vec_norm <- function(x) x / sqrt(sum(x^2))
  vec_sigmoid <- function(pen) 1/(1+exp(-pen))
  create_weights <- function(fp, fA) { fw <- (fp)^(-1); fw[fA == 0] <- (1 - fp[fA == 0])^(-1); fw }
  wAMD_function <- function(DataM, varlist, trt.var, wgt, beta) {
    diff_vec <- rep(NA, length(beta))
    for (jj in 1:length(varlist)) {
      this.var <- paste("w", varlist[jj], sep = "")
      DataM[, this.var] <- DataM[, varlist[jj]] * DataM[, wgt]
      trt <- sum(DataM[DataM[, trt.var] == 1, this.var]) / sum(DataM[DataM[, trt.var] == 1, wgt])
      untrt <- sum(DataM[DataM[, trt.var] == 0, this.var]) / sum(DataM[DataM[, trt.var] == 0, wgt])
      diff_vec[jj] <- abs(trt - untrt)
    }
    sum(diff_vec * abs(beta))
  }
  betaXY <- coef(lm(formula(paste("Y~A+", paste(var.list, collapse = "+"))), data = Data))[var.list]
  dummy_df <- Data
  coeff_XA <- as.data.frame(matrix(NA, nrow = p, ncol = length(lambda_vec)))
  names(coeff_XA) <- names(lambda_vec); rownames(coeff_XA) <- var.list
  wAMD_vec <- rep(NA, length(lambda_vec)); names(wAMD_vec) <- names(lambda_vec)
  for (lil in names(lambda_vec)) {
    il <- lambda_vec[lil]; ig <- gamma_vals[lil]
    prearray <- log(abs(betaXY)^(-ig))
    prearray <- ifelse(prearray < 0, 0, prearray)
    prearray <- vec_sigmoid(vec_norm(prearray) * p)
    prearray <- ifelse(prearray > 0.9, 10^3, prearray)
    penalty <- prearray
    penalty[which(!is.finite(penalty))] <- 10000000000
    wAMD_vec[lil] <- Inf
    fit2 <- try(cv.glmnet(data.matrix(xx), as.matrix(aa), alpha = 0.5, nfolds = 5,
                          penalty.factor = penalty, family = "binomial"), silent = TRUE)
    if ('try-error' %in% class(fit2)) { tmp <- Inf } else {
      dummy_df[, paste("f.pA", lil, sep = "")] <- predict(fit2, data.matrix(xx), s = n^(il), type = ps_type)
      coeff_XA[var.list, lil] <- coef(fit2, s = n^(il))[2:(p+1)]
      dummy_df[, paste("w", lil, sep = "")] <- create_weights(fp = dummy_df[, paste("f.pA", lil, sep = "")], fA = aa)
      tmp <- wAMD_function(DataM = dummy_df, varlist = var.list, trt.var = "A",
                           wgt = paste("w", lil, sep = ""), beta = betaXY)
    }
    wAMD_vec[lil] <- min(wAMD_vec[lil], tmp)
  }
  tt <- which.min(wAMD_vec)
  beta_ps_allvar <- as.matrix(vec_norm(coeff_XA[, names(tt)]))
  sel <- row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.01)]
  list(selected = sel, beta = coeff_XA[, names(tt)], wAMD = wAMD_vec, lam_exp = lambda_vec[names(tt)],
       theta_outcome = betaXY)
}

att_matching_paper <- function(X, Y, T, selected = NULL, caliper = .25) {
  # X: full covariate matrix; selected: column indices (as returned by the selection functions)
  suppressPackageStartupMessages(library(MatchIt))
  d <- paper_data(X, Y, T); xx <- d$xx; Data <- d$Data
  if (is.null(selected)) selected <- seq_len(ncol(xx))
  xx_true <- xx[, selected]                        # one index -> vector without names -> NA, as in the R code
  vl <- names(xx_true)
  treat.form <- try(formula(paste("A~", paste(vl, collapse = "+"))), silent = TRUE)
  if ('try-error' %in% class(treat.form)) return(NA)
  mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                    ratio = 1, caliper = caliper, link = "logit", estimand = "ATT",
                    replace = FALSE, discard = "both"), silent = TRUE)
  if ('try-error' %in% class(mm)) return(NA)
  fit.effect <- lm(formula(paste("Y~A+", paste(vl, collapse = "+"))), data = match.data(mm))
  unname(fit.effect$coefficients[2])
}

# ---------------------------------------------------------------------------
# Public functions: mode = "general" (default, earlier behaviour) or "paper"
# ---------------------------------------------------------------------------
three_stage <- function(X, Y, T, prior = "SVM", smoother = "sigmoid", gamma1 = NULL, gamma2 = NULL,
                        outcome_type = "auto", seed = 1, verbose = TRUE, mode = "general", threshold = 0.1) {
  if (mode == "paper") return(three_stage_paper(X, Y, T, prior, smoother, threshold))
  three_stage_general(X, Y, T, prior, smoother, gamma1, gamma2, outcome_type, seed, verbose)
}

oaenet <- function(X, Y, T, outcome_type = "auto",
                   lam_exps = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49),
                   conv = 2, ridge_ratio = 1, verbose = TRUE, mode = "general", ps_type = "link") {
  if (mode == "paper") return(oaenet_paper(X, Y, T, ps_type, lam_exps, conv))
  oaenet_general(X, Y, T, outcome_type, lam_exps, conv, ridge_ratio, verbose)
}

# general mode: att_matching(X_selected, Y, T); paper mode: att_matching(X_full, Y, T, selected = idx, mode = "paper")
att_matching <- function(X, Y, T, outcome_type = "auto", caliper = 0.25, mode = "general", selected = NULL) {
  if (mode == "paper") return(att_matching_paper(X, Y, T, selected, caliper))
  if (!is.null(selected)) X <- as.matrix(X)[, selected, drop = FALSE]
  att_matching_general(X, Y, T, outcome_type, caliper)
}
