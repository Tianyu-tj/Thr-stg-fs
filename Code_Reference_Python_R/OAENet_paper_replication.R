################################################################################
# OAENet paper replication -- corrected version
#
# Islam, M. S., Shikalgar, S., & Noor-E-Alam, M. (2025)
# A two-stage feature selection approach for robust evaluation of treatment
# effects in high-dimensional observational data.
#
# This script is ONLY for replicating that OAENet paper.
#
# Key correction relative to the previous script:
#   ATT is estimated from the matched sample using the post-matching outcome
#   regression used in the code supplied by the user:
#
#       Y ~ A + selected covariates
#
#   rather than using a raw matched treated-minus-control mean difference.
#
# Paper-supported components:
#   * Stage 1 outcome model: OLS for continuous Y
#   * adaptive penalty w_j = |beta_j^(OLS)|^(-c)
#   * c = 3 in the simulation study
#   * Stage 2 binary-treatment elastic-net model
#   * 5-fold cross-validation, one-standard-error rule
#   * normalization before variable selection
#   * 1:1 nearest-neighbor propensity-score matching
#
# Implementation details not numerically specified in the article text:
#   * glmnet alpha.  We retain alpha = 0.5 from the supplied prior code.
#   * caliper.       We retain caliper = 0.25 from the supplied prior code.
#
# Install packages once if needed:
# install.packages(c("glmnet", "MatchIt", "MASS"))
################################################################################


required_packages <- c("glmnet", "MatchIt", "MASS")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(glmnet)
library(MatchIt)
library(MASS)


################################################################################
# Utility functions
################################################################################

expit <- function(x) {
  1 / (1 + exp(-x))
}


make_sigma <- function(p = 100, variance = 1, rho = 0) {
  S <- matrix(rho * variance, nrow = p, ncol = p)
  diag(S) <- variance
  S
}


################################################################################
# OAENet variable selection
################################################################################

oaenet_select <- function(
    X,
    A,
    Y,
    c_gamma = 3,
    alpha_enet = 0.5,
    nfolds = 5,
    seed = NULL
) {

  X <- as.matrix(X)
  A <- as.numeric(A)
  Y <- as.numeric(Y)

  if (!all(A %in% c(0, 1))) {
    stop("Treatment A must be coded 0/1.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ---------------------------------------------------------------------------
  # Normalize covariates before variable selection.
  # ---------------------------------------------------------------------------
  X_scaled <- scale(X)

  if (any(!is.finite(X_scaled))) {
    stop("Non-finite values after scaling. Check for zero-variance columns.")
  }

  # ---------------------------------------------------------------------------
  # Stage 1: outcome model
  #
  # The paper writes:
  #       beta^(OLS) = argmin ||Y - X beta||^2
  #
  # Therefore treatment A is not added to this first-stage regression here.
  # ---------------------------------------------------------------------------
  fit_y <- lm.fit(
    x = cbind("(Intercept)" = 1, X_scaled),
    y = Y
  )

  beta_y <- fit_y$coefficients[-1]
  names(beta_y) <- colnames(X_scaled)

  # ---------------------------------------------------------------------------
  # Adaptive penalty:
  #       w_j = |beta_j^(OLS)|^(-c)
  # with c = 3 in the paper simulations.
  #
  # epsilon is only a numerical safeguard when OLS returns exactly zero.
  # ---------------------------------------------------------------------------
  epsilon <- 1e-8

  penalty_factor <- pmax(abs(beta_y), epsilon)^(-c_gamma)
  penalty_factor[!is.finite(penalty_factor)] <- 1e12
  names(penalty_factor) <- colnames(X_scaled)

  # Rescaling all penalty factors by a common positive constant does not change
  # their relative adaptive weighting and improves numerical conditioning.
  penalty_factor <- penalty_factor / mean(penalty_factor)

  # ---------------------------------------------------------------------------
  # Stage 2: treatment model
  # 5-fold CV and one-standard-error rule.
  # ---------------------------------------------------------------------------
  fit_a <- cv.glmnet(
    x = X_scaled,
    y = A,
    family = "binomial",
    alpha = alpha_enet,
    penalty.factor = penalty_factor,
    nfolds = nfolds,
    standardize = FALSE,
    type.measure = "deviance"
  )

  coef_1se <- as.matrix(
    coef(fit_a, s = "lambda.1se")
  )

  beta_a <- coef_1se[-1, 1]
  names(beta_a) <- rownames(coef_1se)[-1]

  selected_index <- which(abs(beta_a) > 0)
  selected_vars <- names(beta_a)[selected_index]

  list(
    selected_vars = selected_vars,
    selected_index = selected_index,
    beta_outcome = beta_y,
    penalty_factor = penalty_factor,
    beta_treatment = beta_a,
    lambda_1se = fit_a$lambda.1se,
    lambda_min = fit_a$lambda.min,
    cv_fit = fit_a
  )
}


################################################################################
# Matching + ATT estimation
################################################################################

estimate_att_oaenet <- function(
    Data,
    selected_vars,
    treatment = "A",
    outcome = "Y",
    ratio = 1,
    caliper = 0.25,
    replace = FALSE
) {

  if (length(selected_vars) == 0) {
    return(
      list(
        att = NA_real_,
        matchit = NULL,
        matched_data = NULL,
        n_treated_matched = 0L,
        n_control_matched = 0L
      )
    )
  }

  treat_formula <- reformulate(
    termlabels = selected_vars,
    response = treatment
  )

  mm <- try(
    suppressWarnings(
      matchit(
        treat.form,
        data = Data,
        method = "nearest",
        distance = "glm",
        ratio = 1,
        caliper = .25,
        link = "logit",
        estimand = "ATT",
        replace = FALSE,
        discard = "both"
      )
    ),
    silent = TRUE
  )

  if (inherits(mm, "try-error")) {
    return(
      list(
        att = NA_real_,
        matchit = NULL,
        matched_data = NULL,
        n_treated_matched = 0L,
        n_control_matched = 0L
      )
    )
  }

  md <- MatchIt::match.data(mm)

  # ---------------------------------------------------------------------------
  # IMPORTANT:
  # Preserve the ATT estimation step from the supplied earlier R code:
  #
  #     lm(Y ~ A + selected covariates, data = match.data(mm))
  #
  # The coefficient on A is stored as ATT.
  #
  # This is the key correction from the previous file, where I incorrectly
  # replaced this with the raw difference in matched outcome means.
  # ---------------------------------------------------------------------------
  effect_formula <- reformulate(
    termlabels = c(treatment, selected_vars),
    response = outcome
  )

  fit_effect <- try(
    lm(
      formula = effect_formula,
      data = md
    ),
    silent = TRUE
  )

  if (inherits(fit_effect, "try-error")) {
    att <- NA_real_
  } else {
    att <- unname(coef(fit_effect)[treatment])
  }

  list(
    att = att,
    matchit = mm,
    matched_data = md,
    n_treated_matched = sum(md[[treatment]] == 1),
    n_control_matched = sum(md[[treatment]] == 0),
    effect_model = fit_effect
  )
}


################################################################################
# Complete OAENet fit for one generated dataset
################################################################################

run_oaenet <- function(
    Data,
    covariates,
    treatment = "A",
    outcome = "Y",
    c_gamma = 3,
    alpha_enet = 0.5,
    nfolds = 5,
    ratio = 1,
    caliper = 0.25,
    seed = NULL
) {

  X <- as.matrix(Data[, covariates, drop = FALSE])
  A <- Data[[treatment]]
  Y <- Data[[outcome]]

  t0 <- proc.time()[["elapsed"]]

  selection <- try(
    oaenet_select(
      X = X,
      A = A,
      Y = Y,
      c_gamma = c_gamma,
      alpha_enet = alpha_enet,
      nfolds = nfolds,
      seed = seed
    ),
    silent = TRUE
  )

  if (inherits(selection, "try-error")) {
    return(
      list(
        att = NA_real_,
        selected_vars = character(0),
        selected_index = integer(0),
        n_selected = NA_integer_,
        n_treated_matched = NA_integer_,
        n_control_matched = NA_integer_,
        elapsed_seconds = proc.time()[["elapsed"]] - t0,
        error = as.character(selection)
      )
    )
  }

  att_fit <- estimate_att_oaenet(
    Data = Data,
    selected_vars = selection$selected_vars,
    treatment = treatment,
    outcome = outcome,
    ratio = ratio,
    caliper = caliper,
    replace = FALSE
  )

  list(
    att = att_fit$att,
    selected_vars = selection$selected_vars,
    selected_index = selection$selected_index,
    n_selected = length(selection$selected_vars),
    beta_outcome = selection$beta_outcome,
    penalty_factor = selection$penalty_factor,
    beta_treatment = selection$beta_treatment,
    lambda_1se = selection$lambda_1se,
    lambda_min = selection$lambda_min,
    n_treated_matched = att_fit$n_treated_matched,
    n_control_matched = att_fit$n_control_matched,
    matchit = att_fit$matchit,
    matched_data = att_fit$matched_data,
    effect_model = att_fit$effect_model,
    elapsed_seconds = proc.time()[["elapsed"]] - t0
  )
}


################################################################################
# Paper simulation Scenario 1
################################################################################

simulate_scenario1 <- function(
    n = 1000,
    rho = 0,
    seed = NULL
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- 100

  Sigma <- make_sigma(
    p = p,
    variance = 1,
    rho = rho
  )

  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma,
    empirical = FALSE
  )

  colnames(X) <- paste0("X", 1:p)

  # Treatment equation in the paper:
  #
  # A ~ Bernoulli[
  #   expit(X1 + ... + X10 + X21 + ... + X30)
  # ]
  eta_A <-
    rowSums(X[, 1:10, drop = FALSE]) +
    rowSums(X[, 21:30, drop = FALSE])

  pr_A <- expit(eta_A)

  A <- rbinom(
    n = n,
    size = 1,
    prob = pr_A
  )

  # Outcome equation printed in the paper:
  #
  # Y = 0.5 A + 0.6 X1 + ... + 0.6 X20
  #
  # No additional stochastic error term is added because the article's
  # displayed simulation equation does not specify one.
  Y <-
    0.5 * A +
    0.6 * rowSums(X[, 1:20, drop = FALSE])

  Data <- data.frame(
    Y = Y,
    A = A,
    X,
    check.names = FALSE
  )

  list(
    data = Data,
    true_te = 0.5,
    confounders = paste0("X", 1:10),
    outcome_predictors = paste0("X", 11:20),
    treatment_predictors = paste0("X", 21:30),
    target_vars = paste0("X", 1:20)
  )
}


################################################################################
# Paper simulation Scenario 2
################################################################################

simulate_scenario2 <- function(
    n = 1000,
    rho = 0,
    seed = NULL
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- 100

  Sigma <- make_sigma(
    p = p,
    variance = 4,
    rho = rho
  )

  X <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma,
    empirical = FALSE
  )

  colnames(X) <- paste0("X", 1:p)

  # Treatment equation printed in the paper.
  eta_A <-
    0.5 * X[, 1] -
    1.0 * X[, 2] +
    0.3 * X[, 5] -
    0.3 * X[, 6] +
    0.3 * X[, 7] -
    0.3 * X[, 8]

  pr_A <- expit(eta_A)

  A <- rbinom(
    n = n,
    size = 1,
    prob = pr_A
  )

  # Outcome equation printed in the paper.
  Y <-
    1.0 * A +
    2.0 * X[, 1] +
    2.0 * X[, 2] +
    5.0 * X[, 3] +
    5.0 * X[, 4]

  Data <- data.frame(
    Y = Y,
    A = A,
    X,
    check.names = FALSE
  )

  list(
    data = Data,
    true_te = 1.0,
    confounders = paste0("X", 1:2),
    outcome_predictors = paste0("X", 3:4),
    target_vars = paste0("X", 1:4)
  )
}


################################################################################
# Scenario settings
################################################################################

paper_scenarios <- data.frame(
  scenario = c(
    "1A", "1B", "1C", "1D", "1E", "1F",
    "2A", "2B", "2C", "2D", "2E", "2F"
  ),
  dgp = c(rep(1, 6), rep(2, 6)),
  n = c(
    1000, 1000, 1000, 500, 500, 500,
    1000, 1000, 1000, 500, 500, 500
  ),
  rho = rep(
    c(0, 0.25, 0.50, 0, 0.25, 0.50),
    2
  )
)


################################################################################
# Repeat one scenario
################################################################################

replicate_one_scenario <- function(
    scenario_id,
    n_rep = 1000,
    seed = 100,
    c_gamma = 3,
    alpha_enet = 0.5,
    nfolds = 5,
    caliper = 0.25,
    verbose = TRUE
) {

  setting <- paper_scenarios[
    paper_scenarios$scenario == scenario_id,
    ,
    drop = FALSE
  ]

  if (nrow(setting) != 1) {
    stop("Unknown scenario_id.")
  }

  covariates <- paste0("X", 1:100)

  att <- rep(NA_real_, n_rep)
  runtime <- rep(NA_real_, n_rep)
  n_selected <- rep(NA_real_, n_rep)
  n_matched_treated <- rep(NA_real_, n_rep)
  n_matched_control <- rep(NA_real_, n_rep)

  selection_count <- setNames(
    rep(0, length(covariates)),
    covariates
  )

  for (r in seq_len(n_rep)) {

    this_seed <- seed + r - 1L

    sim <- if (setting$dgp == 1) {

      simulate_scenario1(
        n = setting$n,
        rho = setting$rho,
        seed = this_seed
      )

    } else {

      simulate_scenario2(
        n = setting$n,
        rho = setting$rho,
        seed = this_seed
      )
    }

    fit <- try(
      run_oaenet(
        Data = sim$data,
        covariates = covariates,
        treatment = "A",
        outcome = "Y",
        c_gamma = c_gamma,
        alpha_enet = alpha_enet,
        nfolds = nfolds,
        ratio = 1,
        caliper = caliper,
        seed = this_seed
      ),
      silent = TRUE
    )

    if (!inherits(fit, "try-error")) {

      att[r] <- fit$att
      runtime[r] <- fit$elapsed_seconds
      n_selected[r] <- fit$n_selected
      n_matched_treated[r] <- fit$n_treated_matched
      n_matched_control[r] <- fit$n_control_matched

      if (length(fit$selected_vars) > 0) {
        selection_count[fit$selected_vars] <-
          selection_count[fit$selected_vars] + 1
      }
    }

    if (
      verbose &&
      (r %% max(1, floor(n_rep / 20)) == 0 || r == n_rep)
    ) {
      message(
        sprintf(
          "Scenario %s: %d / %d completed",
          scenario_id,
          r,
          n_rep
        )
      )
    }
  }

  true_te <- if (setting$dgp == 1) 0.5 else 1.0
  valid_att <- att[is.finite(att)]

  summary_row <- data.frame(
    scenario = scenario_id,
    n_rep = n_rep,
    successful_att = length(valid_att),
    true_te = true_te,
    mean_att = if (length(valid_att) > 0) {
      mean(valid_att)
    } else {
      NA_real_
    },
    bias = if (length(valid_att) > 0) {
      mean(valid_att) - true_te
    } else {
      NA_real_
    },
    variance = if (length(valid_att) > 1) {
      var(valid_att)
    } else {
      NA_real_
    },
    mean_runtime_seconds = mean(runtime, na.rm = TRUE),
    mean_n_selected = mean(n_selected, na.rm = TRUE),
    mean_matched_treated = mean(n_matched_treated, na.rm = TRUE),
    mean_matched_control = mean(n_matched_control, na.rm = TRUE)
  )

  list(
    summary = summary_row,
    att = att,
    selection_rate = selection_count / n_rep,
    n_selected = n_selected,
    runtime = runtime
  )
}


################################################################################
# Repeat all 12 paper scenarios
################################################################################

replicate_all_scenarios <- function(
    n_rep = 1000,
    seed = 100,
    c_gamma = 3,
    alpha_enet = 0.5,
    nfolds = 5,
    caliper = 0.25,
    verbose = TRUE
) {

  out <- vector(
    "list",
    length = nrow(paper_scenarios)
  )

  names(out) <- paper_scenarios$scenario

  for (i in seq_len(nrow(paper_scenarios))) {

    sid <- paper_scenarios$scenario[i]

    message("Starting Scenario ", sid)

    out[[sid]] <- replicate_one_scenario(
      scenario_id = sid,
      n_rep = n_rep,
      seed = seed + (i - 1L) * 100000L,
      c_gamma = c_gamma,
      alpha_enet = alpha_enet,
      nfolds = nfolds,
      caliper = caliper,
      verbose = verbose
    )
  }

  summary_table <- do.call(
    rbind,
    lapply(out, function(z) z$summary)
  )

  rownames(summary_table) <- NULL

  list(
    summary = summary_table,
    results = out
  )
}


################################################################################
# Diagnostic comparison to the published Scenario 1A OAENet result
################################################################################

check_scenario_1A <- function(result) {

  if (is.list(result) && !is.null(result$summary)) {
    s <- result$summary
  } else {
    s <- result
  }

  paper_mean_att <- 0.5031
  paper_bias <- 0.0031
  paper_variance <- 0.0090

  comparison <- data.frame(
    metric = c("mean_att", "bias", "variance"),
    replication = c(
      s$mean_att[1],
      s$bias[1],
      s$variance[1]
    ),
    paper = c(
      paper_mean_att,
      paper_bias,
      paper_variance
    )
  )

  comparison$difference <-
    comparison$replication - comparison$paper

  comparison
}


################################################################################
# QUICK TEST
#
# Run this first.
################################################################################

# test_1A <- replicate_one_scenario(
#   scenario_id = "1A",
#   n_rep = 20,
#   seed = 100
# )
#
# print(test_1A$summary)
# print(check_scenario_1A(test_1A))
#
# To inspect selected variables:
# sort(test_1A$selection_rate, decreasing = TRUE)


################################################################################
# FULL SCENARIO 1A REPLICATION
################################################################################

result_1A <- replicate_one_scenario(
  scenario_id = "1A",
  n_rep = 1000,
  seed = 100
)

print(result_1A$summary)
print(check_scenario_1A(result_1A))


################################################################################
# FULL 12-SCENARIO REPLICATION
################################################################################

# full_result <- replicate_all_scenarios(
#   n_rep = 1000,
#   seed = 100
# )
# 
# print(full_result$summary)
# 
# write.csv(
#   full_result$summary,
#   "OAENet_replication_summary.csv",
#   row.names = FALSE
)
################################################################################

# test_1A <- replicate_one_scenario(
#   scenario_id = "1A",
#   n_rep = 20,
#   seed = 100
# )
# 
# test_1A$summary