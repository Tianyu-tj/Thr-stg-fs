# simulate_demo.R - synthetic-data demonstration of the three-stage framework and OAENet.
#
# Data-generating process (Section 6.1 of the three-stage paper; Shortreed & Ertefaie 2017):
#   X ~ N(0, Sigma), Sigma_jk = rho (j != k), p = 20
#   T ~ Bernoulli(expit(X beta))
#   continuous outcome: Y = 0.5 T + X theta + eps,          eps ~ N(0, 1)
#   binary outcome:     Y ~ Bernoulli(expit(-1 + 0.5 T + X theta))
# Confounders X1-X2, pure outcome predictors X3-X4, pure treatment predictors X5-X6, noise X7-X20.
# Desired selection: {X1, X2, X3, X4}.
#
# Usage:  Rscript simulate_demo.R            (edit REPS / N below)

source("causal_fs.R")
suppressPackageStartupMessages(library(MASS))

REPS <- 10; N <- 500; RHOS <- c(0.25, 0.75); SCENARIOS <- 1:4; OUTCOMES <- c("continuous", "binary")

scenario_coefs <- function(s) {
  z <- rep(0, 14)
  switch(s,
    list(theta = c(.6, .6, .6, .6, rep(0, 16)), beta = c(1, 1, 0, 0, 1, 1, z)),
    list(theta = c(.6, .6, .6, .6, rep(0, 16)), beta = c(.4, .4, 0, 0, 1, 1, z)),
    list(theta = c(.2, .2, .6, .6, rep(0, 16)), beta = c(1, 1, 0, 0, 1, 1, z)),
    list(theta = c(.6, .6, .6, .6, rep(0, 16)), beta = c(1, 1, 0, 0, 1.8, 1.8, z)))
}

simulate <- function(n, rho, scenario, seed, outcome = "continuous", p = 20) {
  set.seed(seed)
  cf <- scenario_coefs(scenario)
  S <- matrix(rho, p, p); diag(S) <- 1
  X <- mvrnorm(n, rep(0, p), S)
  T <- rbinom(n, 1, plogis(X %*% cf$beta))
  if (outcome == "continuous") {
    Y <- 0.5 * T + as.numeric(X %*% cf$theta) + rnorm(n); true_att <- 0.5
  } else {
    eta0 <- -1 + as.numeric(X %*% cf$theta)
    Y <- rbinom(n, 1, plogis(eta0 + 0.5 * T))
    true_att <- mean(plogis(eta0[T == 1] + 0.5) - plogis(eta0[T == 1]))
  }
  list(X = X, Y = Y, T = T, true_att = true_att)
}

methods <- list(
  "Enh-ESVMS" = function(d, s) three_stage(d$X, d$Y, d$T, "SVM", "sigmoid", seed = s, verbose = FALSE)$selected,
  "Enh-ELRT"  = function(d, s) three_stage(d$X, d$Y, d$T, "LR", "tanh", seed = s, verbose = FALSE)$selected,
  "OAENet"    = function(d, s) oaenet(d$X, d$Y, d$T, verbose = FALSE)$selected
)

rows <- list()
for (outcome in OUTCOMES) for (sc in SCENARIOS) for (rho in RHOS) for (seed in seq_len(REPS)) {
  d <- simulate(N, rho, sc, seed, outcome)
  target <- att_matching(d$X[, 1:4], d$Y, d$T)
  for (m in names(methods)) {
    sel <- methods[[m]](d, seed)
    att <- if (length(sel)) att_matching(d$X[, sel, drop = FALSE], d$Y, d$T) else NA
    rows[[length(rows) + 1]] <- data.frame(
      outcome = outcome, scenario = sc, rho = rho, seed = seed, method = m,
      X1_X4_all = as.integer(all(1:4 %in% sel)), X5_X6_any = as.integer(any(5:6 %in% sel)),
      n_noise = sum(sel > 6), att = att, target_att = target, true_att = d$true_att,
      abs_bias_vs_target = abs(att - target), abs_bias_vs_truth = abs(att - d$true_att))
  }
  cat(outcome, "scenario", sc, "rho", rho, "seed", seed, "\n")
}
res <- do.call(rbind, rows)
write.csv(res, "simulation_results_R.csv", row.names = FALSE)
summ <- aggregate(cbind(X1_X4_all, X5_X6_any, n_noise, abs_bias_vs_target, abs_bias_vs_truth) ~ outcome + method + scenario,
                  data = res, FUN = mean)
print(summ, digits = 3)
write.csv(summ, "simulation_summary_R.csv", row.names = FALSE)
