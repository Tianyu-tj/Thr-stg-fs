## simulate_demo.R - paper-mode replication of the synthetic experiments (FS_CI_Final.R, Section 6.1)
##
## X ~ N(0, Sigma), Sigma_jk = rho, p = 20, X standardised after generation
## A = 1{U < expit(X alpha)},  Y = X beta + N(0, 1) + 0.5 A
## coefficients multiplied by COEF_MULT (2 in FS_CI_Final.R; 1 = values printed in Section 6.1)
## methods: Enh-ELRT, Enh-ESVMS (three_stage, mode = "paper") and ONet (oaenet, mode = "paper")
## metric : percentage selection bias |ATT - ATT_target| / 0.5 * 100 and per-variable selection rates
##
## EXPORT_DIR (optional) writes every data set and R's selections, so that
##   python simulate_demo.py --from-r <EXPORT_DIR>
## can compare the Python implementation with R on identical data.

source("causal_fs.R")
suppressPackageStartupMessages({ library(MASS); library(msaenet); library(MatchIt) })

COEF_MULT  <- 2
SCENARIOS  <- 1:4
RHO_LIST   <- c(0, 0.25, 0.5, 0.75)
ITR_MAX    <- 30
N          <- 200
P          <- 20
OUT_FILE   <- "simulation_results_paper_R.csv"
EXPORT_DIR <- NULL            # e.g. "r_export"

scen <- list(
  `1` = list(beta = c(.6, .6, .6, .6), alpha = c(1, 1, 0, 0, 1, 1)),
  `2` = list(beta = c(.6, .6, .6, .6), alpha = c(.4, .4, 0, 0, 1, 1)),
  `3` = list(beta = c(.2, .2, .6, .6), alpha = c(1, 1, 0, 0, 1, 1)),
  `4` = list(beta = c(.6, .6, .6, .6), alpha = c(1, 1, 0, 0, 1.8, 1.8)))
pad <- function(v) c(v, rep(0, P - length(v)))
expit <- function(x) exp(x) / (1 + exp(x))

if (!is.null(EXPORT_DIR)) dir.create(EXPORT_DIR, showWarnings = FALSE)
rows <- list(); rsel <- list()
for (sc in SCENARIOS) {
  beta_v <- COEF_MULT * pad(scen[[as.character(sc)]]$beta)
  alpha_v <- COEF_MULT * pad(scen[[as.character(sc)]]$alpha)
  for (rho in RHO_LIST) {
    for (itr in 1:ITR_MAX) {
      set.seed(itr)
      Sigma <- matrix(rho, P, P); diag(Sigma) <- 1
      X <- mvrnorm(n = N, mu = rep(0, P), Sigma = Sigma)
      A <- as.numeric(runif(N) < expit(as.numeric(X %*% alpha_v)))
      Y <- as.numeric(X %*% beta_v) + rnorm(N, sd = 1) + 0.5 * A
      X <- scale(X); colnames(X) <- paste0("X", 1:P)
      key <- sprintf("mult%g_sc%d_rho%g_itr%d", COEF_MULT, sc, rho, itr)
      if (!is.null(EXPORT_DIR))
        write.csv(data.frame(X, A = A, Y = Y), file.path(EXPORT_DIR, paste0("data_", key, ".csv")), row.names = FALSE)

      att_target <- att_matching(X, Y, A, mode = "paper", selected = 1:4)
      sels <- list(
        `Enh-ELRT`  = three_stage(X, Y, A, prior = "LR",  smoother = "tanh",    mode = "paper")$selected,
        `Enh-ESVMS` = three_stage(X, Y, A, prior = "SVM", smoother = "sigmoid", mode = "paper")$selected,
        `ONet`      = oaenet(X, Y, A, mode = "paper")$selected)
      for (m in names(sels)) {
        sel <- sels[[m]]
        att <- att_matching(X, Y, A, mode = "paper", selected = sel)
        r <- data.frame(scenario = sc, rho = rho, seed = itr, method = m, n_selected = length(sel),
                        att = att, att_target = att_target, pct_bias = 100 * abs(att - att_target) / 0.5)
        for (j in 1:P) r[[paste0("X", j)]] <- as.integer(j %in% sel)
        rows[[length(rows) + 1]] <- r
        rsel[[length(rsel) + 1]] <- data.frame(dataset = key, method = m, selected = paste(sel, collapse = ";"))
      }
    }
    cat(sprintf("scenario %d rho %g done\n", sc, rho))
    write.csv(do.call(rbind, rows), OUT_FILE, row.names = FALSE)
    if (!is.null(EXPORT_DIR)) write.csv(do.call(rbind, rsel), file.path(EXPORT_DIR, "r_selections.csv"), row.names = FALSE)
  }
}

res <- do.call(rbind, rows)
res$XC <- rowMeans(res[, c("X1", "X2")]); res$XP <- rowMeans(res[, c("X3", "X4")])
res$XT <- rowMeans(res[, c("X5", "X6")]); res$XS <- rowMeans(res[, paste0("X", 7:P)])
f <- function(x) round(mean(x, na.rm = TRUE), 3)
print(aggregate(cbind(XC, XP, XT, XS, pct_bias) ~ scenario + rho + method, data = res, FUN = f, na.action = na.pass), row.names = FALSE)
cat("\nPooled over rho (Figure 3 / Figure 4 layout):\n")
print(aggregate(cbind(XC, XP, XT, XS, pct_bias) ~ scenario + method, data = res, FUN = f, na.action = na.pass), row.names = FALSE)
