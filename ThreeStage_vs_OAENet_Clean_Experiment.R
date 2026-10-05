
################################################################################
# Clean comparison experiment: Published Three-Stage Framework vs OAENet
#
# Methods:
#   Enh-ELRT  : exposure logistic regression -> tanh penalty smoothing
#               -> outcome elastic net -> adaptive elastic net
#   Enh-ESVMS : exposure linear SVM -> sigmoid penalty smoothing
#               -> outcome elastic net -> adaptive elastic net
#   OAENet    : outcome OLS -> outcome-adaptive weights
#               -> adaptive elastic-net treatment model
#   Target    : oracle X_C U X_P = {X1,X2,X3,X4}
#
# Published main synthetic design:
#   p = 20
#   N = {200,500,1000}
#   rho = {0,.25,.5,.75}
#   scenarios = 1:4
#   seeds = 1:30
#   true treatment effect = 0.5
#
# NOTE 1:
# The published treatment DGP is typeset with epsilon_treat ~ N(0,1)
# and an indicator involving expit(X beta), which is unusual.
# Default here is the standard Bernoulli-logit interpretation:
#     T ~ Bernoulli(expit(X beta))
# Set treatment_mode="paper_literal" to run the printed equation literally.
#
# NOTE 2:
# The papers parameterize L1 and L2 separately. The exact glmnet alpha
# is not stated in the paper text, so alpha_enet=0.5 is configurable.
################################################################################

packages <- c("glmnet","MatchIt","MASS","e1071")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) install.packages(missing)

library(glmnet)
library(MatchIt)
library(MASS)
library(e1071)

CONFIG <- list(
  p=20,
  true_te=.5,
  alpha_enet=.5,
  nfolds=5,
  gamma1_elrt=.5,
  gamma1_esvms=1,
  gamma2=1,
  coef_tol=1e-8,
  treatment_mode="bernoulli_logit",
  match_caliper=NULL
)

make_sigma <- function(p,rho) {
  S <- matrix(rho,p,p)
  diag(S) <- 1
  S
}

scenario_coefficients <- function(scenario,p=20) {
  theta <- rep(0,p); beta <- rep(0,p)
  if (scenario==1) {
    theta[1:4] <- c(.6,.6,.6,.6)
    beta[c(1,2,5,6)] <- c(1,1,1,1)
  } else if (scenario==2) {
    theta[1:4] <- c(.6,.6,.6,.6)
    beta[c(1,2,5,6)] <- c(.4,.4,1,1)
  } else if (scenario==3) {
    theta[1:4] <- c(.2,.2,.6,.6)
    beta[c(1,2,5,6)] <- c(1,1,1,1)
  } else if (scenario==4) {
    theta[1:4] <- c(.6,.6,.6,.6)
    beta[c(1,2,5,6)] <- c(1,1,1.8,1.8)
  } else stop("scenario must be 1,2,3,4")
  list(theta=theta,beta=beta)
}

simulate_scenario <- function(n,rho,scenario,seed,
                              treatment_mode=CONFIG$treatment_mode) {
  set.seed(seed)
  pars <- scenario_coefficients(scenario,CONFIG$p)
  X <- MASS::mvrnorm(
    n=n, mu=rep(0,CONFIG$p),
    Sigma=make_sigma(CONFIG$p,rho)
  )
  colnames(X) <- paste0("X",1:CONFIG$p)
  eta <- as.vector(X %*% pars$beta)

  if (treatment_mode=="bernoulli_logit") {
    T <- rbinom(n,1,plogis(eta))
  } else if (treatment_mode=="paper_literal") {
    eps_treat <- rnorm(n)
    T <- as.integer(plogis(eta) < eps_treat)
  } else stop("unknown treatment_mode")

  Y <- CONFIG$true_te*T + as.vector(X %*% pars$theta) + rnorm(n)
  data.frame(Y=Y,T=T,X,check.names=FALSE)
}

safe_inverse_power <- function(x,gamma,eps=1e-8,cap=1e12) {
  z <- pmax(abs(x),eps)^(-gamma)
  z[!is.finite(z)] <- cap
  pmin(z,cap)
}

sigmoid_penalty <- function(beta,gamma=1)
  (1/(1+exp(-abs(beta))))^gamma

tanh_penalty <- function(beta,gamma=.5)
  pmax(tanh(abs(beta)),1e-8)^gamma

exposure_logistic_coef <- function(X,T) {
  fit <- glm(T~.,data=data.frame(T=T,X),family=binomial())
  b <- coef(fit)[-1]
  b[!is.finite(b)] <- 0
  as.numeric(b)
}

exposure_svm_coef <- function(X,T) {
  fit <- e1071::svm(
    x=as.matrix(X), y=as.factor(T),
    kernel="linear", scale=FALSE, type="C-classification"
  )
  as.numeric(drop(t(fit$coefs) %*% fit$SV))
}

weighted_gaussian_enet <- function(X,Y,penalty_factor,
                                   alpha=CONFIG$alpha_enet) {
  fit <- cv.glmnet(
    as.matrix(X),Y,
    family="gaussian",
    alpha=alpha,
    penalty.factor=penalty_factor,
    nfolds=CONFIG$nfolds,
    standardize=TRUE
  )
  b <- as.matrix(coef(fit,s="lambda.1se"))[-1,1]
  list(beta=b,fit=fit)
}

enhanced_three_stage <- function(
    Data,
    prior=c("LR","SVM"),
    smoother=c("tanh","sigmoid"),
    gamma1=NULL,
    gamma2=CONFIG$gamma2
) {
  prior <- match.arg(prior)
  smoother <- match.arg(smoother)
  X <- as.matrix(Data[,grep("^X",names(Data)),drop=FALSE])
  Y <- Data$Y; T <- Data$T
  vnames <- colnames(X)

  beta_prior <- if (prior=="LR") {
    exposure_logistic_coef(X,T)
  } else {
    exposure_svm_coef(X,T)
  }

  if (is.null(gamma1))
    gamma1 <- if (smoother=="tanh") CONFIG$gamma1_elrt else CONFIG$gamma1_esvms

  w1 <- if (smoother=="tanh")
    tanh_penalty(beta_prior,gamma1)
  else
    sigmoid_penalty(beta_prior,gamma1)

  # Stage 2: outcome elastic net informed by exposure model
  f2 <- weighted_gaussian_enet(X,Y,w1)
  theta <- f2$beta

  # Stage 3: adaptive weights from Stage-2 elastic-net coefficients
  psi <- safe_inverse_power(theta,gamma2)

  # Final adaptive elastic net
  f3 <- weighted_gaussian_enet(X,Y,psi)
  final <- f3$beta

  selected <- vnames[abs(final)>CONFIG$coef_tol]

  list(
    selected=selected,
    beta_prior=setNames(beta_prior,vnames),
    penalty_stage2=setNames(w1,vnames),
    theta_stage2=setNames(theta,vnames),
    penalty_stage3=setNames(psi,vnames),
    theta_final=setNames(final,vnames)
  )
}

fit_enh_elrt <- function(Data)
  enhanced_three_stage(Data,"LR","tanh",
                       CONFIG$gamma1_elrt,CONFIG$gamma2)

fit_enh_esvms <- function(Data)
  enhanced_three_stage(Data,"SVM","sigmoid",
                       CONFIG$gamma1_esvms,CONFIG$gamma2)

fit_oaenet <- function(Data,gamma=3) {
  X <- as.matrix(Data[,grep("^X",names(Data)),drop=FALSE])
  Y <- Data$Y; T <- Data$T
  vnames <- colnames(X)
  Xs <- scale(X)

  # OAENet Stage 1: OLS outcome coefficients
  fy <- lm.fit(cbind(1,Xs),Y)
  theta <- fy$coefficients[-1]
  theta[!is.finite(theta)] <- 0

  w <- safe_inverse_power(theta,gamma)

  # OAENet Stage 2: adaptive elastic-net treatment model
  ft <- cv.glmnet(
    Xs,T,
    family="binomial",
    alpha=CONFIG$alpha_enet,
    penalty.factor=w,
    nfolds=CONFIG$nfolds,
    standardize=FALSE
  )
  b <- as.matrix(coef(ft,s="lambda.1se"))[-1,1]
  selected <- vnames[abs(b)>CONFIG$coef_tol]

  list(
    selected=selected,
    theta_ols=setNames(theta,vnames),
    penalty=setNames(w,vnames),
    beta_treatment=setNames(b,vnames)
  )
}

estimate_att <- function(Data,selected_vars) {
  if (length(selected_vars)==0) return(NA_real_)

  form <- reformulate(selected_vars,response="T")
  args <- list(
    formula=form,data=Data,method="nearest",
    distance="glm",ratio=1,estimand="ATT",replace=FALSE
  )
  if (!is.null(CONFIG$match_caliper))
    args$caliper <- CONFIG$match_caliper

  mm <- suppressWarnings(do.call(MatchIt::matchit,args))
  md <- MatchIt::match.data(mm)
  if (nrow(md)==0 || length(unique(md$T))<2) return(NA_real_)

  # Published experiment: matching, then regression on matched data
  f <- reformulate(c("T",selected_vars),response="Y")
  unname(coef(lm(f,data=md))["T"])
}

run_one_comparison <- function(
    n=500,rho=.25,scenario=1,seed=1,
    treatment_mode=CONFIG$treatment_mode
) {
  D <- simulate_scenario(n,rho,scenario,seed,treatment_mode)

  fits <- list(
    Enh_ELRT=fit_enh_elrt(D),
    Enh_ESVMS=fit_enh_esvms(D),
    OAENet=fit_oaenet(D)
  )

  target <- paste0("X",1:4)
  att_target <- estimate_att(D,target)

  ans <- lapply(names(fits),function(m) {
    s <- fits[[m]]$selected
    att <- estimate_att(D,s)
    data.frame(
      method=m,n=n,rho=rho,scenario=scenario,seed=seed,
      n_selected=length(s),
      selected=paste(s,collapse=","),
      att=att,
      target_att=att_target,
      abs_selection_bias=abs(att-att_target),
      pct_selection_bias=100*abs(att-att_target)/CONFIG$true_te
    )
  })

  ans[[length(ans)+1]] <- data.frame(
    method="Target",n=n,rho=rho,scenario=scenario,seed=seed,
    n_selected=4,selected=paste(target,collapse=","),
    att=att_target,target_att=att_target,
    abs_selection_bias=0,pct_selection_bias=0
  )
  do.call(rbind,ans)
}

run_full_experiment <- function(
    Ns=c(200,500,1000),
    rhos=c(0,.25,.5,.75),
    scenarios=1:4,
    seeds=1:30,
    treatment_mode=CONFIG$treatment_mode,
    verbose=TRUE
) {
  grid <- expand.grid(
    scenario=scenarios,n=Ns,rho=rhos,seed=seeds,
    KEEP.OUT.ATTRS=FALSE
  )
  out <- vector("list",nrow(grid))

  for (i in seq_len(nrow(grid))) {
    g <- grid[i,]
    if (verbose)
      cat(sprintf("[%d/%d] s=%d n=%d rho=%.2f seed=%d\n",
                  i,nrow(grid),g$scenario,g$n,g$rho,g$seed))
    out[[i]] <- tryCatch(
      run_one_comparison(
        g$n,g$rho,g$scenario,g$seed,treatment_mode
      ),
      error=function(e) {
        data.frame(
          method=c("Enh_ELRT","Enh_ESVMS","OAENet","Target"),
          n=g$n,rho=g$rho,scenario=g$scenario,seed=g$seed,
          n_selected=NA,selected=NA,att=NA,target_att=NA,
          abs_selection_bias=NA,pct_selection_bias=NA
        )
      }
    )
  }
  do.call(rbind,out)
}

summarize_results <- function(results) {
  aggregate(
    cbind(att,abs_selection_bias,pct_selection_bias,n_selected) ~
      method+scenario+n+rho,
    data=results,
    FUN=function(x)c(mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE))
  )
}

selection_frequency <- function(
    n=500,rho=.25,scenario=1,seeds=1:30,
    treatment_mode=CONFIG$treatment_mode
) {
  vars <- paste0("X",1:CONFIG$p)
  methods <- c("Enh_ELRT","Enh_ESVMS","OAENet")
  cnt <- matrix(0,length(methods),length(vars),
                dimnames=list(methods,vars))

  for (s in seeds) {
    D <- simulate_scenario(n,rho,scenario,s,treatment_mode)
    fits <- list(
      Enh_ELRT=fit_enh_elrt(D),
      Enh_ESVMS=fit_enh_esvms(D),
      OAENet=fit_oaenet(D)
    )
    for (m in methods)
      cnt[m,fits[[m]]$selected] <- cnt[m,fits[[m]]$selected]+1
  }
  as.data.frame(cnt/length(seeds))
}

################################################################################
# QUICK START
################################################################################

# one <- run_one_comparison(n=500,rho=.25,scenario=1,seed=1)
# print(one)

# sf <- selection_frequency(n=500,rho=.25,scenario=1)
# print(sf)

# Full published grid = 4*3*4*30 = 1440 generated datasets:
# full <- run_full_experiment()
# write.csv(full,"three_stage_vs_oaenet_full_results.csv",row.names=FALSE)
# summary <- summarize_results(full)
# write.csv(summary,"three_stage_vs_oaenet_summary.csv",row.names=FALSE)
################################################################################
