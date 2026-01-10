# install.packages("fastDummies")
library(fastDummies)
library(msaenet)
library(glmnet)
library(MatchIt)
library(lmtest)
library(phonTools)
library(MASS)
library(sandwich)
library(CovSel)
library(BCEE)
library(mlbench)
library(Boruta)
library(bacr)
library(mltools)
library(data.table)
library("matrixStats")
library(dplyr)
library(tictoc)
library(boot)
library(e1071)
library(caTools)
library(caret)
library(interp)
library(WeightIt)
library(tictoc)
library(CovSel)
library(ctmle)
library(BART)
library(tmle)
library(SuperLearner)
library(caret)
library(splines)
library(grf)
library(rpart)
library(FLAME)
library(tictoc)

######### general settings ########
# define path
path = "~/RProjects/FSProject_Extra/"

#Set the number of simulation
itr_max = 30
# sample size
n = 1000
# total number of predictors
p = 20

######### define scenarios #########
# scenario = 0 # can be 0 ~ 4
for (scenario in c(0:4)){

# set information for simulating coviariates
mean_x = 0 
sig_x = 1
rho_list = c(0, 0.25, 0.5, 0.75)
# rho_list = c(0.5, 0.75)

#### detail descriptions of scenarios ######
if (scenario == 0){ 
  gd.var = 9 # number of good variables (consistent with beta_v and alpha_v)

  # Set strength of relationship between covariates and outcome
  beta_v =  2*c( 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0, 0, 0, rep(0, p - gd.var) ) 
  # Set strength of relationship between covariates and treatment
  alpha_v = 2*c( 1.0, 1.0, 1.0, 0, 0, 0, 1.0, 1.0, 1.0, rep(0, p - gd.var) )
  pC = 3 # confounders 1:3
  pI = 3 # treatment predictors 7:9
  pP = 3 # outcome predictors 4:6
}
if (scenario == 1){ # normal
  gd.var = 7
  beta_v =  2*c( 0.6, 0.6, 0.6, 0.6, 0, 0, 0, rep(0, p - gd.var) )
  alpha_v = 2*c( 1.0, 1.0, 0, 0, 1.0, 1.0, 0, rep(0, p - gd.var) )
  pC = 2
  pI = 2 
  pP = 2
}
if (scenario == 2){ # weak treat pred
  gd.var = 7
  beta_v =  2*c( 0.6, 0.6, 0.6, 0.6, 0, 0, 0, rep(0, p - gd.var) )
  alpha_v = 2*c( 0.4, 0.4, 0, 0, 1.0, 1.0, 0, rep(0, p - gd.var) )
  pC = 2
  pI = 2 
  pP = 2
}
if (scenario == 3){ # weak confounders
  gd.var = 7
  beta_v =  2*c( 0.2, 0.2, 0.6, 0.6, 0, 0, 0, rep(0, p - gd.var) )
  alpha_v = 2*c( 1.0, 1.0, 0, 0, 1.0, 1.0, 0, rep(0, p - gd.var) )
  pC = 2
  pI = 2 
  pP = 2
}
if (scenario == 4){ # strong treat pred
  gd.var = 7
  beta_v =  2*c( 0.6,0.6,0.6,0.6, 0, 0, 0, rep(0, p - gd.var) )
  alpha_v = 2*c( 1.0, 1.0, 0, 0, 1.8, 1.8, 0, rep(0, p - gd.var) ) 
  pC = 2
  pI = 2 
  pP = 2
}


########## create target list ###########
pS = p - (pC+pI+pP)
#all variable list
var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))
#Some other list we will use
target.var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""))
potconf.var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""))
conf.var.list = c(paste("Xc",1:pC,sep=""))

names(beta_v) = names(alpha_v) = var.list
### set true average treatment effect
bA = 0.5

### define some functions for generating data, ATE estimates, and the wAMD,
expit = function(x){ 
  pr = ( exp(x) / (1+exp(x)) ) 
  return(pr)
}

vec_norm = function(x){
  x / sqrt(sum(x^2))
}

######## define att and var selection matrix##########

### method_list
method_list = c("BART", "BART_FS", 
                "TMLE_SL","TMLE_SL_FS", 
                "CaFo", "CaFo_FS",
                "XL", "XL_FS",
                "RL", "RL_FS",
                "FLAME", "FLAME_FS")

#declare a zero matrix to hold att
att.target = zeros(itr_max)
att.BART = zeros(itr_max)
att.BART_FS = zeros(itr_max)
att.TMLE_SL = zeros(itr_max)
att.TMLE_SL_FS = zeros(itr_max)
att.CaFo = zeros(itr_max)
att.CaFo_FS = zeros(itr_max)
att.XL = zeros(itr_max)
att.XL_FS = zeros(itr_max)
att.RL = zeros(itr_max)
att.RL_FS = zeros(itr_max)
att.FLAME = zeros(itr_max)
att.FLAME_FS = zeros(itr_max)

#Time Elaspsed
time.BART = zeros(itr_max)
time.BART_FS = zeros(itr_max)
time.TMLE_SL = zeros(itr_max)
time.TMLE_SL_FS = zeros(itr_max)
time.CaFo = zeros(itr_max)
time.CaFo_FS = zeros(itr_max)
time.XL = zeros(itr_max)
time.XL_FS = zeros(itr_max)
time.RL = zeros(itr_max)
time.RL_FS = zeros(itr_max)
time.FLAME = zeros(itr_max)
time.FLAME_FS = zeros(itr_max)

####### Penalty Smooth Function ##########
vec_sigmoid <- function(pen) {
  pen <- 1/(1+exp(-pen))
}

vec_tanh <- function(pen){
  pen <- tanh(pen)
}

###### Running the simulation ###########
for (k in 1:length(rho_list)){
  rho = as.numeric(rho_list[k])
  
  # Running simulation for each time
  for (i in 1:itr_max) {
    itr = i
    set.seed(itr) # random seed set to agree with itr
    label = 0 # method label
    
    print(itr)
    
    #### Simulated data ####
    Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list)) 
    diag(Sigma_x) = sig_x^2
    Mean_x = rep(mean_x,length(var.list))
    Data = as.data.frame(mvrnorm(n = n,mu=Mean_x,Sigma = Sigma_x,empirical = FALSE))
    names(Data) = var.list
    
    #Creating treatment
    gA_x = rowSums(Data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))
    pA = expit( gA_x )
    Data$A = as.numeric( runif(n=length(pA)) < pA) # simulate A 
    
    #Creating outcome
    gY_xA = rowSums(Data[,var.list]*matrix(beta_v,nrow=n,ncol=length(var.list),byrow=TRUE))   
    Data$Y = gY_xA + rnorm(n=n,sd=sig_x) #adding random noise
    Data$Y = Data$Y + Data$A*bA
    
    # Normlize coviarates to have mean 0 and standard deviation 1
    temp.mean = colMeans(Data[,var.list])
    Temp.mean = matrix(temp.mean,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
    Data[,var.list] = Data[,var.list] - Temp.mean
    temp.sd = apply(Data[var.list],FUN=sd,MARGIN=2)
    Temp.sd = matrix(temp.sd,ncol=length(var.list),nrow=nrow(Data),byrow=TRUE)
    Data[var.list] = Data[,var.list] / Temp.sd
    rm(list=c("temp.mean","Temp.mean","temp.sd","Temp.sd"))
    
    X_ori=Data[,1:p]
    Y=Data[,(p+2)]
    A=Data[,(p+1)]
    
    ######## target model
    target.form = formula(paste("A~",paste(target.var.list,collapse="+")))
    mm = matchit(target.form, data = Data, method = "nearest", distance = "glm",
                 ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                 replace = FALSE, discard = "both")
    
    effect.form = formula(paste("Y~A+",paste(target.var.list,collapse="+")))
    fit.effect = lm(effect.form, data = match.data(mm))
    coeftest(fit.effect, vcov. = vcovCL)["A",,drop=FALSE]
    final.coef = fit.effect$coefficients
    att.target[itr] = final.coef[2]
    sel.var.index = c(1:(pC+pP))
    # Percent var select
    
    print("target")
    print(att.target[itr])
    
    #### Enh-lr-sigmoid ####
    ts <- Sys.time()
    
    fit1 = glm(A ~.,family = binomial(link='logit'), data = cbind(A,X_ori))
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)]
    weight_logit
    
    factor = 0.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(as.matrix(X_ori), as.matrix(Y), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     penalty.factor=penalty)
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    
    beta_el
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    fit3 = aenet(data.matrix(X_ori), data.matrix(Y), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_adnet_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    # beta_adnet_raw = as.matrix(coef(fit3, s=0.2))
    beta_adnet_allvar = as.matrix(beta_adnet_raw[1:p])
    beta_adnet_non_zero = row(beta_adnet_allvar)[which(abs(beta_adnet_allvar) >= 0.01)]
    sel.var.index = c(beta_adnet_non_zero)
    
    sel_time <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
    print(sel.var.index)
    print(sel_time)
    
    X_new <- X_ori[sel.var.index]
    
    ####### BART before ########
    ts <- Sys.time()
    
    X_treat <- X_ori[A == 1, ]
    Y_treat <- Y[A == 1]
    X_control <- X_ori[A == 0, ]
    Y_control <- Y[A == 0]
    
    bart_control <- wbart(x.train = X_control, y.train = Y_control)
    ave(bart_control$varcount[,1])[1]
    
    pred_tmp <- try(mean(predict(bart_control, X_treat)))
    Y_treat_pred <- mean(predict(bart_control, X_treat)) # use control model to predict treatment
    ATT_before_treated <- mean(Y_treat - Y_treat_pred) # treated(obs) - control(pred)
    
    att.BART[i] <- ATT_before_treated
    time.BART[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
    
    print("BART ATT")
    print(ATT_before_treated)
    print("BART time")
    print(time.BART[i])
    
    ####### BART after ########
    ts <- Sys.time()
    
    X_treat <- X_new[A == 1, ]
    X_control <- X_new[A == 0, ]
    
    bart_control <- wbart(x.train = X_control, y.train = Y_control)
    Y_treat_pred <- mean(predict(bart_control, X_treat))
    ATT_after_treated <- mean(Y_treat - Y_treat_pred)
    
    att.BART_FS[i] <- ATT_after_treated
    time.BART_FS[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs")) + sel_time
    
    print("BART FS ATT")
    print(ATT_after_treated)
    print("BART FS time")
    print(time.BART_FS[i])
    
    
    # #### TMLE with super learner before ####
    ts = Sys.time()
    
    SL.library <- c("SL.glm", "SL.glmnet", "SL.gam", "SL.mean")
    tmle_result <- tmle(Y = Y, A = A, W = X_ori,
                        family = "gaussian",
                        Q.SL.library = SL.library,
                        g.SL.library = SL.library)
    att.TMLE_SL[i] <- mean(tmle_result$estimates$ATT$CI)
    time.TMLE_SL[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
    
    print("TMLE ATT")
    print(att.TMLE_SL[i])
    
    ##### TMLE with super learner after ####
    ts = Sys.time()
    
    tmle_result <- tmle(Y = Y, A = A, W = X_new,
                        family = "gaussian",
                        Q.SL.library = SL.library,
                        g.SL.library = SL.library)
    
    att.TMLE_SL_FS[i] <- mean(tmle_result$estimates$ATT$CI)
    time.TMLE_SL_FS[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs")) + sel_time
    
    print("TMLE FS ATT")
    print(att.TMLE_SL_FS[i])
    
    ##### causal forests before #####
    ts = Sys.time()
    
    X_treat <- X_ori[A == 1, ]
    Y_treat <- Y[A == 1]
    
    cf <- causal_forest(X_ori, Y, A)
    treated_effects <- predict(cf, newdata = X_treat)$predictions
    
    att.CaFo[i] <- mean(treated_effects)
    time.CaFo[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
    
    print("CAFO ATT")
    print(att.CaFo[i])
    
    ##### causal forests after #####
    ts = Sys.time()
    
    X_treat <- X_new[A == 1, ]
    cf <- causal_forest(X_new, Y, A)
    treated_effects <- predict(cf, newdata = X_treat)$predictions
    
    att.CaFo_FS[i] <- mean(treated_effects)
    time.CaFo_FS[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs")) + sel_time
    
    print("CAFO FS ATT")
    print(att.CaFo_FS[i])
    
    
    #### X-learner before ####
    ts = Sys.time()
    
    treat <- X_ori[A == 1, ]
    control <- X_ori[A == 0, ]
    Y_treat <- Y[A == 1]
    Y_control <- Y[A == 0]
    outcome_model_treated <- regression_forest(treat, Y_treat)
    outcome_model_control <- regression_forest(control, Y_control)
    mu_hat_treated <- predict(outcome_model_treated, X_ori)$predictions
    mu_hat_control <- predict(outcome_model_control, X_ori)$predictions
    D_treated <- Y[A == 1] - mu_hat_control[A == 1]
    D_control <- mu_hat_treated[A == 0] - Y[A == 0]
    tau_treated_model <- regression_forest(treat, D_treated)
    tau_control_model <- regression_forest(control, D_control)
    tau_treated_hat <- predict(tau_treated_model, X_ori)$predictions
    tau_control_hat <- predict(tau_control_model, X_ori)$predictions
    p_hat <- mean(A)
    tau_hat <- p_hat * tau_control_hat + (1 - p_hat) * tau_treated_hat
    
    att.XL[i] <- mean(tau_hat[A == 1])
    time.XL[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
    
    print("XL ATT")
    print(att.XL[i])
    
    #### X-learner after ####
    ts = Sys.time()
    
    treat <- X_new[A == 1, ]
    control <- X_new[A == 0, ]
    Y_treat <- Y[A == 1]
    Y_control <- Y[A == 0]
    outcome_model_treated <- regression_forest(treat, Y_treat)
    outcome_model_control <- regression_forest(control, Y_control)
    mu_hat_treated <- predict(outcome_model_treated, X_new)$predictions
    mu_hat_control <- predict(outcome_model_control, X_new)$predictions
    D_treated <- Y[A == 1] - mu_hat_control[A == 1]
    D_control <- mu_hat_treated[A == 0] - Y[A == 0]
    tau_treated_model <- regression_forest(treat, D_treated)
    tau_control_model <- regression_forest(control, D_control)
    tau_treated_hat <- predict(tau_treated_model, X_new)$predictions
    tau_control_hat <- predict(tau_control_model, X_new)$predictions
    p_hat <- mean(A)
    tau_hat <- p_hat * tau_control_hat + (1 - p_hat) * tau_treated_hat
    
    att.XL_FS[i] <- mean(tau_hat[A == 1])
    time.XL_FS[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs")) + sel_time
    
    print("XL FS ATT")
    print(att.XL_FS[i])
    
    
    #### R-learner before #####
    ts = Sys.time()
    
    data = data.frame(Y, A, X_ori)
    outcome_forest <- regression_forest(X_ori, Y)
    propensity_forest <- regression_forest(X_ori, A)
    Y_hat <- predict(outcome_forest)$predictions
    Y_residuals <- Y - Y_hat
    A_hat <- predict(propensity_forest)$predictions
    A_residuals <- A - A_hat
    tau_forest <- causal_forest(X_ori, Y_residuals, A_residuals)
    tau_hat <- predict(tau_forest)$predictions
    
    att.RL[i] <- mean(tau_hat[A == 1])
    time.RL[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
    
    print("RL ATT")
    print(att.RL[i])
    
    #### R-learner after #####
    ts = Sys.time()
    
    data = data.frame(Y, A, X_new)
    outcome_forest <- regression_forest(X_new, Y)
    propensity_forest <- regression_forest(X_new, A)
    Y_hat <- predict(outcome_forest)$predictions
    Y_residuals <- Y - Y_hat
    A_hat <- predict(propensity_forest)$predictions
    A_residuals <- A - A_hat
    tau_forest <- causal_forest(X_new, Y_residuals, A_residuals)
    tau_hat <- predict(tau_forest)$predictions
    
    att.RL_FS[i] <- mean(tau_hat[A == 1])
    time.RL_FS[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs")) + sel_time
    
    print("RL FS ATT")
    print(att.RL_FS[i])
    
    ###### FLAME before #####
    ts = Sys.time()

    data <- data.frame(X_ori, A, Y)
    cat_data = copy(data)
    cat_data <- cat_data[, sapply(cat_data, function(col) length(unique(col)) > 1)]
    colnames(cat_data)[which(colnames(cat_data) == "Y")] <- 'outcome'
    colnames(cat_data)[which(colnames(cat_data) == "A")] <- 'treated'

    for (cat_col in 1:p){
      # cat_data[,cat_col] = cut(as.numeric(as.character(cat_data[,cat_col])),
      #                          breaks = c(-Inf, mean_x-3*sig_x, mean_x-2*sig_x, mean_x-sig_x, mean_x, mean_x+sig_x, mean_x+2*sig_x, mean_x+3*sig_x, Inf),
      #                          labels = c(1,2,3,4,5,6,7,8))
      cat_data[, cat_col] <- as.numeric(as.character(cat_data[, cat_col]))
      cat_data[,cat_col] <- cut(cat_data[,cat_col],
                                breaks = quantile(cat_data[,cat_col], probs = seq(0, 1, by = 0.1)),
                                labels = 0:9,include.lowest = TRUE)
    }

    flame_obj <- try(FLAME(data = cat_data, holdout = 0.1))
    
    if (inherits(flame_obj, "try-error")){
      att.FLAME[i] = NA
      time.FLAME[i] <- NA
    }else{

      att.FLAME[i] <- try(ATT(flame_obj))
      if (inherits(att.FLAME[i], "try-error")){
        att.FLAME[i] = NA
        time.FLAME[i] <- NA
      }
      time.FLAME[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs"))
  
      print("FLAME ATT")
      print(att.FLAME[i])
    }

    ###### FLAME after #####
    ts = Sys.time()

    data <- data.frame(X_new, A, Y)
    cat_data = copy(data)
    cat_data <- cat_data[, sapply(cat_data, function(col) length(unique(col)) > 1)]
    colnames(cat_data)[which(colnames(cat_data) == "Y")] <- 'outcome'
    colnames(cat_data)[which(colnames(cat_data) == "A")] <- 'treated'

    for (cat_col in 1:(length(sel.var.index))){
      # cat_data[,cat_col] = cut(as.numeric(as.character(cat_data[,cat_col])),
      #                          breaks = c(-Inf, mean_x-3*sig_x, mean_x-2*sig_x, mean_x-sig_x, mean_x, mean_x+sig_x, mean_x+2*sig_x, mean_x+3*sig_x, Inf),
      #                          labels = c(1,2,3,4,5,6,7,8))
      cat_data[,cat_col] <- cut(cat_data[,cat_col],
                                breaks = quantile(cat_data[,cat_col], probs = seq(0, 1, by = 0.1)),
                                labels = 0:9,include.lowest = TRUE)
    }

    flame_obj <- try(FLAME(data = cat_data, holdout = 0.1))

    if (inherits(flame_obj, "try-error")){
      att.FLAME_FS[i] = NA
      time.FLAME_FS[i] <- NA
    }else{
      att.FLAME_FS[i] <- try(ATT(flame_obj))
      if (inherits(att.FLAME_FS[i], "try-error")){
        att.FLAME_FS[i] = NA
        time.FLAME_FS[i] <- NA
      }
      time.FLAME_FS[i] <- as.numeric(difftime(Sys.time(), ts, units = "secs")) + sel_time
  
      print("FLAME FS ATT")
      print(att.FLAME_FS[i])
    }

  }
  
  title = paste("s_",scenario,"_att_",bA,"_rho_",rho,"_n_",n,"_p_",p,"_itr_",itr,"_comp_thrs", sep="")
  
  df_att = data.frame(att.target, att.BART, att.BART_FS, 
                      att.TMLE_SL, att.TMLE_SL_FS,
                      att.CaFo, att.CaFo_FS,
                      att.XL, att.XL_FS,
                      att.RL, att.RL_FS,
                      att.FLAME, att.FLAME_FS)
  df_time = data.frame(time.BART, time.BART_FS, 
                       time.TMLE_SL, time.TMLE_SL_FS,
                       time.CaFo, time.CaFo_FS,
                       time.XL, time.XL_FS,
                       time.RL, time.RL_FS,
                       time.FLAME, time.FLAME_FS)
  
  write.csv(df_att, paste(path, "att/", title,"att.csv",sep=""))
  write.csv(df_time, paste(path, "time/",title,"time.csv",sep=""))
}


}




