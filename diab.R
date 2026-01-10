library(caret)
library(data.table)
library(caTools)
library(dplyr)
library("msaenet")
library(randomForest)
library(pROC)
library(e1071)
library(boot)
library(rlearner)
library(glmnet)
library(BCEE)
library(ctmle)
library(pracma)

path = "~/RProjects/FSReal/diab.csv"

df_ori = read.csv(path)
sample = 5000
itr_max = 30

var_select<- zeros(5, ncol(df_ori))
att <- zeros(6, itr_max)

for (i in 1:itr_max){
  set.seed(i)
  
  df <- df_ori %>% slice_sample(n = sample, replace = TRUE)
  
  df$Diabetes_012 <- ifelse(df$Diabetes_012 == 0, 0, 1)
  
  dummy <- dummyVars(" ~ .", data = df)
  df <- data.frame(predict(dummy, newdata = df))
  treat_index = (which(names(df) == "HighBP"))
  
  vec_sigmoid <- function(pen) {
    pen <- 1/(1+exp(-pen))
  }
  
  vec_tanh <- function(pen){
    pen <- tanh(pen)
  }
  
  # treat & outcome
  a <- df$HighBP
  y <- df$Diabetes_012
  x <- subset(df, select = -c(Diabetes_012, HighBP))
  
  
  

  ### ESVM-S
  # fit1 = glm(a ~.,family=binomial(link='logit'),data=x)
  fit1 = svm(a ~.,
             data = x,
             type = 'C-classification',
             kernel = 'linear', scale = TRUE)
  
  weight_logit = coef(fit1)
  weight_logit = weight_logit[-c(1)] # exclude the intercept
  
  factor = 0.5
  penalty = abs(weight_logit)
  penalty[is.na(penalty)] <- 10000000000
  
  penalty <- vec_sigmoid(penalty)^factor
  # penalty <- vec_tanh(penalty)^factor
  
  # elastic net estimator
  fit2 = cv.glmnet(data.matrix(x), data.matrix(y), alpha = 0.5,
                   nfolds = 3L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                   relax = TRUE, parallel=TRUE, standardize=TRUE,
                   type.measure='mse', penalty.factor=penalty)
  # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
  beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
  
  ###fit2_mse can be optimized through gamma 1###
  # fit2_mse = fit2$cvm[fit2$lambda == fit2$lambda.1se]
  # print(fit2_mse)
  
  beta_el_allvar = as.matrix(beta_el[2:(length(beta_el))])
  
  gamma = 2
  penalty_adel = (1/abs(beta_el_allvar))^gamma
  penalty_adel[is.na(penalty_adel)] <- 10000000000
  
  # adaptive elastic net estimator
  fit3 = aenet(data.matrix(x), data.matrix(y), alphas = seq(0.05, 0.95, 0.05),
               nfolds = 5L, rule = "lambda.1se",
               penalty.factor = penalty_adel, parallel = TRUE)
  
  beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
  
  
  ###fit3_mse can be optimized through gamma###
  # fit3.pred = predict(fit3, data.matrix(x))
  # fit3_mse = msaenet.rmse(y, fit3.pred)
  # print(fit3_mse)
  
  beta_aden_allvar = as.matrix(beta_aden_raw) 
  
  beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.01)] # check
  sel.var.index = beta_aden_non_zero
  sel.var.index
  
  indices <- which(names(df) %in% names(x)[sel.var.index])
  
  for (v in indices){
    var_select[1, v] = var_select[1, v] + 1
  }
  
  print(names(df)[indices])
  # added = c(10, 11, 22, 23, 24, 26)
  
  
  # causal inference of SCC to the obesity level
  new_df = df[, c(treat_index, indices, 1)]
  cas_model = glm(as.factor(Diabetes_012) ~.,family=binomial(link='logit'), data=new_df)
  print(coef(cas_model)[2]) # 0.381
  
  
  old_df = df
  old_model = glm(as.factor(Diabetes_012) ~.,family=binomial(link='logit'), data=old_df)
  print(coef(old_model)) # 0.028 counter intuitive
  
  
  # before selection
  rlasso_model <- rlasso(as.matrix(x), a, y)
  tau_hat <- predict(rlasso_model, as.matrix(x))
  mean(tau_hat) # 0.0821
  summary(tau_hat)
  att[6, i] = att[6, i] + mean(tau_hat)
  
  # after selection
  rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
  tau_hat <- predict(rlasso_model, as.matrix(df[, c(indices)]))
  mean(tau_hat) # 0.1542
  summary(tau_hat)
  print("ESVMS")
  print(indices)
  print(mean(tau_hat)) # 0.1156
  att[1, i] = att[1, i] + mean(tau_hat)
  print(att[1, i])
  
  # comparison with OAL OAENet BCEE CTMLE
  #### OAL #####
  my_logit_mod = glm(Diabetes_012 ~ ., data = df, family = "binomial")
  odds_ratio = (coef(my_logit_mod))
  odds_ratio = odds_ratio[-c(1, treat_index)]
  gamma = 5
  
  penalty = 1/(abs(odds_ratio))^gamma
  
  penalty[is.na(penalty)] <- 10000000000
  penalty
  K.cv.fit = cv.glmnet(data.matrix(x), data.matrix(a), penalty.factor = penalty, type.measure = "mse",
                       nfolds = 5, gamma = 1,
                       relax = FALSE, family="binomial")
  fit2=glmnet(as.matrix(x), as.matrix(a), alpha =1, lambda = K.cv.fit$lambda.1se, penalty.factor = penalty, family="binomial")
  beta_ps_raw = as.matrix(coef(fit2, s = 0.02))
  beta_ps_allvar = as.matrix(beta_ps_raw[2:length(beta_ps_raw)])
  beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.1)]
  sel.var.index  = beta_ps_non_zero
  
  indices <- which(names(df) %in% names(x)[sel.var.index])
  
  for (v in indices){
    var_select[2, v] = var_select[2, v] + 1
  }
  
  
  if (length(sel.var.index)!= 0 & length(sel.var.index)!= 1){
    rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
    tau_hat <- predict(rlasso_model, as.matrix(df[, c(indices)]))
    mean(tau_hat)
    att[2, i] = att[2, i] + mean(tau_hat)
    print("OAL")
    print(indices)
    print(mean(tau_hat))
  }
  
  #### OAENet #####
  my_logit_mod = glm(Diabetes_012 ~ ., data = df, family = "binomial")
  odds_ratio = (coef(my_logit_mod))
  odds_ratio = odds_ratio[-c(1, treat_index)]
  gamma = 1
  
  penalty = 1/(abs(odds_ratio))^gamma
  
  penalty[is.na(penalty)] <- 10000000000
  penalty
  K.cv.fit = cv.glmnet(data.matrix(x), data.matrix(a),  type.measure = "mse",
                       nfolds = 5, gamma = c(0, 0.25, 0.5, 0.75, 1), relax = TRUE, family="binomial")
  
  cv.alpha = K.cv.fit$relaxed
  alpha.opt = cv.alpha$gamma.1se
  fit2=glmnet(as.matrix(x), as.matrix(a), alpha =alpha.opt, lambda = K.cv.fit$lambda.1se, penalty.factor = penalty, family="binomial")
  beta_ps_raw = as.matrix(coef(fit2, s = 0.02))
  beta_ps_allvar = as.matrix(beta_ps_raw[2:length(beta_ps_raw)])
  beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.1)]
  sel.var.index  = beta_ps_non_zero
  
  indices <- which(names(df) %in% names(x)[sel.var.index])
  
  for (v in indices){
    var_select[3, v] = var_select[3, v] + 1
  }
  
  rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
  tau_hat <- predict(rlasso_model, as.matrix(df[, c(indices)]))
  
  print("OAENet")
  print(indices)
  print(mean(tau_hat)) # 0.1033
  att[3, i] = att[3, i] + mean(tau_hat)
  
  
  #### BCEE #####
  bcee <- try(GBCEE(as.matrix(a), as.matrix(y), as.matrix(x), omega = 300*sqrt(nrow(df)),niter = 5000, 
                    family.X = "binomial", family.Y = "gaussian",
                    X1 = 1, X0 = 0, priorX = NA, priorY = NA, maxsize = NA, OR = 20,
                    truncation = c(0.01, 0.99), var.comp = "asymptotic", B = 200))
  if ('try-error' %in% class(bcee)){}else{
    new_df1<-data.frame(bcee$models.X)
    new_df2<-data.frame(bcee$models.Y)
    new_df<-rbind(new_df1[1,1:dim(x)[2]],new_df2[1,1:dim(x)[2]])
    
    sel.var.index = which(new_df[2,] > 0)
    sel.var.index 
    
    indices <- which(names(df) %in% names(x)[sel.var.index])
    
    for (v in indices){
      var_select[4, v] = var_select[4, v] + 1
    }
    
    rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
    tau_hat <- predict(rlasso_model, as.matrix(df[, c(indices)]))
    print("BCEE")
    print(indices)
    print(mean(tau_hat)) # 0.08115329
    att[4, i] = att[4, i] + mean(tau_hat)
  }
  
  
  
  #### CTMLE #####
  Q <- cbind(rep(mean(y[a == 0]), nrow(df)), rep(mean(y[a == 1]), nrow(df)))
  fit1 <- ctmleDiscrete(Y = y, A = a, W = x, Q = Q, preOrder = TRUE, detailed = TRUE)
  ctmle.cov = fit1$candidates$terms[-c(1)]
  sel.var.index = which(colnames(x)==ctmle.cov)
  sel.var.index
  
  indices <- which(names(df) %in% names(x)[sel.var.index])
  
  for (v in indices){
    var_select[5, v] = var_select[5, v] + 1
  }
  
  rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
  tau_hat <- predict(rlasso_model, as.matrix(df[, c(indices)]))
  print("CTMLE")
  print(indices)
  print(mean(tau_hat))
  att[5, i] = att[5, i] + mean(tau_hat)
}
write.csv(data.frame(var_select), "~/Desktop/dia.csv")
write.csv(data.frame(att / itr_max), "~/Desktop/dia_att.csv")

