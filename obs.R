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

path = "~/RProjects/FSReal/obs.csv"

df_ori = read.csv(path)
df = copy(df_ori)



# treatment and outcome
df$SCC <- ifelse(df$SCC == c("yes"), 1, 0) # monitor calories daily
df$NObeyesdad <- ifelse(df$NObeyesdad == c("Normal_Weight"), 1, 0)


# dummy variables
df$SMOKE <- ifelse(df$SMOKE == c("yes"), 1, 0)
df$Gender <- ifelse(df$Gender == c("Male"), 1, 0)
df$family_history_with_overweight <- ifelse(df$family_history_with_overweight == c("yes"), 1, 0)
df$FAVC <- ifelse(df$FAVC == c("yes"), 1, 0)

dummy <- dummyVars(" ~ .", data = df)
df <- data.frame(predict(dummy, newdata = df))
treat_index = (which(names(df) == "SCC"))

vec_sigmoid <- function(pen) {
  pen <- 1/(1+exp(-pen))
}

vec_tanh <- function(pen){
  pen <- tanh(pen)
}


# treat & outcome
a <- df$SCC
y <- df$NObeyesdad
x <- subset(df, select = -c(SCC, NObeyesdad))

var_select<- zeros(5, ncol(df))


for (i in 1:1){
  ### ESVM-S
  # fit1 = glm(a ~.,family=binomial(link='logit'),data=x)
  fit1 = svm(a ~.,
            data = x,
            type = 'C-classification',
            kernel = 'linear', scale = TRUE)
  
  weight_logit = coef(fit1)
  weight_logit = weight_logit[-c(1)] # exclude the intercept
  
  factor = 1
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
  
  gamma = 1
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
  
  print(names(x)[sel.var.index])
  # added = c(10, 11, 22, 23, 24, 26)
  
  
  # causal inference of SCC to the obesity level
  new_df = df[, c(treat_index, indices, ncol(df))]
  cas_model = glm(as.factor(NObeyesdad) ~.,family=binomial(link='logit'), data=new_df)
  print(coef(cas_model)[2]) # 0.381
  
  
  old_df = df
  old_model = glm(as.factor(NObeyesdad) ~.,family=binomial(link='logit'), data=old_df)
  print(coef(old_model)) # 0.028 counter intuitive
  
  
  # before selection
  rlasso_model <- rlasso(as.matrix(x), a, y)
  tau_hat_before <- predict(rlasso_model, as.matrix(x))
  mean(tau_hat) # 0.046
  summary(tau_hat_before)
  hist(tau_hat_before, main = "Before CATE Distribution")
  
  # after selection
  rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
  tau_hat_ESVMS <- predict(rlasso_model, as.matrix(df[, c(sel.var.index)]))
  mean(tau_hat_ESVMS) # 0.1542
  summary(tau_hat_ESVMS)
  hist(tau_hat_ESVMS, main = "After CATE Distribution")
  print("ESVMS")
  print(indices)
  print(mean(tau_hat_ESVMS))
  
  # comparison with OAL OAENet BCEE CTMLE
  #### OAL #####
  my_logit_mod = glm(NObeyesdad ~ ., data = df, family = "binomial")
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
  beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.05)]
  sel.var.index  = beta_ps_non_zero
  
  print(names(df)[sel.var.index])
  indices <- which(names(df) %in% names(x)[sel.var.index])
  
  for (v in indices){
    var_select[2, v] = var_select[2, v] + 1
  }
  
  if (length(indices )!= 0){
    rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
    tau_hat_OAL <- predict(rlasso_model, as.matrix(df[, c(indices )]))
    mean(tau_hat_OAL)
    print("OAL")
    print(indices)
    print(mean(tau_hat_OAL))
  }
  
  #### OAENet #####
  my_logit_mod = glm(NObeyesdad ~ ., data = df, family = "binomial")
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
  beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.01)]
  sel.var.index  = beta_ps_non_zero
  
  
  print(names(df)[sel.var.index])
  indices <- which(names(df) %in% names(x)[sel.var.index])
  for (v in indices){
    var_select[3, v] = var_select[3, v] + 1
  }
  
  rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
  tau_hat_OAENet <- predict(rlasso_model, as.matrix(df[, c(indices)]))
  
  print("OAENet")
  print(indices)
  print(mean(tau_hat_OAENet))
  
  
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
    
    print(names(df)[sel.var.index])
    indices <- which(names(df) %in% names(x)[sel.var.index])
    for (v in indices){
      var_select[4, v] = var_select[4, v] + 1
    }
    
    rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
    tau_hat_BCEE <- predict(rlasso_model, as.matrix(df[, c(indices)]))
    print("BCEE")
    print(indices)
    print(mean(tau_hat_BCEE))
  }
  
  
  
  #### CTMLE #####
  Q <- cbind(rep(mean(y[a == 0]), nrow(df)), rep(mean(y[a == 1]), nrow(df)))
  fit1 <- ctmleDiscrete(Y = y, A = a, W = x, Q = Q, preOrder = TRUE, detailed = TRUE)
  ctmle.cov = fit1$candidates$terms[-c(1)]
  sel.var.index = which(colnames(x)==ctmle.cov)
  sel.var.index
  print(names(df)[sel.var.index])
  indices <- which(names(df) %in% names(x)[sel.var.index])
  
  rlasso_model <- rlasso(as.matrix(df[, c(indices)]), a, y)
  tau_hat_CTMLE <- predict(rlasso_model, as.matrix(df[, c(indices)]))
  print("CTMLE")
  print(indices)
  print(mean(tau_hat_CTMLE))
  
  
  for (v in indices){
    var_select[5, v] = var_select[5, v] + 1
  }

}
write.csv(data.frame(var_select), "~/s.csv")

df_att = data.frame(cbind(tau_hat_before,tau_hat_ESVMS, tau_hat_OAENet, tau_hat_CTMLE))
write.csv(df_att,"~/Desktop/df_obs.csv")

# # LR no selection
# train_index <- sample(1:nrow(df), 0.6 * nrow(df))
# train_data <- df[train_index, ]
# test_data <- df[-train_index, ]

# logit_model = glm(as.factor(NObeyesdad) ~.,family=binomial(link='logit'), data=train_data)
# 
# lr_pred_prob <- predict(logit_model, test_data, type = "response")
# lr_predictions <- ifelse(lr_pred_prob > 0.5, 1, 0)
# 
# lr_conf_matrix <- table(test_data$NObeyesdad, lr_predictions)
# roc_lr_ns <- roc(as.numeric(test_data$NObeyesdad), as.numeric(lr_predictions))
# 
# print("lr_conf_matrix_ns")
# print(lr_conf_matrix)

# # RF no selection
# rf_model_ns <- randomForest(as.factor(NObeyesdad) ~ .,
#                             data = train_data,
#                             ntree = 1000,
#                             mtry = 3,
#                             importance = TRUE)
# rf_pred <- predict(rf_model_ns, test_data)
# rf_conf_matrix <- table(test_data$NObeyesdad, rf_pred)
# print("rf_conf_matrix_ns")
# print(rf_conf_matrix)
# 
# # LR selection
# df_selected = df[sel.var.index]
# train_data <- df[train_index, c(sel.var.index, treat_index, ncol(df))]
# test_data <- df[-train_index, c(sel.var.index, treat_index, ncol(df))]
# 
# logit_model = glm(as.factor(NObeyesdad) ~.,family=binomial(link='logit'), data=train_data)
# 
# lr_pred_prob <- predict(logit_model, test_data, type = "response")
# lr_predictions <- ifelse(lr_pred_prob > 0.5, 1, 0)
# 
# lr_conf_matrix <- table(test_data$NObeyesdad, lr_predictions)
# roc_lr_s <- roc(as.numeric(test_data$NObeyesdad), as.numeric(lr_predictions))
# 
# print("lr_conf_matrix_s")
# print(lr_conf_matrix)
# 
# 
# # RF selection
# rf_model_s <- randomForest(as.factor(NObeyesdad) ~ .,
#                             data = train_data,
#                             ntree = 1000,
#                             mtry = 3,
#                             importance = TRUE)
# rf_pred <- predict(rf_model_s, test_data)
# rf_conf_matrix <- table(test_data$NObeyesdad, rf_pred)
# 
# print("rf_conf_matrix_s")
# print(rf_conf_matrix)

# plot(roc_lr_ns, col = "blue", lwd = 2, main = "ROC curve")
# plot(roc_lr_s, col = "red", lwd = 2, add = TRUE)

