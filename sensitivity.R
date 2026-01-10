####Data are simulated in such a way that only good variables are correlated, noisy variables are not correlated
########################################
library("msaenet")
library("glmnet")
library("MatchIt")
library("lmtest")
library("phonTools")
library("MASS")
library("sandwich")
#install.packages("BCEE")
library(BCEE)
#install.packages("mlbench")
library(mlbench)
#install.packages("Boruta")
library(Boruta)
#install.packages("bacr")
library(bacr)
#library(rJava)
#library(CovSelHigh)

library(mltools)
library(data.table)
library("matrixStats")
#install.packages("dplyr")
library(dplyr)
library(tictoc)

# install.packages("boot")
library(boot)

# install.packages('e1071')
library(e1071)

# install.packages('caTools')
library(caTools)

# install.packages('CovSel')
library(CovSel)

# install.packages('ctmle')
library(ctmle)

# install.packages('dplyr')
library(dplyr)




######### general settings ########
# define path
path = "~/RProjects/FSProject_Sensitivity/"

install.packages(paste(path,"lqa_1.0-3.tar.gz",sep=""), repos = NULL, type="source")
library(lqa)


#Set the number of simulation
itr_max = 30
# sample size
n = 1000
# total number of predictors
p = 20

######### define scenarios #########
# scenario = 0 # can be 0 ~ 4
for (scenario in c(1)){

# set information for simulating coviariates
mean_x = 0
sig_x = 1
# rho_list = c(0, 0.25, 0.5, 0.75)
rho_list = c(0, 0.25, 0.5, 0.75)





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
method_list = c("Targ", 
                "Enh_lr_tanh", "Enh_lr_tanh_2","Enh_lr_tanh_3","Enh_lr_tanh_4",
                "Enh_lr_sigmoid", "Enh_lr_sigmoid_2","Enh_lr_sigmoid_3","Enh_lr_sigmoid_4",
                "Enh_svm_tanh", "Enh_svm_tanh_2","Enh_svm_tanh_3","Enh_svm_tanh_4",
                "Enh_svm_sigmoid","Enh_svm_sigmoid_2","Enh_svm_sigmoid_3","Enh_svm_sigmoid_4"
                )

#declare a zero matrix to hold att
att.target = zeros(itr_max)

att.Enh_lr_tanh = zeros(itr_max)
att.Enh_lr_tanh_2 = zeros(itr_max)
att.Enh_lr_tanh_3 = zeros(itr_max)
att.Enh_lr_tanh_4 = zeros(itr_max)
att.Enh_lr_sigmoid = zeros(itr_max)
att.Enh_lr_sigmoid_2 = zeros(itr_max)
att.Enh_lr_sigmoid_3 = zeros(itr_max)
att.Enh_lr_sigmoid_4 = zeros(itr_max)
att.Enh_svm_tanh = zeros(itr_max) 
att.Enh_svm_tanh_2 = zeros(itr_max) 
att.Enh_svm_tanh_3 = zeros(itr_max) 
att.Enh_svm_tanh_4 = zeros(itr_max)
att.Enh_svm_sigmoid = zeros(itr_max) 
att.Enh_svm_sigmoid_2 = zeros(itr_max) 
att.Enh_svm_sigmoid_3 = zeros(itr_max) 
att.Enh_svm_sigmoid_4 = zeros(itr_max) 


#declare zero matrices to store the indexes of selected variables
select.var.list.Enh_lr_tanh = zeros(itr_max, p)
select.var.list.Enh_lr_tanh_2 = zeros(itr_max, p)
select.var.list.Enh_lr_tanh_3 = zeros(itr_max, p)
select.var.list.Enh_lr_tanh_4 = zeros(itr_max, p)
select.var.list.Enh_lr_sigmoid = zeros(itr_max, p)
select.var.list.Enh_lr_sigmoid_2 = zeros(itr_max, p)
select.var.list.Enh_lr_sigmoid_3 = zeros(itr_max, p)
select.var.list.Enh_lr_sigmoid_4 = zeros(itr_max, p)
select.var.list.Enh_svm_tanh = zeros(itr_max, p) 
select.var.list.Enh_svm_tanh_2 = zeros(itr_max, p) 
select.var.list.Enh_svm_tanh_3 = zeros(itr_max, p) 
select.var.list.Enh_svm_tanh_4 = zeros(itr_max, p) 
select.var.list.Enh_svm_sigmoid = zeros(itr_max, p) 
select.var.list.Enh_svm_sigmoid_2 = zeros(itr_max, p) 
select.var.list.Enh_svm_sigmoid_3 = zeros(itr_max, p) 
select.var.list.Enh_svm_sigmoid_4 = zeros(itr_max, p) 

#Declare zero matrices to store the total number of variables selected
total.num.var.Enh_lr_tanh = zeros(itr_max)
total.num.var.Enh_lr_tanh_2 = zeros(itr_max)
total.num.var.Enh_lr_tanh_3 = zeros(itr_max)
total.num.var.Enh_lr_tanh_4 = zeros(itr_max)
total.num.var.Enh_lr_sigmoid = zeros(itr_max)
total.num.var.Enh_lr_sigmoid_2 = zeros(itr_max)
total.num.var.Enh_lr_sigmoid_3 = zeros(itr_max)
total.num.var.Enh_lr_sigmoid_4 = zeros(itr_max)
total.num.var.Enh_svm_tanh = zeros(itr_max)
total.num.var.Enh_svm_tanh_2 = zeros(itr_max)
total.num.var.Enh_svm_tanh_3 = zeros(itr_max)
total.num.var.Enh_svm_tanh_4 = zeros(itr_max)
total.num.var.Enh_svm_sigmoid = zeros(itr_max) 
total.num.var.Enh_svm_sigmoid_2 = zeros(itr_max) 
total.num.var.Enh_svm_sigmoid_3 = zeros(itr_max) 
total.num.var.Enh_svm_sigmoid_4 = zeros(itr_max)

#Declare total percentage of selected variables for each method
total_sel_var = zeros(length(method_list), p)

#Declare total failure times
total_failure_times = zeros(length(method_list))

#Comparison for svm and lr, record total negative, confounder negative, outcome negative
Comp_svm_lr_Percent = zeros(length(rho_list), 3)

#Time Elaspsed
t_elapse = zeros(length(method_list)-3, length(rho_list)) # exclude "tar","conf"&"pot.conf"

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

  #Declare percentage of the variables selected
  Var_Select = zeros(length(method_list), p)

  # Failure times of when apply matching
  fail_times = zeros(length(method_list))

  non_zero_T = 0

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

    xx=Data[,1:p]
    xa=Data[,1:(p+1)]
    yy=Data[,(p+2)]
    aa = Data[,(p+1)]

    #print(p)
    #### estimate outcome model ####
    t0s = Sys.time()

    y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
    lm.Y = lm(y.form,data=Data) #(with OLS)
    betaXY = coef(lm.Y)[var.list]

    gt = Sys.time()-t0s



    ######### 1. Target model #########
    label = label+1

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
    for (v in sel.var.index){
      Var_Select[label, v] <- Var_Select[label, v]+1
    }

    print("target")
    print(att.target[itr])

    ########### Enh_lr_tanh (0.5 & 1) (0.5 & 1) ###################
    label = label+1

    # time recording
    ts = Sys.time()
    
    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)

    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept

    factor = 0.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000

    penalty <- vec_tanh(penalty)^factor
    # penalty <- append(0, penalty)

    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)

    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])

    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    penalty_adel[is.na(penalty_adel)] <- 10000000000

    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)

    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item

    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    xx_true = xx[, beta_aden_non_zero]
    Enh_lr_tanh.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_lr_tanh.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_lr_tanh[itr] = NA
      select.var.list.Enh_lr_tanh[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_tanh[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))

      if ('try-error' %in% class(mm)){
        att.Enh_lr_tanh[itr] = NA
        select.var.list.Enh_lr_tanh[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_tanh[itr] = NA
        fail_times[label] = fail_times[label]+1
      } else{
        effect.form = formula(paste("Y~A+",paste(Enh_lr_tanh.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))

        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_lr_tanh[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_lr_tanh[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_lr_tanh[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }

        print("select.var.list.Enh_lr_tanh")
        print(sel.var.index)
        print("att.Enh_lr_tanh")
        print(att.Enh_lr_tanh[itr])
      }
    }

    
    ########### Enh_lr_tanh_2 (0.5 & 2) (0.8 & 1) ###################
    label = label+1
    
    # time recording
    ts = Sys.time()
    
    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 0.5
    factor = 0.8
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    # penalty <- append(0, penalty)
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 2
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    Enh_lr_tanh_2.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_lr_tanh_2.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_lr_tanh_2[itr] = NA
      select.var.list.Enh_lr_tanh_2[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_tanh_2[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_lr_tanh_2[itr] = NA
        select.var.list.Enh_lr_tanh_2[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_tanh_2[itr] = NA
        fail_times[label] = fail_times[label]+1
      } else{
        effect.form = formula(paste("Y~A+",paste(Enh_lr_tanh_2.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))
        
        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_lr_tanh_2[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_lr_tanh_2[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_lr_tanh_2[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }
        
        print("select.var.list.Enh_lr_tanh_2")
        print(sel.var.index)
        print("att.Enh_lr_tanh_2")
        print(att.Enh_lr_tanh_2[itr])
      }
    }
    
    
    ########### Enh_lr_tanh_3 (0.5 & 3) (1 & 1) ###################
    label = label+1
    
    # time recording
    ts = Sys.time()
    
    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 0.5
    factor = 1
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    # penalty <- append(0, penalty)
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 3
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 <- try(aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE))
    if ('try-error' %in% class(fit3)){
      att.Enh_lr_tanh_3[itr] = NA
      select.var.list.Enh_lr_tanh_3[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_tanh_3[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
    
      beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
      beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
      
      beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
      sel.var.index = beta_aden_non_zero
      
      t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)
      
      xx_true = xx[, beta_aden_non_zero]
      Enh_lr_tanh_3.true.var.list = names(xx_true)
      treat.form <- try(formula(paste("A~",paste(Enh_lr_tanh_3.true.var.list,collapse="+"))))
      if ('try-error' %in% class(treat.form)){
        att.Enh_lr_tanh_3[itr] = NA
        select.var.list.Enh_lr_tanh_3[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_tanh_3[itr] = NA
        fail_times[label] = fail_times[label]+1
      }else{
        # 1 to 1 nearest matching
        mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                          ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                          replace = FALSE, discard = "both"))
        
        if ('try-error' %in% class(mm)){
          att.Enh_lr_tanh_3[itr] = NA
          select.var.list.Enh_lr_tanh_3[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_lr_tanh_3[itr] = NA
          fail_times[label] = fail_times[label]+1
        } else{
          effect.form = formula(paste("Y~A+",paste(Enh_lr_tanh_3.true.var.list,collapse="+")))
          fit.effect = lm(effect.form, data = match.data(mm))
          
          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.Enh_lr_tanh_3[itr] = final.coef[2]
          # Storing the Indexes of selected variables
          select.var.list.Enh_lr_tanh_3[itr, 1:length(sel.var.index)] = sel.var.index
          # Storing the total number of variables selected
          total.num.var.Enh_lr_tanh_3[itr] = length(sel.var.index)
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }
          
          print("select.var.list.Enh_lr_tanh_3")
          print(sel.var.index)
          print("att.Enh_lr_tanh_3")
          print(att.Enh_lr_tanh_3[itr])
        }
      }
    }
    
    ########### Enh_lr_tanh_4 (0.5 & 4) (1.5 & 1) ###################
    label = label+1
    
    # time recording
    ts = Sys.time()
    
    fit1 = glm(A~.,family=binomial(link='logit'),data=xa)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 0.5
    factor = 1.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    # penalty <- append(0, penalty)
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 4
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    att.Enh_lr_tanh_4[itr] = NA
    select.var.list.Enh_lr_tanh_4[itr, 1:length(sel.var.index)] = 0
    total.num.var.Enh_lr_tanh_4[itr] = NA
    fail_times[label] = fail_times[label]+1
  
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    Enh_lr_tanh_4.true.var.list = names(xx_true)
    
    if (is.null(Enh_lr_tanh_4.true.var.list)){
      att.Enh_lr_tanh_4[itr] = NA
      select.var.list.Enh_lr_tanh_4[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_tanh_4[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      treat.form <- try(formula(paste("A~",paste(Enh_lr_tanh_4.true.var.list,collapse="+"))))
      if ('try-error' %in% class(treat.form)){
        att.Enh_lr_tanh_4[itr] = NA
        select.var.list.Enh_lr_tanh_4[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_tanh_4[itr] = NA
        fail_times[label] = fail_times[label]+1
      }else{
        # 1 to 1 nearest matching
        mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                          ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                          replace = FALSE, discard = "both"))
        
        if ('try-error' %in% class(mm)){
          att.Enh_lr_tanh_4[itr] = NA
          select.var.list.Enh_lr_tanh_4[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_lr_tanh_4[itr] = NA
          fail_times[label] = fail_times[label]+1
        }else{
          if (is.null(Enh_lr_tanh_4.true.var.list)){
            att.Enh_lr_tanh_4[itr] = NA
            select.var.list.Enh_lr_tanh_4[itr, 1:length(sel.var.index)] = 0
            total.num.var.Enh_lr_tanh_4[itr] = NA
            fail_times[label] = fail_times[label]+1
          }else{
            effect.form <- try(formula(paste("Y~A+",paste(Enh_lr_tanh_4.true.var.list,collapse="+"))))
            if ('try-error' %in% class(effect.form)){
              att.Enh_lr_tanh_4[itr] = NA
              select.var.list.Enh_lr_tanh_4[itr, 1:length(sel.var.index)] = 0
              total.num.var.Enh_lr_tanh_4[itr] = NA
              fail_times[label] = fail_times[label]+1
            }else{
              fit.effect = lm(effect.form, data = match.data(mm))
              
              ### post record
              final.coef = fit.effect$coefficients
              # att
              att.Enh_lr_tanh_4[itr] = final.coef[2]
              # Storing the Indexes of selected variables
              select.var.list.Enh_lr_tanh_4[itr, 1:length(sel.var.index)] = sel.var.index
              # Storing the total number of variables selected
              total.num.var.Enh_lr_tanh_4[itr] = length(sel.var.index)
              # Percent var select
              for (v in sel.var.index){
                Var_Select[label, v] <- Var_Select[label, v]+1
              }
              
              print("select.var.list.Enh_lr_tanh_4")
              print(sel.var.index)
              print("att.Enh_lr_tanh_4")
              print(att.Enh_lr_tanh_4[itr])
            }
          }
        }
      }
    }
    
    
    ########### Enh_lr_sigmoid (1 & 1) (0.5 & 1) ###################
    label = label+1
    # time recording
    ts = Sys.time()

    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)

    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept

    # factor = 1
    factor = 0.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000

    penalty <- vec_sigmoid(penalty)^factor

    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)

    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])

    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma

    penalty_adel[is.na(penalty_adel)] <- 10000000000

    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)

    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item

    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    xx_true = xx[, beta_aden_non_zero]
    Enh_lr_sigmoid.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_lr_sigmoid.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_lr_sigmoid[itr] = NA
      select.var.list.Enh_lr_sigmoid[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_sigmoid[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))

      if ('try-error' %in% class(mm)){
        att.Enh_lr_sigmoid[itr] = NA
        select.var.list.Enh_lr_sigmoid[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_sigmoid[itr] = NA
        fail_times[label] = fail_times[label]+1

      } else{
        effect.form <- try(formula(paste("Y~A+",paste(Enh_lr_tanh_4.true.var.list,collapse="+"))))
        if ('try-error' %in% class(effect.form)){
          att.Enh_lr_sigmoid[itr] = NA
          select.var.list.Enh_lr_sigmoid[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_lr_sigmoid[itr] = NA
          fail_times[label] = fail_times[label]+1
        }else{
          fit.effect = lm(effect.form, data = match.data(mm))
  
          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.Enh_lr_sigmoid[itr] = final.coef[2]
          # Storing the Indexes of selected variables
          select.var.list.Enh_lr_sigmoid[itr, 1:length(sel.var.index)] = sel.var.index
          # Storing the total number of variables selected
          total.num.var.Enh_lr_sigmoid[itr] = length(sel.var.index)
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }
  
          print("select.var.list.Enh_lr_sigmoid")
          print(sel.var.index)
          print("att.Enh_lr_sigmoid")
          print(att.Enh_lr_sigmoid[itr])
        }
      }
    }

    
    
    ########### Enh_lr_sigmoid_2 (1 & 2) (1 & 1) ###################
    label = label+1
    # time recording
    ts = Sys.time()
    
    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    factor = 1
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_sigmoid(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 2
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    Enh_lr_sigmoid_2.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_lr_sigmoid_2.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_lr_sigmoid_2[itr] = NA
      select.var.list.Enh_lr_sigmoid_2[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_sigmoid_2[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_lr_sigmoid_2[itr] = NA
        select.var.list.Enh_lr_sigmoid_2[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_sigmoid_2[itr] = NA
        fail_times[label] = fail_times[label]+1
        
      } else{
        effect.form <- try(formula(paste("Y~A+",paste(Enh_lr_tanh_4.true.var.list,collapse="+"))))
        if ('try-error' %in% class(effect.form)){
          att.Enh_lr_sigmoid_2[itr] = NA
          select.var.list.Enh_lr_sigmoid_2[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_lr_sigmoid_2[itr] = NA
          fail_times[label] = fail_times[label]+1
        }else{
          fit.effect = lm(effect.form, data = match.data(mm))
          
          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.Enh_lr_sigmoid_2[itr] = final.coef[2]
          # Storing the Indexes of selected variables
          select.var.list.Enh_lr_sigmoid_2[itr, 1:length(sel.var.index)] = sel.var.index
          # Storing the total number of variables selected
          total.num.var.Enh_lr_sigmoid_2[itr] = length(sel.var.index)
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }
          
          print("select.var.list.Enh_lr_sigmoid_2")
          print(sel.var.index)
          print("att.Enh_lr_sigmoid_2")
          print(att.Enh_lr_sigmoid_2[itr])
        }
      }
    }
    
    
    ########### Enh_lr_sigmoid_3 (1 & 3) (1.5 & 1) ###################
    label = label+1
    # time recording
    ts = Sys.time()
    
    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 1
    factor = 1.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_sigmoid(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 3
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    Enh_lr_sigmoid_3.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_lr_sigmoid_3.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_lr_sigmoid_3[itr] = NA
      select.var.list.Enh_lr_sigmoid_3[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_sigmoid_3[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_lr_sigmoid_3[itr] = NA
        select.var.list.Enh_lr_sigmoid_3[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_sigmoid_3[itr] = NA
        fail_times[label] = fail_times[label]+1
        
      } else{
        effect.form <- try(formula(paste("Y~A+",paste(Enh_lr_tanh_4.true.var.list,collapse="+"))))
        if ('try-error' %in% class(effect.form)){
          att.Enh_lr_sigmoid_3[itr] = NA
          select.var.list.Enh_lr_sigmoid_3[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_lr_sigmoid_3[itr] = NA
          fail_times[label] = fail_times[label]+1
        }else{
          fit.effect = lm(effect.form, data = match.data(mm))
          
          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.Enh_lr_sigmoid_3[itr] = final.coef[2]
          # Storing the Indexes of selected variables
          select.var.list.Enh_lr_sigmoid_3[itr, 1:length(sel.var.index)] = sel.var.index
          # Storing the total number of variables selected
          total.num.var.Enh_lr_sigmoid_3[itr] = length(sel.var.index)
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }
          
          print("select.var.list.Enh_lr_sigmoid_3")
          print(sel.var.index)
          print("att.Enh_lr_sigmoid_3")
          print(att.Enh_lr_sigmoid_3[itr])
        }
      }
    }
    
    ########### Enh_lr_sigmoid_4 (1 & 4) (2 & 1) ###################
    label = label+1
    # time recording
    ts = Sys.time()
    
    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 1
    factor = 2
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_sigmoid(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 4
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    Enh_lr_sigmoid_4.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_lr_sigmoid_4.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_lr_sigmoid_4[itr] = NA
      select.var.list.Enh_lr_sigmoid_4[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_lr_sigmoid_4[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_lr_sigmoid_4[itr] = NA
        select.var.list.Enh_lr_sigmoid_4[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_lr_sigmoid_4[itr] = NA
        fail_times[label] = fail_times[label]+1
        
      } else{
        effect.form <- try(formula(paste("Y~A+",paste(Enh_lr_tanh_4.true.var.list,collapse="+"))))
        if ('try-error' %in% class(effect.form)){
          att.Enh_lr_sigmoid_4[itr] = NA
          select.var.list.Enh_lr_sigmoid_4[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_lr_sigmoid_4[itr] = NA
          fail_times[label] = fail_times[label]+1
        }else{
          fit.effect = lm(effect.form, data = match.data(mm))
          
          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.Enh_lr_sigmoid_4[itr] = final.coef[2]
          # Storing the Indexes of selected variables
          select.var.list.Enh_lr_sigmoid_4[itr, 1:length(sel.var.index)] = sel.var.index
          # Storing the total number of variables selected
          total.num.var.Enh_lr_sigmoid_4[itr] = length(sel.var.index)
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }
          
          print("select.var.list.Enh_lr_sigmoid_4")
          print(sel.var.index)
          print("att.Enh_lr_sigmoid_4")
          print(att.Enh_lr_sigmoid_4[itr])
        }
      }
    }
    
    
    
    ########### Enh_svm_tanh (0.5 & 1) (0.5 & 1) ######
    label = label+1
    
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    factor = 0.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    Enh_svm_tanh.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_svm_tanh.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_svm_tanh[itr] = NA
      select.var.list.Enh_svm_tanh[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_tanh[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_svm_sigmoid[itr] = NA
        select.var.list.Enh_svm_tanh[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_tanh[itr] = NA
        fail_times[label] = fail_times[label]+1
      } else{
        effect.form = formula(paste("Y~A+",paste(Enh_svm_tanh.true.var.list, collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))
        
        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_svm_tanh[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_svm_tanh[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_svm_tanh[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }
        
        print("select.var.list.Enh_svm_tanh")
        print(sel.var.index)
        print("att.Enh_svm_tanh")
        print(att.Enh_svm_tanh[itr])
      }
    }
    
    
    
    ########### Enh_svm_tanh (0.5 & 2) (0.8 & 1) ######
    label = label+1
    
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 0.5
    factor = 0.8
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 2
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    Enh_svm_tanh_2.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(Enh_svm_tanh_2.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_svm_tanh_2[itr] = NA
      select.var.list.Enh_svm_tanh_2[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_tanh_2[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_svm_sigmoid[itr] = NA
        select.var.list.Enh_svm_tanh_2[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_tanh_2[itr] = NA
        fail_times[label] = fail_times[label]+1
      } else{
        effect.form = formula(paste("Y~A+",paste(Enh_svm_tanh_2.true.var.list, collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))
        
        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_svm_tanh_2[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_svm_tanh_2[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_svm_tanh_2[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }
        
        print("select.var.list.Enh_svm_tanh_2")
        print(sel.var.index)
        print("att.Enh_svm_tanh_2")
        print(att.Enh_svm_tanh_2[itr])
      }
    }
    
    
    ########### Enh_svm_tanh_3 (0.5 & 3) (1.0 & 1)######
    label = label+1
    
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 0.5
    factor = 1
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 3
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 <- try(aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE))
    if ('try-error' %in% class(fit3)){
      att.Enh_svm_tanh_3[itr] = NA
      select.var.list.Enh_svm_tanh_3[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_tanh_3[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
      beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
      
      beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
      sel.var.index = beta_aden_non_zero
      
      t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
      
      xx_true = xx[, beta_aden_non_zero]
      Enh_svm_tanh_3.true.var.list = names(xx_true)
      treat.form <- try(formula(paste("A~",paste(Enh_svm_tanh_3.true.var.list,collapse="+"))))
      if ('try-error' %in% class(treat.form)){
        att.Enh_svm_tanh_3[itr] = NA
        select.var.list.Enh_svm_tanh_3[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_tanh_3[itr] = NA
        fail_times[label] = fail_times[label]+1
      }else{
        # 1 to 1 nearest matching
        mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                          ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                          replace = FALSE, discard = "both"))
        
        if ('try-error' %in% class(mm)){
          att.Enh_svm_sigmoid[itr] = NA
          select.var.list.Enh_svm_tanh_3[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_svm_tanh_3[itr] = NA
          fail_times[label] = fail_times[label]+1
        } else{
          effect.form = formula(paste("Y~A+",paste(Enh_svm_tanh_3.true.var.list, collapse="+")))
          fit.effect = lm(effect.form, data = match.data(mm))
          
          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.Enh_svm_tanh_3[itr] = final.coef[2]
          # Storing the Indexes of selected variables
          select.var.list.Enh_svm_tanh_3[itr, 1:length(sel.var.index)] = sel.var.index
          # Storing the total number of variables selected
          total.num.var.Enh_svm_tanh_3[itr] = length(sel.var.index)
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }
          
          print("select.var.list.Enh_svm_tanh_3")
          print(sel.var.index)
          print("att.Enh_svm_tanh_3")
          print(att.Enh_svm_tanh_3[itr])
        }
      }
    }
    
    ########### Enh_svm_tanh_4 (0.5 & 4) (1.5 & 1)######
    label = label+1
    
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 0.5
    factor = 1.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_tanh(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 4
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 <- try(aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE))
    if ('try-error' %in% class(fit3)){
      att.Enh_svm_tanh_4[itr] = NA
      select.var.list.Enh_svm_tanh_4[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_tanh_4[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
    
      beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
      beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
      
      beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
      sel.var.index = beta_aden_non_zero
      
      t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
      
      xx_true = xx[, beta_aden_non_zero]
      Enh_svm_tanh_4.true.var.list = names(xx_true)
      treat.form <- try(formula(paste("A~",paste(Enh_svm_tanh_4.true.var.list,collapse="+"))))
      if ('try-error' %in% class(treat.form)){
        att.Enh_svm_tanh_4[itr] = NA
        select.var.list.Enh_svm_tanh_4[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_tanh_4[itr] = NA
        fail_times[label] = fail_times[label]+1
      }else{
        # 1 to 1 nearest matching
        mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                          ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                          replace = FALSE, discard = "both"))
        
        if ('try-error' %in% class(mm)){
          att.Enh_svm_sigmoid[itr] = NA
          select.var.list.Enh_svm_tanh_4[itr, 1:length(sel.var.index)] = 0
          total.num.var.Enh_svm_tanh_4[itr] = NA
          fail_times[label] = fail_times[label]+1
        } else{
          effect.form = formula(paste("Y~A+",paste(Enh_svm_tanh_4.true.var.list, collapse="+")))
          fit.effect = lm(effect.form, data = match.data(mm))
          
          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.Enh_svm_tanh_4[itr] = final.coef[2]
          # Storing the Indexes of selected variables
          select.var.list.Enh_svm_tanh_4[itr, 1:length(sel.var.index)] = sel.var.index
          # Storing the total number of variables selected
          total.num.var.Enh_svm_tanh_4[itr] = length(sel.var.index)
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }
          
          print("select.var.list.Enh_svm_tanh_4")
          print(sel.var.index)
          print("att.Enh_svm_tanh_4")
          print(att.Enh_svm_tanh_4[itr])
        }
      }
    }
    
    
    
    ########### Enh_svm_sigmoid (1 & 1) (0.5 & 1) ###################
    label = label+1
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 1
    factor = 0.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_sigmoid(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    svm_sigmoid.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(svm_sigmoid.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_svm_sigmoid[itr] = NA
      select.var.list.Enh_svm_sigmoid[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_sigmoid[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_svm_sigmoid[itr] = NA
        select.var.list.Enh_svm_sigmoid[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_sigmoid[itr] = NA
        fail_times[label] = fail_times[label]+1
        
      }else{
        effect.form = formula(paste("Y~A+",paste(svm_sigmoid.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))
        
        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_svm_sigmoid[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_svm_sigmoid[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_svm_sigmoid[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }
        
        print("select.var.list.svm_sigmoid")
        print(sel.var.index)
        print("att.svm_sigmoid")
        print(att.Enh_svm_sigmoid[itr])
      }
    }
    
    
    ########### Enh_svm_sigmoid_2 (1 & 2) (1 & 1)###################
    label = label+1
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    factor = 1
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_sigmoid(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 2
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    svm_sigmoid.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(svm_sigmoid.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_svm_sigmoid_2[itr] = NA
      select.var.list.Enh_svm_sigmoid_2[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_sigmoid_2[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_svm_sigmoid_2[itr] = NA
        select.var.list.Enh_svm_sigmoid_2[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_sigmoid_2[itr] = NA
        fail_times[label] = fail_times[label]+1
        
      }else{
        effect.form = formula(paste("Y~A+",paste(svm_sigmoid.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))
        
        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_svm_sigmoid_2[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_svm_sigmoid_2[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_svm_sigmoid_2[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }
        
        print("select.var.list.svm_sigmoid")
        print(sel.var.index)
        print("att.svm_sigmoid")
        print(att.Enh_svm_sigmoid_2[itr])
      }
    }
    
    ########### Enh_svm_sigmoid_3 (1 & 3) (1.5 & 1) ###################
    label = label+1
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 1
    factor = 1.5
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_sigmoid(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 3
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    svm_sigmoid.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(svm_sigmoid.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_svm_sigmoid_3[itr] = NA
      select.var.list.Enh_svm_sigmoid_3[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_sigmoid_3[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_svm_sigmoid_3[itr] = NA
        select.var.list.Enh_svm_sigmoid_3[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_sigmoid_3[itr] = NA
        fail_times[label] = fail_times[label]+1
        
      }else{
        effect.form = formula(paste("Y~A+",paste(svm_sigmoid.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))
        
        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_svm_sigmoid_3[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_svm_sigmoid_3[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_svm_sigmoid_3[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }
        
        print("select.var.list.svm_sigmoid")
        print(sel.var.index)
        print("att.svm_sigmoid")
        print(att.Enh_svm_sigmoid_3[itr])
      }
    }
    
    ########### Enh_svm_sigmoid_4 (1 & 4) (2 & 1)###################
    label = label+1
    # time recording
    ts = Sys.time()
    
    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    
    # factor = 1
    factor = 2
    penalty = abs(weight_logit)
    penalty[is.na(penalty)] <- 10000000000
    
    penalty <- vec_sigmoid(penalty)^factor
    
    # elastic net estimator
    fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                     nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                     relax = TRUE, parallel=TRUE, standardize=TRUE,
                     type.measure='mse', penalty.factor=penalty)
    # fit2 = glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5, penalty.factor = penalty)
    
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    
    # gamma = 4
    gamma = 1
    penalty_adel = (1/abs(beta_el_allvar))^gamma
    
    penalty_adel[is.na(penalty_adel)] <- 10000000000
    
    # adaptive elastic net estimator
    fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty_adel, parallel = TRUE)
    
    beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    
    beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
    sel.var.index = beta_aden_non_zero
    
    t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)
    
    xx_true = xx[, beta_aden_non_zero]
    svm_sigmoid.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(svm_sigmoid.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.Enh_svm_sigmoid_4[itr] = NA
      select.var.list.Enh_svm_sigmoid_4[itr, 1:length(sel.var.index)] = 0
      total.num.var.Enh_svm_sigmoid_4[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      # 1 to 1 nearest matching
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))
      
      if ('try-error' %in% class(mm)){
        att.Enh_svm_sigmoid_4[itr] = NA
        select.var.list.Enh_svm_sigmoid_4[itr, 1:length(sel.var.index)] = 0
        total.num.var.Enh_svm_sigmoid_4[itr] = NA
        fail_times[label] = fail_times[label]+1
        
      }else{
        effect.form = formula(paste("Y~A+",paste(svm_sigmoid.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))
        
        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.Enh_svm_sigmoid_4[itr] = final.coef[2]
        # Storing the Indexes of selected variables
        select.var.list.Enh_svm_sigmoid_4[itr, 1:length(sel.var.index)] = sel.var.index
        # Storing the total number of variables selected
        total.num.var.Enh_svm_sigmoid_4[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }
        
        print("select.var.list.svm_sigmoid")
        print(sel.var.index)
        print("att.svm_sigmoid")
        print(att.Enh_svm_sigmoid_4[itr])
      }
    }
    
    
  }

  title = paste("s_",scenario,"_att_",bA,"_rho_",rho,"_n_",n,"_p_",p,"_itr_",itr,"_comp_thrs", sep="")

  # save boxplot
  # png(filename=paste(path, "figure/", title,".png",sep=""))
  # boxplot(att.target, att.conf, att.potconf,
  #         att.Enh_lr_tanh, att.Enh_lr_sigmoid, att.svm_sigmoid,
  #         att.lr_sigmoid, att.lr_bare,
  #         att.adl,att.onet,att.BACR,att.BCEE,
  #         names = method_list)
  # abline(h=bA, lty = 2)
  # dev.off()


  # save att
  df_att = data.frame(att.target,
                      att.Enh_lr_tanh, att.Enh_lr_tanh_2, 
                      att.Enh_lr_tanh_3, att.Enh_lr_tanh_4,
                      att.Enh_lr_sigmoid, att.Enh_lr_sigmoid_2,
                      att.Enh_lr_sigmoid_3, att.Enh_lr_sigmoid_4,
                      att.Enh_svm_tanh, att.Enh_svm_tanh_2,
                      att.Enh_svm_tanh_3, att.Enh_svm_tanh_4,
                      att.Enh_svm_sigmoid, att.Enh_svm_sigmoid_2,
                      att.Enh_svm_sigmoid_3, att.Enh_svm_sigmoid_4
                      )

  # save num of vars selected
  df_nsel.var = data.frame(total.num.var.Enh_lr_tanh,
                           total.num.var.Enh_lr_tanh_2,
                           total.num.var.Enh_lr_tanh_3,
                           total.num.var.Enh_lr_tanh_4,
                           total.num.var.Enh_lr_sigmoid,
                           total.num.var.Enh_lr_sigmoid_2,
                           total.num.var.Enh_lr_sigmoid_3,
                           total.num.var.Enh_lr_sigmoid_4,
                           total.num.var.Enh_svm_tanh, 
                           total.num.var.Enh_svm_tanh_2, 
                           total.num.var.Enh_svm_tanh_3, 
                           total.num.var.Enh_svm_tanh_4, 
                           total.num.var.Enh_svm_sigmoid,
                           total.num.var.Enh_svm_sigmoid_2,
                           total.num.var.Enh_svm_sigmoid_3,
                           total.num.var.Enh_svm_sigmoid_4)

  # save times of vars selected
  df_selt.var = data.frame(Var_Select)
  total_sel_var = total_sel_var + df_selt.var
  row.names(df_selt.var) <- method_list

  # save time recording
  success_time = rep(itr_max, length(method_list))
  t_elapse[,k] = t_elapse[,k]/(success_time-fail_times)[4:length(method_list)]
  df_time = data.frame(t_elapse[,k])
  row.names(df_time) <- method_list[4:length(method_list)]

  # writing how many times fail
  df_fail = data.frame(fail_times)
  total_failure_times = total_failure_times+fail_times
  row.names(df_fail) <- method_list

  # svm and lr comparison but do not write until all rhos are done
  Comp_svm_lr_Percent[k, 1] = Comp_svm_lr_Percent[k, 1]/(pI*itr_max)# total negative
  Comp_svm_lr_Percent[k, 2] = Comp_svm_lr_Percent[k, 2]/(pC*itr_max) # confounder chosen negative
  Comp_svm_lr_Percent[k, 3] = Comp_svm_lr_Percent[k, 3]/(pP*itr_max) # outcome predictor negative


  write.csv(df_att, paste(path, "att/", title,"att.csv",sep=""))
  write.csv(df_nsel.var, paste(path, "nvar/",title,"nvar.csv",sep=""))
  write.csv(df_selt.var, paste(path, "tvar/",title,"vars_test.csv",sep=""))
  write.csv(df_time, paste(path, "time/",title,"time.csv",sep=""))
  write.csv(df_fail, paste(path, "fail/",title,"fail.csv",sep=""))
}

# # svm and lr comparison but do not write until all rhos are done
# df_comp = data.frame(Comp_svm_lr_Percent)
# row.names(df_comp) <- rho_list
# colnames(df_comp) <- c("%p negative", "%pC negative", "%pP negative")

df_total_sel_var = data.frame(total_sel_var)
row.names(df_total_sel_var) <- method_list

df_total_failure = data.frame(total_failure_times)
row.names(df_total_failure) <- method_list

df_total_time = data.frame(t_elapse/itr_max)
colnames(df_total_time) <- rho_list
row.names(df_total_time) <- method_list[4:length(method_list)]

# write.csv(df_comp, paste(path, "comp/s",scenario, "_comp_",n,"_",p,".csv",sep=""))
write.csv(df_total_sel_var, paste(path, "tvar/s",scenario, "_varsel_",n,"_",p,".csv",sep=""))
write.csv(df_total_failure, paste(path, "fail/s",scenario, "_failure_",n,"_",p,".csv",sep=""))
write.csv(df_total_time, paste(path, "time/s",scenario, "_time_",n,"_",p,".csv",sep=""))

}




