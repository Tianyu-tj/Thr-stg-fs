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
path = "~/RProjects/FSProject_Final/"

install.packages(paste(path,"lqa_1.0-3.tar.gz",sep=""), repos = NULL, type="source")
library(lqa)


#Set the number of simulation
itr_max = 30
# sample size
n = 200
# total number of predictors
p = 10

######### define scenarios #########
# scenario = 0 # can be 0 ~ 4
for (scenario in c(0:4)){

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
method_list = c("Targ", "Conf", "PotConf",
                "Enh_lr_tanh", "Enh_lr_sigmoid", "Enh_svm_tanh", "Enh_svm_sigmoid",
                "SVM(sigmoid)", "LogR(sigmoid)", "LogR(bare)",
                "Adl","ONet", "BACR", "BCEE",
                "Boruta_T", "Boruta_Y", "DWR", "CTMLE"
                )

#declare a zero matrix to hold att
att.target = zeros(itr_max)
att.potconf = zeros(itr_max)
att.conf = zeros(itr_max)

att.Enh_lr_tanh = zeros(itr_max)
att.Enh_lr_sigmoid = zeros(itr_max)
att.Enh_svm_tanh = zeros(itr_max) 
att.Enh_svm_sigmoid = zeros(itr_max) 

att.svm_sigmoid = zeros(itr_max)
att.lr_sigmoid = zeros(itr_max)
att.lr_bare = zeros(itr_max)

att.adl = zeros(itr_max)
att.onet = zeros(itr_max)
att.BACR = zeros(itr_max)
att.BCEE = zeros(itr_max)
att.Boruta_T = zeros(itr_max) 
att.Boruta_Y = zeros(itr_max)
att.DWR = zeros(itr_max)
att.CTMLE = zeros(itr_max)



#declare zero matrices to store the indexes of selected variables
select.var.list.Enh_lr_tanh = zeros(itr_max, p)
select.var.list.Enh_lr_sigmoid = zeros(itr_max, p)
select.var.list.Enh_svm_tanh = zeros(itr_max, p) 
select.var.list.Enh_svm_sigmoid = zeros(itr_max, p) 

select.var.list.svm_sigmoid = zeros(itr_max, p)
select.var.list.lr_sigmoid = zeros(itr_max, p)
select.var.list.lr_bare = zeros(itr_max, p)

select.var.list.adl = zeros(itr_max, p)
select.var.list.onet = zeros(itr_max, p)
select.var.list.bacr = zeros(itr_max, p)
select.var.list.bcee = zeros(itr_max, p)
select.var.list.Boruta_T = zeros(itr_max, p) 
select.var.list.Boruta_Y = zeros(itr_max, p)
select.var.list.DWR = zeros(itr_max, p)
select.var.list.CTMLE = zeros(itr_max, p)

#Declare zero matrices to store the total number of variables selected
total.num.var.Enh_lr_tanh = zeros(itr_max)
total.num.var.Enh_lr_sigmoid = zeros(itr_max)
total.num.var.Enh_svm_tanh = zeros(itr_max) 
total.num.var.Enh_svm_sigmoid = zeros(itr_max) 

total.num.var.svm_sigmoid = zeros(itr_max)
total.num.var.lr_sigmoid = zeros(itr_max)
total.num.var.lr_bare = zeros(itr_max)

total.num.var.adl = zeros(itr_max)
total.num.var.onet = zeros(itr_max)
total.num.var.bacr = zeros(itr_max)
total.num.var.bcee = zeros(itr_max)
total.num.var.Boruta_T = zeros(itr_max) 
total.num.var.Boruta_Y = zeros(itr_max)
total.num.var.DWR = zeros(itr_max)
total.num.var.CTMLE = zeros(itr_max)

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

###### WAMD function and relevant settings ###########
lambda_vec = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
names(lambda_vec) = as.character(lambda_vec)
gamma_convergence_factor = 2 # n^(lam)*n^(gamma/2-1)=n^(conv) -> lam=conv+1-gamma/2 -> gamma=2*(conv-lambda+1)
gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
names(gamma_vals) = names(lambda_vec)

wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta))
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){
    this.var = paste("w",varlist[jj],sep="")
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt]
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt])
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt])
    diff_vec[jj] = abs( trt[jj] - untrt[jj] )
  }
  wdiff_vec = diff_vec * abs(beta)
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret)
}

create_weights = function(fp,fA,fw){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
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



    ####### C1. Only adpative elastic net in outcome model########
    # fit1 = glm(yy ~.,data=xx)
    fit1 <- cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                                   nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                                   relax = TRUE, parallel=TRUE, standardize=TRUE,
                                   type.measure='mse')
    weight_logit = coef(fit1)
    weight_logit = weight_logit[-c(1)]
    penalty = 1/abs(weight_logit)

    fit2 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 3L, rule = "lambda.1se",
                 penalty.factor = penalty, parallel = TRUE)
    beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_ps_allvar = as.matrix(beta_el[1:p])
    beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.01)]
    sel.var.index = beta_ps_non_zero
    print(sel.var.index)
    ######################################################

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


# 
    ######### 2. Potential confounder model #######
    label = label+1
    potconf.form = formula(paste("A~",paste(potconf.var.list,collapse="+")))
    mm <- try(matchit(potconf.form, data = Data, method = "nearest", distance = "glm",
                 ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                 replace = FALSE, discard = "both"))
    if ('try-error' %in% class(mm)){
      att.potconf[itr] = NA
    } else{
      effect.form = formula(paste("Y~A+",paste(potconf.var.list,collapse="+")))
      fit.effect = lm(effect.form, data = match.data(mm))
      #coeftest(fit.effect, vcov. = vcovCL)["A",,drop=FALSE]
      final.coef = fit.effect$coefficients
      att.potconf[itr] = final.coef[2]
    }

    sel.var.index = c(1:(pC+pP+pI))
    # Percent var select
    for (v in sel.var.index){
      Var_Select[label, v] <- Var_Select[label, v]+1
    }






    ######### 3. Confounder model #######
    label = label+1
    conf.form = formula(paste("A~",paste(conf.var.list,collapse="+")))
    mm = matchit(conf.form, data = Data, method = "nearest", distance = "glm",
                 ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                 replace = FALSE, discard = "both")
    effect.form = formula(paste("Y~A+",paste(conf.var.list,collapse="+")))
    fit.effect = lm(effect.form, data = match.data(mm))
    #coeftest(fit.effect, vcov. = vcovCL)["A",,drop=FALSE]
    final.coef = fit.effect$coefficients
    att.conf[itr] = final.coef[2]
    sel.var.index = c(1:pC)
    # Percent var select
    for (v in sel.var.index){
      Var_Select[label, v] <- Var_Select[label, v]+1
    }







    ########### 4. Enh_lr_tanh###################
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

    ########### 5. Enh_lr_sigmoid ###################
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
        effect.form = formula(paste("Y~A+",paste(Enh_lr_tanh.true.var.list,collapse="+")))
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


    ########### 6. Enh_svm_tanh ######
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


    ########### 7. Enh_svm_sigmoid ###################
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

    # ########## comp_svm. Enh_svm_sigmoid ###################
    # fit1 = svm(A ~.,
    #            data = xa,
    #            type = 'C-classification',
    #            kernel = 'linear', scale = TRUE)
    # # fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)
    # 
    # weight_logit = coef(fit1)
    # weight_logit = weight_logit[-c(1)] # exclude the intercept
    # 
    # factor = 1
    # penalty = abs(weight_logit)
    # penalty[is.infinite(penalty)] <- 10000000000
    # 
    # penalty <- vec_sigmoid(penalty)^factor
    # 
    # # elastic net estimator
    # fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
    #                  nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
    #                  relax = TRUE, parallel=TRUE, standardize=TRUE,
    #                  type.measure='mse', penalty.factor=penalty)
    # # fit2 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
    # #              nfolds = 5L, rule = "lambda.1se",
    # #              penalty.factor = penalty, parallel = TRUE)
    # 
    # beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    # beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    # 
    # gamma = 1
    # penalty_adel = (1/abs(beta_el_allvar))^gamma
    # 
    # penalty_adel[is.infinite(penalty_adel)] <- 10000000000
    # svm_pen = penalty_adel
    # 
    # # adaptive elastic net estimator
    # fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
    #              nfolds = 5L, rule = "lambda.1se",
    #              penalty.factor = penalty_adel, parallel = TRUE)
    # 
    # beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    # beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    # 
    # svm_coef = beta_aden_allvar
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # ########## comp_svm. Enh_lr_sigmoid ###################
    # # fit1 = svm(A ~.,
    # #            data = xa,
    # #            type = 'C-classification',
    # #            kernel = 'linear', scale = TRUE)
    # fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)
    # 
    # weight_logit = coef(fit1)
    # weight_logit = weight_logit[-c(1)] # exclude the intercept
    # 
    # factor = 1
    # penalty = abs(weight_logit)
    # penalty[is.infinite(penalty)] <- 10000000000
    # 
    # penalty <- vec_sigmoid(penalty)^factor
    # 
    # # elastic net estimator
    # fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
    #                  nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
    #                  relax = TRUE, parallel=TRUE, standardize=TRUE,
    #                  type.measure='mse', penalty.factor=penalty)
    # # fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
    # #              nfolds = 5L, rule = "lambda.1se",
    # #              penalty.factor = penalty_adel, parallel = TRUE)
    # 
    # beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
    # beta_el_allvar = as.matrix(beta_el[2:(p+1)])
    # 
    # gamma = 1
    # penalty_adel = (1/abs(beta_el_allvar))^gamma
    # 
    # penalty_adel[is.infinite(penalty_adel)] <- 10000000000
    # lr_pen = penalty_adel
    # 
    # # adaptive elastic net estimator
    # fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
    #              nfolds = 5L, rule = "lambda.1se",
    #              penalty.factor = penalty_adel, parallel = TRUE)
    # 
    # beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
    # beta_aden_allvar = as.matrix(beta_aden_raw[1:p]) #remove intercept to the last item
    # 
    # lr_coef = beta_aden_allvar
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 

    ########### 8. svm_sigmoid ##########
    label = label+1
    ts = Sys.time()

    fit1 = svm(A ~.,
               data = xa,
               type = 'C-classification',
               kernel = 'linear', scale = TRUE)
    # fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)

    weight_svm = coef(fit1)
    abs(weight_svm)
    weight_svm = weight_svm[-c(1)] # exclude the intercept
    oddgamma = 1
    penalty = abs(weight_svm) ^ oddgamma
    penalty[is.na(penalty)] <- 10000000000
    penalty = vec_sigmoid(penalty)
    svm_penalty = penalty

    fit2 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = svm_penalty , parallel = TRUE)
    # fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
    #                 nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
    #                 relax = TRUE, parallel=TRUE, standardize=TRUE,
    #                 type.measure='mse', penalty.factor=svm_penalty)

    beta_svm_raw = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_svm_allvar = as.matrix(beta_svm_raw[1:p]) #remove intercept to the last item


    svm_coef = beta_svm_allvar

    beta_svm_non_zero = row(beta_svm_allvar)[which(abs(beta_svm_allvar) >= 0.1)] # check
    sel.var.index  = beta_svm_non_zero

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    xx_true = xx[, beta_svm_non_zero]
    svm_sigmoid.true.var.list = names(xx_true)
    treat.form = formula(paste("A~",paste(svm_sigmoid.true.var.list,collapse="+")))
    mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                      ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                      replace = FALSE, discard = "both"))

    if ('try-error' %in% class(mm)){
      att.svm_sigmoid[itr] = NA
      select.var.list.svm_sigmoid[itr, 1:length(sel.var.index)] = 0
      total.num.var.svm_sigmoid [itr] = NA
      fail_times[label] = fail_times[label]+1
    } else{
      effect.form = formula(paste("Y~A+",paste(svm_sigmoid.true.var.list,collapse="+")))
      fit.effect = lm(effect.form, data = match.data(mm))

      ### post record
      final.coef = fit.effect$coefficients
      # att
      att.svm_sigmoid[itr] = final.coef[2]
      #Storing the Indexes of selected variables
      select.var.list.svm_sigmoid[itr, 1:length(sel.var.index)] = sel.var.index
      #Storing the totoal number of variables selected
      total.num.var.svm_sigmoid [itr] = length(sel.var.index)
      # Percent var select
      for (v in sel.var.index){
        Var_Select[label, v] <- Var_Select[label, v]+1
      }

      print("att.svm_sigmoid")
      print(att.svm_sigmoid[itr])
      print("select.var.list.svm_sigmoid")
      print(sel.var.index)
    }




    ########### 9. lr_sigmoid ##########
    label = label+1
    ts = Sys.time()

    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)

    weight_logit = coef(fit1)
    abs(weight_logit)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    oddgamma = 1
    penalty = abs(weight_logit) ^ oddgamma
    penalty[is.na(penalty)] <- 10000000000
    penalty = vec_sigmoid(penalty)
    logit_penalty = penalty

    fit2 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = logit_penalty, parallel = TRUE)
    # fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
    #                  nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
    #                  relax = TRUE, parallel=TRUE, standardize=TRUE,
    #                  type.measure='mse', penalty.factor=logit_penalty)

    beta_logit_raw = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_logit_allvar = as.matrix(beta_logit_raw[1:p]) #remove intercept to the last item

    lr_coef = beta_logit_allvar

    beta_logit_non_zero = row(beta_logit_allvar)[which(abs(beta_logit_allvar) >= 0.1)] # check
    sel.var.index  = beta_logit_non_zero

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    xx_true = xx[, beta_logit_non_zero]
    lr_sigmoid.true.var.list = names(xx_true)
    treat.form = formula(paste("A~",paste(lr_sigmoid.true.var.list,collapse="+")))
    mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                            ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                            replace = FALSE, discard = "both"))

    if ('try-error' %in% class(mm)){
      att.lr_sigmoid[itr] = NA
      select.var.list.lr_sigmoid[itr, 1:length(sel.var.index)] = 0
      total.num.var.lr_sigmoid [itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      effect.form = formula(paste("Y~A+",paste(lr_sigmoid.true.var.list,collapse="+")))
      fit.effect = lm(effect.form, data = match.data(mm))

      ### post record
      final.coef = fit.effect$coefficients
      # att
      att.lr_sigmoid[itr] = final.coef[2]
      #Storing the Indexes of selected variables
      select.var.list.lr_sigmoid[itr, 1:length(sel.var.index)] = sel.var.index
      #Storing the totoal number of variables selected
      total.num.var.lr_sigmoid[itr] = length(sel.var.index)
      # Percent var select
      for (v in sel.var.index){
        Var_Select[label, v] <- Var_Select[label, v]+1
      }

      print("att.lr_sigmoid")
      print(att.lr_sigmoid[itr])
      print("select.var.list.lr_sigmoid")
      print(sel.var.index)
    }



#     # ########### 7.1 Comparing the penalty attribute to each variable of svm and lr ######
#     # # number of negative coefficients selected by svm
#     # # penalty_comp = vec_norm(abs(svm_penalty)) - vec_norm(abs(logit_penalty)) # note: record together with rho
#     # # Comp_svm_lr_Percent[k, 1] = Comp_svm_lr_Percent[k, 1]+sum(penalty_comp < 0) # total negative
#     # # Comp_svm_lr_Percent[k, 2] = Comp_svm_lr_Percent[k, 2]+sum(penalty_comp[1:pC] < 0) # confounder chosen negative
#     # # Comp_svm_lr_Percent[k, 3] = Comp_svm_lr_Percent[k, 3]+sum(penalty_comp[pC+1:pC+pP] < 0) # outcome predictor negative
#     #
#     # coef_comp = lr_coef - svm_coef # note: record together with rho
#     # # negative -> lr_coef smaller than svm-coef
#     #
#     # # if for pI, over 50%, svm is better; below 50%, lr is better
#     # # if for pC and pP, over 50%,  svm is better; below 50%, lr is better
#     # Comp_svm_lr_Percent[k, 1] = Comp_svm_lr_Percent[k, 1]+sum(coef_comp[pC+pP+1:pC+pP+pI] >= 0) # treatment penalty chosen negative
#     # # non_zero_T = non_zero_T + sum(coef_comp[pC+pP+1:pC+pP+pI] != 0)
#     # Comp_svm_lr_Percent[k, 2] = Comp_svm_lr_Percent[k, 2]+sum(coef_comp[1:pC] < 0) # confounder chosen negative
#     # Comp_svm_lr_Percent[k, 3] = Comp_svm_lr_Percent[k, 3]+sum(coef_comp[pC+1:pC+pP] < 0) # outcome predictor negative
# 
# 
    ########### 10. lr_bare  ##########
    label = label+1
    ts = Sys.time()

    fit1 = glm(A ~.,family=binomial(link='logit'),data=xa)

    weight_logit = coef(fit1)
    abs(weight_logit)
    weight_logit = weight_logit[-c(1)] # exclude the intercept
    oddgamma = 1
    penalty = abs(weight_logit) ^ oddgamma
    penalty[is.na(penalty)] <- 10000000000
    logit_penalty = penalty

    fit2 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
                 nfolds = 5L, rule = "lambda.1se",
                 penalty.factor = penalty, parallel = TRUE)

    beta_logit_raw = as.matrix(coef(fit2, s=fit2$lambda.1se))
    beta_logit_allvar = as.matrix(beta_logit_raw[1:p]) #remove intercept to the last item
    beta_logit_non_zero = row(beta_logit_allvar)[which(abs(beta_logit_allvar) >= 0.1)] # check
    sel.var.index  = beta_logit_non_zero

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    xx_true = xx[, beta_logit_non_zero]
    lr_bare.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(lr_bare.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.lr_bare[itr] = NA
      select.var.list.lr_bare[itr, 1:length(sel.var.index)] = 0
      total.num.var.lr_bare[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))

      if ('try-error' %in% class(mm)){
        att.lr_bare[itr] = NA
        select.var.list.lr_bare[itr, 1:length(sel.var.index)] = 0
        total.num.var.lr_bare[itr] = NA
        fail_times[label] = fail_times[label]+1
      }else{
        effect.form = formula(paste("Y~A+",paste(lr_bare.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))

    #     ### post record
    #     final.coef = fit.effect$coefficients
    #     # att
    #     att.lr_bare[itr] = final.coef[2]
    #     #Storing the Indexes of selected variables
    #     select.var.list.lr_bare[itr, 1:length(sel.var.index)] = sel.var.index
    #     #Storing the totoal number of variables selected
    #     total.num.var.lr_bare[itr] = length(sel.var.index)
    #     # Percent var select
    #     for (v in sel.var.index){
    #       Var_Select[label, v] <- Var_Select[label, v]+1
    #     }
    # 
    #     print("att.lr_bare")
    #     print(att.lr_bare[itr])
    #     print("select.var.list.lr_bare")
    #     print(sel.var.index)
    #   }
    # }
    # # 
    # # 
    # # 
    # # 
    # # 
    # # 
    # # 
    # # 
    # # 

    ########### 11. adaptive lasso ##########
    label = label+1
    ts = Sys.time()

    # create a dummy dataframe
    dummy_df = Data

    coeff_XA = as.data.frame(matrix(NA,nrow=p,ncol=length(lambda_vec)))
    names(coeff_XA) = names(lambda_vec)
    rownames(coeff_XA) = var.list
    wAMD_vec = rep(NA, length(lambda_vec))

    w.full.form = formula(paste("A~",paste(var.list,collapse="+")))

    for (lil in names(lambda_vec) ){
      il = lambda_vec[lil]
      ig = gamma_vals[lil]

      penalty = (abs(betaXY))^(-ig)
      penalty[which(!is.finite(penalty))] <- 10000000000

      ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
      # fit2 = glmnet(data.matrix(xx), as.matrix(aa), alpha = 1,
      #               gamma = c(0, 0.25, 0.5, 0.75, 1),
      #               penalty.factor=penalty,
      #               lambda=n^(il),
      #               family="binomial")


      fit2 = lqa.formula( w.full.form, data=data.frame(xa), penalty=oal_pen, family=binomial(logit) )
      # logit_oal=lqa(x=data.matrix(xx), y=as.matrix(aa), family = binomial,penalty = oal_pen)

      # generate propensity score
      # dummy_df[,paste("f.pA",lil,sep="")] = predict(fit2, data.matrix(xx))
      dummy_df[,paste("f.pA",lil,sep="")] = predict(fit2)
      # dummy_df = subset(dummy_df, select = -c("id") )
      coeff_XA[var.list,lil] = coef(fit2)[2:(p+1)]
      dummy_df[,paste("w",lil,sep="")] = create_weights(fp=dummy_df[,paste("f.pA",lil,sep="")],fA=aa)
      wAMD_vec[lil] = wAMD_function(DataM=dummy_df,
                                    varlist=var.list,
                                    trt.var="A",
                                    wgt=paste("w",lil,sep=""),
                                    beta=betaXY)$wAMD
    }

    wAMD_vec
    # find the lambda value that creates the smallest wAMD
    tt = which.min(wAMD_vec)

    beta_ps_zero = coeff_XA[,names(tt)]
    beta_ps_allvar = as.matrix(vec_norm(coeff_XA[,names(tt)]))

    beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.01)]
    sel.var.index = beta_ps_non_zero
    sel.var.index

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)+gt

    xx_true = xx[,beta_ps_non_zero]
    adl.true.var.list = names(xx_true)

    # change here
    treat.form <- try(formula(paste("A~",paste(adl.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.adl[itr] = NA
      select.var.list.adl[itr, 1:length(sel.var.index)] = 0
      total.num.var.adl[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))

      if ('try-error' %in% class(mm)){
        att.adl[itr] = NA
        select.var.list.adl[itr, 1:length(sel.var.index)] = 0
        total.num.var.adl[itr] = NA
        fail_times[label] = fail_times[label]+1
      }
      else{
        effect.form = formula(paste("Y~A+",paste(adl.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))

        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.adl[itr] = final.coef[2]
        #Storing the Indexes of selected variables
        select.var.list.adl[itr, 1:length(sel.var.index)] = sel.var.index
        #Storing the totoal number of variables selected
        total.num.var.adl[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }

        print("select.var.list.adl")
        print(sel.var.index)
        print("att.adl")
        print(att.adl[itr])

      }
    }






    ########## 12. outcome adaptive elastic net ##########
    label = label+1
    ts = Sys.time()


    # create a dummy dataframe
    dummy_df = Data

    coeff_XA = as.data.frame(matrix(NA,nrow=p,ncol=length(lambda_vec)))
    names(coeff_XA) = names(lambda_vec)
    rownames(coeff_XA) = var.list
    wAMD_vec = rep(NA, length(lambda_vec))

    w.full.form = formula(paste("A~",paste(var.list,collapse="+")))

    for ( lil in names(lambda_vec) ){
      il = lambda_vec[lil]
      ig = gamma_vals[lil]

      # penalty = (abs(betaXY))^(-ig)
      prearray = log(abs(betaXY)^(-ig))
      prearray = ifelse(prearray < 0, 0, prearray)
      prearray = vec_sigmoid(vec_norm(prearray)*p)
      prearray = ifelse(prearray > 0.9, 10^3, prearray)

      penalty=prearray

      penalty[which(!is.finite(penalty))] <- 10000000000

      # lambda2_seq = c(0, 0.9, by = 0.1) # search overall possible lambda_2

      bl2 = 0
      wAMD_vec[lil] = Inf

      # lambda_seq = n^seq(0.9, 0.1, by=-0.1)
      ### create the outcome adaptive elastic net penalty with coefficient specific weights determined by outcome model
      fit2 <- try(cv.glmnet(data.matrix(xx), as.matrix(aa), alpha = 0.5,
              # lambda = lambda_seq,
              nfolds = 5,
              penalty.factor=penalty,
              family="binomial"))

      # fit2 <- try(aenet(data.matrix(xx), data.matrix(yy), alphas = 0.5,
      #              nfolds = 5L, rule = "lambda.1se",
      #              penalty.factor = penalty, parallel = TRUE))

      if ('try-error' %in% class(fit2)){ # convergence error occurs
        tmp = Inf
      }
      else{
        # generate propensity score
        dummy_df[,paste("f.pA",lil,sep="")] = predict(fit2, data.matrix(xx), s=n^(il))
        coeff_XA[var.list,lil] = coef(fit2, s=n^(il))[2:(p+1)]

        dummy_df[,paste("w",lil,sep="")] = create_weights(fp=dummy_df[,paste("f.pA",lil,sep="")],fA=aa)
        tmp = wAMD_function(DataM=dummy_df,
                            varlist=var.list,
                            trt.var="A",
                            wgt=paste("w",lil,sep=""),
                            beta=betaXY)$wAMD
      }
      wAMD_vec[lil] = min(wAMD_vec[lil], tmp)

    }

    # print(wAMD_vec)
    # find the lambda value that creates the smallest wAMD
    tt = which.min(wAMD_vec)

    beta_ps_zero = coeff_XA[,names(tt)]
    beta_ps_allvar = as.matrix(vec_norm(coeff_XA[,names(tt)]))

    beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.01)]
    sel.var.index = beta_ps_non_zero
    sel.var.index

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)+gt

    xx_true = xx[,beta_ps_non_zero]
    onet.true.var.list = names(xx_true)
    treat.form <- try(formula(paste("A~",paste(onet.true.var.list,collapse="+"))))
    if ('try-error' %in% class(treat.form)){
      att.onet[itr] = NA
      select.var.list.onet[itr, 1:length(sel.var.index)] = 0
      total.num.var.onet[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                        ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                        replace = FALSE, discard = "both"))

      if ('try-error' %in% class(mm)){
        att.onet[itr] = NA
        select.var.list.onet[itr, 1:length(sel.var.index)] = 0
        total.num.var.onet[itr] = NA
        fail_times[label] = fail_times[label]+1
      }
      else{
        effect.form = formula(paste("Y~A+",paste(onet.true.var.list,collapse="+")))
        fit.effect = lm(effect.form, data = match.data(mm))

        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.onet[itr] = final.coef[2]
        #Storing the Indexes of selected variables
        select.var.list.onet[itr, 1:length(sel.var.index)] = sel.var.index
        #Storing the totoal number of variables selected
        total.num.var.onet[itr] = length(sel.var.index)
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }

        print("select.var.list.onet")
        print(sel.var.index)
        print("att.onet")
        print(att.onet[itr])
      }
    }





    ############### 13. BACR ##########
    label = label+1

    ts = Sys.time()

    bacr = bac(data=Data, exposure="A", outcome="Y", confounders=paste(var.list, sep = ""),
               interactors=NULL, familyX="binomial", familyY="gaussian", omega=Inf,
               num_its=200, burnM=1, burnB=1, thin=1)
    bacr_df<-data.frame(bacr$models)
    bacr_df<-bacr_df[1,1:dim(xx)[2]]
    sel.var.index = which(bacr_df>0)
    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    bacr.cov = Data[ , sel.var.index]
    n_bacr<-dim(bacr.cov)[2]
    bacr.cov = as.matrix(bacr.cov)
    bacr.data = data.frame((bacr.cov), aa, yy)

    mm <- try(matchit(aa ~ bacr.cov, data = bacr.data, method = "nearest",
                            distance = "glm", ratio = 1, caliper = .25, link = "logit",
                            estimand = "ATT", replace = FALSE, discard = "both"))

    if ('try-error' %in% class(mm)){
      att.BACR[itr] = NA
      select.var.list.bacr[itr, 1:length(sel.var.index)] = 0
      total.num.var.bacr [itr] = NA
      fail_times[label] = fail_times[label]+1

    } else {
      matched_data = (match.data(mm))

      df_b = matched_data[,1:(n_bacr + 2)]
      fit.effect = lm(yy ~ aa + ., data = df_b)

      ### post record
      final.coef = fit.effect$coefficients
      att.BACR[itr] = final.coef[2]
      # Storing the Indexes of selected variables
      select.var.list.bacr[itr, 1:length(sel.var.index)] = sel.var.index
      # Storing the total number of variables selected
      total.num.var.bacr [itr] = n_bacr
      # Percent var select
      for (v in sel.var.index){
        Var_Select[label,v] <- Var_Select[label,v]+1
      }
    }





    ########## 14. BCEE ########
    label = label+1
    ts = Sys.time()

    xxx = unname(as.matrix(xx))
    bcee <- try(GBCEE(as.matrix(aa), as.matrix(yy), xxx, omega = 300*sqrt(n),
                    niter = 5000, family.X = "binomial", family.Y = "gaussian",
                    X1 = 1, X0 = 0, priorX = NA, priorY = NA, maxsize = NA, OR = 20,
                    truncation = c(0.01, 0.99), var.comp = "asymptotic", B = 200), silent=TRUE)
    if ('try-error' %in% class(bcee)){
      att.BCEE[itr] = NA
      select.var.list.bcee[itr, 1:length(sel.var.index)] = 0
      total.num.var.bcee [itr] = NA
      fail_times[label] = fail_times[label]+1
    } else{
      new_df1<-data.frame(bcee$models.X)
      new_df2<-data.frame(bcee$models.Y)
      new_df<-rbind(new_df1[1,1:dim(xxx)[2]],new_df2[1,1:dim(xxx)[2]])

      sel.var.index = which(new_df[2,] > 0)

      before = t_elapse[label-3, k]
      ela = Sys.time()-ts
      t_elapse[label-3, k] = before+ela
      print(paste("BCEE time:", ela))

      bcee.cov = Data[, sel.var.index]

      n_bcee = dim(bcee.cov)[2]
      bcee.cov = as.matrix(bcee.cov)
      df_bcee = data.frame(bcee.cov, aa, yy)
      mm = matchit(aa ~ bcee.cov, data = df_bcee, method = "nearest",
                   distance = "glm", ratio = 1, caliper = .25, link = "logit",
                   estimand = "ATT", replace = FALSE, discard = "both")

      matched_data = (match.data(mm))

      df_b = matched_data[,1:(n_bcee + 2)]
      fit.effect = lm(yy ~ aa + ., data = df_b)
      final.coef = fit.effect$coefficients
      att.BCEE[itr] = final.coef[2]

      #Storing the Indexes of selected variables
      select.var.list.bcee[itr, 1:length(sel.var.index)] = sel.var.index
      #Storing the total number of variables selected
      total.num.var.bcee[itr] = n_bcee
      #att.BCEE[itr] = bcee$beta
      # Percent var select
      for (v in sel.var.index){
        Var_Select[label, v] <- Var_Select[label, v]+1
      }
    }







    ######### 15. Boruta_T ##########
    label = label+1
    ts = Sys.time()

    bor<-Boruta(xx,aa,pValue = 0.10)
    bor_df<-data.frame(bor$finalDecision)
    bor_df$Index<-1:dim.data.frame(bor_df)[1]
    bor_df<-bor_df[bor_df$bor.finalDecision=="Confirmed",]
    sel.var.index = bor_df$Index

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    bort.cov = Data[ , sel.var.index]
    n_bort<-dim(bort.cov)[2]
    bort.cov = as.matrix(bort.cov)
    bor.data = data.frame(bort.cov, aa, yy)
    mm <- try(matchit(aa ~ bort.cov, data = bor.data, method = "nearest",
                      distance = "glm", ratio = 1, caliper = .25, link = "logit",
                      estimand = "ATT", replace = FALSE, discard = "both"))

    if ('try-error' %in% class(mm)){
      att.Boruta_T[itr] = NA
      select.var.list.Boruta_T[itr, 1:length(sel.var.index)] = 0
      total.num.var.Boruta_T[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      matched_data = (match.data(mm))

      df_b = matched_data[,1:(n_bort + 2)]
      fit.effect = lm(yy ~ aa + . , data = df_b)

      ### post record
      final.coef = fit.effect$coefficients
      att.Boruta_T[itr] = final.coef[2]
      #Storing the Indexes of selected variables
      select.var.list.Boruta_T[itr, 1:length(sel.var.index)] = sel.var.index
      #Storing the totol number of variables selected
      total.num.var.Boruta_T[itr] = n_bort
      for (v in sel.var.index){
        Var_Select[label, v] <- Var_Select[label, v]+1
      }
    }


    ######### 16. Boruta_Y #######
    label = label+1
    ts = Sys.time()
    bor<-Boruta(xx,yy,pValue = 0.10)

    bor_df<-data.frame(bor$finalDecision)
    bor_df$Index<-1:dim.data.frame(bor_df)[1]
    bor_df<-bor_df[bor_df$bor.finalDecision=="Confirmed",]
    sel.var.index = bor_df$Index

    t_elapse[label-3, k] = t_elapse[label-3, k]+as.numeric(Sys.time()-ts)

    bory.cov = Data[ , sel.var.index]
    n_bory<-dim(bory.cov)[2]
    bory.cov = as.matrix(bory.cov)
    bor.y.data = data.frame(bory.cov, aa, yy)
    mm <- try(matchit(aa ~ bory.cov, data = bor.y.data, method = "nearest",
                      distance = "glm", ratio = 1, caliper = .25, link = "logit",
                      estimand = "ATT", replace = FALSE, discard = "both"))

    if ('try-error' %in% class(mm)){
      att.Boruta_Y[itr] = NA
      select.var.list.Boruta_Y[itr, 1:length(sel.var.index)] = 0
      total.num.var.Boruta_Y[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      matched_data = (match.data(mm))

      df_by = matched_data[,1:(n_bory + 2)]
      fit.effect = lm(yy ~ aa + . , data = df_by)

      ### post record
      final.coef = fit.effect$coefficients
      # att
      att.Boruta_Y[itr] = final.coef[2]
      #Storing the Indexes of selected variables
      select.var.list.Boruta_Y[itr, 1:length(sel.var.index)] = sel.var.index
      #Storing the total number of variables selected
      total.num.var.Boruta_Y[itr] = n_bory
      # Percent var select
      for (v in sel.var.index){
        Var_Select[label, v] <- Var_Select[label, v]+1
      }
    }

    print(t_elapse[,k])









    ######### 17. Deluna(DWR) #######
    label = label+1
    ts = Sys.time()

    sink("file")
    fit1 <- try(cov.sel(aa, yy, xx))
    sink()
    if ('try-error' %in% class(fit1)){
      att.DWR[itr] = NA
      select.var.list.DWR[itr, 1:length(sel.var.index)] = 0
      total.num.var.DWR [itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      dwr_decision = fit1$covar
      sel.var.index = which(colnames(xx)==dwr_decision)

      t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)

      dwr.cov = Data[ , sel.var.index]
      n_dwr<-dim(dwr.cov)[2]
      dwr.cov = as.matrix(dwr.cov)
      dwr.data = data.frame(dwr.cov, aa, yy)
      mm <- try(matchit(aa ~ dwr.cov, data = dwr.data, method = "nearest",
                        distance = "glm", ratio = 1, caliper = .25, link = "logit",
                        estimand = "ATT", replace = FALSE, discard = "both"))

      if ('try-error' %in% class(mm)){
        att.DWR[itr] = NA
        select.var.list.DWR[itr, 1:length(sel.var.index)] = 0
        total.num.var.DWR [itr] = NA
        fail_times[label] = fail_times[label]+1
      }else{
        matched_data = (match.data(mm))

        df_dwr = matched_data[,1:(n_dwr + 2)]
        fit.effect = lm(yy ~ aa + . , data = df_dwr)

        ### post record
        final.coef = fit.effect$coefficients
        # att
        att.DWR[itr] = final.coef[2]
        #Storing the Indexes of selected variables
        select.var.list.DWR[itr, 1:length(sel.var.index)] = sel.var.index
        #Storing the total number of variables selected
        total.num.var.DWR[itr] = n_dwr
        # Percent var select
        for (v in sel.var.index){
          Var_Select[label, v] <- Var_Select[label, v]+1
        }

        print("select.var.list.DWR")
        print(sel.var.index)
        print("att.dwr")
        print(att.DWR[itr])
      }
    }

    ######### 18. CTMLE(Collaborative Targeted Maximum Likelihood Estimation) #######
    label = label+1
    ts = Sys.time()

    Q <- cbind(rep(mean(yy[aa == 0]), n), rep(mean(yy[aa == 1]), n))
    fit1 <- try(ctmleDiscrete(Y = yy, A = aa, W = xx, Q = Q, preOrder = TRUE, detailed = TRUE))
    if ('try-error' %in% class(fit1)){
      att.CTMLE[itr] = NA
      select.var.list.CTMLE[itr, 1:length(sel.var.index)] = 0
      total.num.var.CTMLE[itr] = NA
      fail_times[label] = fail_times[label]+1
    }else{
      ctmle.cov = fit1$candidates$terms[-c(1)]
      sel.var.index = which(colnames(xx)==ctmle.cov)
      t_elapse[label-3, k] = t_elapse[label-3, k] +as.numeric(Sys.time()-ts)

      ctmle.cov = Data[ , sel.var.index]
      n_ctmle<-dim(ctmle.cov)[2]
      ctmle.cov = as.matrix(ctmle.cov)
      ctmle.data = data.frame(ctmle.cov, aa, yy)
      mm <- try(matchit(aa ~ ctmle.cov, data = ctmle.data, method = "nearest",
                        distance = "glm", ratio = 1, caliper = .25, link = "logit",
                        estimand = "ATT", replace = FALSE, discard = "both"))

      if ('try-error' %in% class(mm)){
        att.CTMLE[itr] = NA
        select.var.list.CTMLE[itr, 1:length(sel.var.index)] = 0
        total.num.var.CTMLE[itr] = NA
        fail_times[label] = fail_times[label]+1
      }else{
        matched_data = (match.data(mm))

        if (is.null(n_ctmle)){
          att.CTMLE[itr] = NA
          select.var.list.CTMLE[itr, 1:length(sel.var.index)] = 0
          total.num.var.CTMLE[itr] = NA
          fail_times[label] = fail_times[label]+1
        }
        else{
          df_ctmle <- matched_data[,1:(n_ctmle + 2)]
          fit.effect = lm(yy ~ aa + . , data = df_ctmle)

          ### post record
          final.coef = fit.effect$coefficients
          # att
          att.CTMLE[itr] = final.coef[2]
          #Storing the Indexes of selected variables
          select.var.list.CTMLE[itr, 1:length(sel.var.index)] = sel.var.index
          #Storing the total number of variables selected
          total.num.var.CTMLE[itr] = n_ctmle
          # Percent var select
          for (v in sel.var.index){
            Var_Select[label, v] <- Var_Select[label, v]+1
          }

          print("select.var.list.CTMLE")
          print(sel.var.index)
          print("att.ctmle")
          print(att.CTMLE[itr])
        }
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
  df_att = data.frame(att.target, att.conf, att.potconf,
                      att.Enh_lr_tanh, att.Enh_lr_sigmoid,
                      att.Enh_svm_tanh, att.Enh_svm_sigmoid,
                      att.svm_sigmoid,
                      att.lr_sigmoid, att.lr_bare,
                      att.adl,att.onet,
                      att.BACR,
                      att.BCEE,
                      att.Boruta_T,
                      att.Boruta_T,
                      att.DWR,
                      att.CTMLE)

  # save num of vars selected
  df_nsel.var = data.frame(total.num.var.Enh_lr_tanh,
                           total.num.var.Enh_lr_sigmoid,
                           total.num.var.Enh_svm_tanh,
                           total.num.var.Enh_svm_sigmoid,
                           total.num.var.svm_sigmoid,
                           total.num.var.lr_sigmoid,
                           total.num.var.lr_bare,
                           total.num.var.adl,
                           total.num.var.onet,
                           total.num.var.bacr,
                           total.num.var.bcee,
                           total.num.var.Boruta_T,
                           total.num.var.Boruta_Y,
                           total.num.var.DWR,
                           total.num.var.CTMLE)

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




