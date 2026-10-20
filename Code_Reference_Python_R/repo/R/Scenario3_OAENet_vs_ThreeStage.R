## check_scenario3_x1_vs_x2.R
## Scenario 3 (weak confounders), same seeds and data for every method.
## The data generation, Enh_lr_tanh, Enh_svm_sigmoid and ONet blocks are copied from FS_CI_Final.R;
## the only changes are (1) the coefficient multiplier MULT is a loop variable instead of a hard-coded 2,
## (2) p = 20 as in the paper, and (3) an extra ONet variant with predict(type = "response").
## Output: selection rate of every variable and the percentage selection bias |ATT - ATT_target| / 0.5.

library(msaenet); library(glmnet); library(MatchIt); library(MASS); library(e1071)

## ---------------- settings ----------------
MULT_LIST = c(1, 2)                 # 1 = coefficients as printed in Section 6.1, 2 = as in FS_CI_Final.R
rho_list  = c(0, 0.25, 0.5, 0.75)
itr_max   = 30
n = 200
p = 20
OUT_FILE  = "scenario3_x1_vs_x2.csv"

## ---------------- scenario 3 (from FS_CI_Final.R) ----------------
gd.var = 7; pC = 2; pI = 2; pP = 2
beta_base  = c( 0.2, 0.2, 0.6, 0.6, 0, 0, 0, rep(0, p - gd.var) )
alpha_base = c( 1.0, 1.0, 0, 0, 1.0, 1.0, 0, rep(0, p - gd.var) )
pS = p - (pC+pI+pP)
var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""),paste("Xi",1:pI,sep=""),paste("Xs",1:pS,sep=""))
target.var.list = c(paste("Xc",1:pC,sep=""),paste("Xp",1:pP,sep=""))
bA = 0.5; mean_x = 0; sig_x = 1

expit = function(x) exp(x) / (1+exp(x))
vec_norm = function(x) x / sqrt(sum(x^2))
vec_sigmoid <- function(pen) 1/(1+exp(-pen))
vec_tanh <- function(pen) tanh(pen)

lambda_vec = c(-10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
names(lambda_vec) = as.character(lambda_vec)
gamma_convergence_factor = 2
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
  list(wAMD = sum(diff_vec * abs(beta)))
}
create_weights = function(fp,fA){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  fw
}

## matching + ATT exactly as in FS_CI_Final.R (returns NA on the same failures)
match_att = function(sel.idx, Data, xx){
  xx_true = xx[, sel.idx]
  vl = names(xx_true)
  treat.form <- try(formula(paste("A~",paste(vl,collapse="+"))), silent = TRUE)
  if ('try-error' %in% class(treat.form)) return(NA)
  mm <- try(matchit(treat.form, data = Data, method = "nearest", distance = "glm",
                    ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
                    replace = FALSE, discard = "both"), silent = TRUE)
  if ('try-error' %in% class(mm)) return(NA)
  effect.form = formula(paste("Y~A+",paste(vl,collapse="+")))
  lm(effect.form, data = match.data(mm))$coefficients[2]
}

## ---------------- methods (blocks from FS_CI_Final.R) ----------------
enh = function(xx, xa, yy, prior, smoother){
  if (prior == "lr"){
    fit1 = glm(A ~., family=binomial(link='logit'), data=xa)
  } else {
    fit1 = svm(A ~., data = xa, type = 'C-classification', kernel = 'linear', scale = TRUE)
  }
  weight_logit = coef(fit1)[-c(1)]
  penalty = abs(weight_logit)
  penalty[is.na(penalty)] <- 10000000000
  penalty <- if (smoother == "tanh") vec_tanh(penalty)^0.5 else vec_sigmoid(penalty)^1
  fit2 = cv.glmnet(data.matrix(xx), data.matrix(yy), alpha = 0.5,
                   nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                   relax = TRUE, parallel=TRUE, standardize=TRUE,
                   type.measure='mse', penalty.factor=penalty)
  beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
  beta_el_allvar = as.matrix(beta_el[2:(p+1)])
  penalty_adel = (1/abs(beta_el_allvar))^1
  penalty_adel[is.na(penalty_adel)] <- 10000000000
  fit3 = aenet(data.matrix(xx), data.matrix(yy), alphas = seq(0.05, 0.95, 0.05),
               nfolds = 5L, rule = "lambda.1se",
               penalty.factor = penalty_adel, parallel = TRUE)
  beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
  beta_aden_allvar = as.matrix(beta_aden_raw[1:p])
  row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)]
}

onet = function(Data, xx, xa, aa, betaXY, ps_type){      # ps_type = "link" (as in FS_CI_Final.R) or "response"
  dummy_df = Data
  coeff_XA = as.data.frame(matrix(NA,nrow=p,ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec); rownames(coeff_XA) = var.list
  wAMD_vec = rep(NA, length(lambda_vec)); names(wAMD_vec) = names(lambda_vec)
  for (lil in names(lambda_vec)){
    il = lambda_vec[lil]; ig = gamma_vals[lil]
    prearray = log(abs(betaXY)^(-ig))
    prearray = ifelse(prearray < 0, 0, prearray)
    prearray = vec_sigmoid(vec_norm(prearray)*p)
    prearray = ifelse(prearray > 0.9, 10^3, prearray)
    penalty = prearray
    penalty[which(!is.finite(penalty))] <- 10000000000
    wAMD_vec[lil] = Inf
    fit2 <- try(cv.glmnet(data.matrix(xx), as.matrix(aa), alpha = 0.5, nfolds = 5,
                          penalty.factor=penalty, family="binomial"), silent = TRUE)
    if ('try-error' %in% class(fit2)) { tmp = Inf } else {
      dummy_df[,paste("f.pA",lil,sep="")] = predict(fit2, data.matrix(xx), s=n^(il), type = ps_type)
      coeff_XA[var.list,lil] = coef(fit2, s=n^(il))[2:(p+1)]
      dummy_df[,paste("w",lil,sep="")] = create_weights(fp=dummy_df[,paste("f.pA",lil,sep="")],fA=aa)
      tmp = wAMD_function(DataM=dummy_df, varlist=var.list, trt.var="A",
                          wgt=paste("w",lil,sep=""), beta=betaXY)$wAMD
    }
    wAMD_vec[lil] = min(wAMD_vec[lil], tmp)
  }
  tt = which.min(wAMD_vec)
  beta_ps_allvar = as.matrix(vec_norm(coeff_XA[,names(tt)]))
  row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.01)]
}

## ---------------- simulation ----------------
methods = c("Enh_lr_tanh", "Enh_svm_sigmoid", "ONet_link(as_in_file)", "ONet_response")
rows = list()
for (MULT in MULT_LIST){
  beta_v = MULT * beta_base; alpha_v = MULT * alpha_base
  names(beta_v) = names(alpha_v) = var.list
  for (rho in rho_list){
    for (itr in 1:itr_max){
      set.seed(itr)
      Sigma_x = matrix(rho*sig_x^2,nrow=length(var.list),ncol=length(var.list)); diag(Sigma_x) = sig_x^2
      Data = as.data.frame(mvrnorm(n = n, mu=rep(mean_x,length(var.list)), Sigma = Sigma_x))
      names(Data) = var.list
      gA_x = rowSums(Data[,var.list]*matrix(alpha_v,nrow=n,ncol=length(var.list),byrow=TRUE))
      Data$A = as.numeric( runif(n=n) < expit(gA_x) )
      gY_xA = rowSums(Data[,var.list]*matrix(beta_v,nrow=n,ncol=length(var.list),byrow=TRUE))
      Data$Y = gY_xA + rnorm(n=n,sd=sig_x) + Data$A*bA
      Data[,var.list] = scale(Data[,var.list])
      xx = Data[,1:p]; xa = Data[,1:(p+1)]; yy = Data[,(p+2)]; aa = Data[,(p+1)]
      betaXY = coef(lm(formula(paste("Y~A+",paste(var.list,collapse="+"))), data=Data))[var.list]
      
      att_target = match_att(1:(pC+pP), Data, xx)
      sels = list(
        enh(xx, xa, yy, "lr", "tanh"),
        enh(xx, xa, yy, "svm", "sigmoid"),
        onet(Data, xx, xa, aa, betaXY, "link"),
        onet(Data, xx, xa, aa, betaXY, "response"))
      for (m in seq_along(methods)){
        sel = sels[[m]]
        att = match_att(sel, Data, xx)
        r = data.frame(mult = MULT, rho = rho, itr = itr, method = methods[m],
                       n_selected = length(sel), att = att, att_target = att_target,
                       pct_bias = 100 * abs(att - att_target) / bA)
        for (j in 1:p) r[[var.list[j]]] = as.integer(j %in% sel)
        rows[[length(rows)+1]] = r
      }
      cat(sprintf("mult=%g rho=%g itr=%d done\n", MULT, rho, itr))
    }
    write.csv(do.call(rbind, rows), OUT_FILE, row.names = FALSE)
  }
}

## ---------------- summary ----------------
res = do.call(rbind, rows)
res$XC = rowMeans(res[, c("Xc1","Xc2")]); res$XP = rowMeans(res[, c("Xp1","Xp2")])
res$XT = rowMeans(res[, c("Xi1","Xi2")]); res$XS = rowMeans(res[, paste("Xs",1:pS,sep="")])
summ = aggregate(cbind(XC, XP, XT, XS, pct_bias) ~ mult + rho + method, data = res,
                 FUN = function(x) round(mean(x, na.rm = TRUE), 3), na.action = na.pass)
fails = aggregate(att ~ mult + rho + method, data = res, FUN = function(x) sum(is.na(x)), na.action = na.pass)
names(fails)[4] = "n_failed"
summ = merge(summ, fails)
summ = summ[order(summ$mult, summ$rho, summ$method), ]
print(summ, row.names = FALSE)
cat("\nPooled over rho:\n")
print(aggregate(cbind(XC, XT, pct_bias) ~ mult + method, data = res,
                FUN = function(x) round(mean(x, na.rm = TRUE), 3), na.action = na.pass), row.names = FALSE)