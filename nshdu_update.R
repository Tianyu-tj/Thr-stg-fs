########################################
library(msaenet)
library(glmnet)
library(MatchIt)
library(lmtest)
library(phonTools)
library(MASS)
library(sandwich)
library(CovSel)
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
library(dplyr)

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

########################################################
workspace = "~/RProjects/"
raw_data = read.csv(paste(workspace,"2015-19_OUD_data1.csv",sep=""))


##### data preprocessing #####

#converting the categories into numeric factors

raw_data$year = as.numeric(factor(raw_data$year))
raw_data$year = factor(raw_data$year)


raw_data = subset(raw_data, select = -c(PNRLWD3SX_flag) )

dat_NA_removed = subset(raw_data, dstworst_flag != "Not Available")
dat_NA_removed = subset(dat_NA_removed, c(pnrlndmor_flag, ANYHLTI2_flag,pnrwygamt_flag,pnrwyoftn_flag,pnrwylngr_flag,pnrrspain_flag,
                          pnrrsrelx_flag,pnrrsexpt_flag,pnrrshigh_flag,pnrrsslep_flag,pnrrsemot_flag,pnrrsdgfx_flag,
                          pnrrshook_flag,pnrrssor_flag,pnrllottm_flag,pnrllimit_flag,pnrlcutdn_flag)!= "Not Available")

dat_NA_removed = subset(dat_NA_removed, pnrlndmor_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, ANYHLTI2_flag != "Not Available") 
dat_NA_removed = subset(dat_NA_removed, pnrwygamt_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrwyoftn_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrwylngr_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrspain_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrsrelx_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrsexpt_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrshigh_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrsslep_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrsemot_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrsdgfx_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrshook_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrrssor_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrllottm_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrllimit_flag!= "Not Available")
dat_NA_removed = subset(dat_NA_removed, pnrlcutdn_flag!= "Not Available")

dat_NA_removed = subset(dat_NA_removed, pnrlndmor_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, ANYHLTI2_flag != "Not Availaible") 
dat_NA_removed = subset(dat_NA_removed, pnrwygamt_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrwyoftn_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrwylngr_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrspain_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrsrelx_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrsexpt_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrshigh_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrsslep_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrsemot_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrsdgfx_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrshook_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrrssor_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrllottm_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrllimit_flag!= "Not Availaible")
dat_NA_removed = subset(dat_NA_removed, pnrlcutdn_flag!= "Not Availaible")

dat_NA_removed$CATAG3 = factor(dat_NA_removed$CATAG3)
dat_NA_removed$NEWRACE2 = factor(dat_NA_removed$NEWRACE2)
dat_NA_removed$eduhighcat = factor(dat_NA_removed$eduhighcat)
dat_NA_removed$irwrkstat = factor(dat_NA_removed$irwrkstat)
dat_NA_removed$income = factor(dat_NA_removed$income)
dat_NA_removed$PDEN10 = factor(dat_NA_removed$PDEN10)
dat_NA_removed$HEALTH2 = factor(dat_NA_removed$HEALTH2)

dat_NA_removed$ANYHLTI2_flag = as.numeric(factor(dat_NA_removed$ANYHLTI2_flag))-1
#dat_NA_removed$ANYHLTI2_flag = factor(dat_NA_removed$ANYHLTI2_flag)

dat_NA_removed$pnrwygamt_flag = as.numeric(factor(dat_NA_removed$pnrwygamt_flag))-1
#dat_NA_removed$pnrwygamt_flag = factor(dat_NA_removed$pnrwygamt_flag)

dat_NA_removed$pnrwyoftn_flag = as.numeric(factor(dat_NA_removed$pnrwyoftn_flag))-1
#dat_NA_removed$pnrwyoftn_flag = factor(dat_NA_removed$pnrwyoftn_flag)

dat_NA_removed$pnrwylngr_flag = as.numeric(factor(dat_NA_removed$pnrwylngr_flag))-1
#dat_NA_removed$pnrwylngr_flag = factor(dat_NA_removed$pnrwylngr_flag)

dat_NA_removed$pnrrspain_flag = as.numeric(factor(dat_NA_removed$pnrrspain_flag))-1
#dat_NA_removed$pnrrspain_flag = factor(dat_NA_removed$pnrrspain_flag)

dat_NA_removed$pnrrsrelx_flag = as.numeric(factor(dat_NA_removed$pnrrsrelx_flag))-1
#dat_NA_removed$pnrrsrelx_flag= factor(dat_NA_removed$pnrrsrelx_flag)

dat_NA_removed$pnrrsexpt_flag = as.numeric(factor(dat_NA_removed$pnrrsexpt_flag))-1
#dat_NA_removed$pnrrsexpt_flag = factor(dat_NA_removed$pnrrsexpt_flag)

dat_NA_removed$pnrrshigh_flag = as.numeric(factor(dat_NA_removed$pnrrshigh_flag))-1
#dat_NA_removed$pnrrshigh_flag = factor(dat_NA_removed$pnrrshigh_flag)

dat_NA_removed$pnrrsslep_flag = as.numeric(factor(dat_NA_removed$pnrrsslep_flag))-1
#dat_NA_removed$pnrrsslep_flag = factor(dat_NA_removed$pnrrsslep_flag)

dat_NA_removed$pnrrsemot_flag = as.numeric(factor(dat_NA_removed$pnrrsemot_flag))-1
#dat_NA_removed$pnrrsemot_flag = factor(dat_NA_removed$pnrrsemot_flag)

dat_NA_removed$pnrrsdgfx_flag = as.numeric(factor(dat_NA_removed$pnrrsdgfx_flag))-1
#dat_NA_removed$pnrrsdgfx_flag = factor(dat_NA_removed$pnrrsdgfx_flag)

dat_NA_removed$pnrrshook_flag = as.numeric(factor(dat_NA_removed$pnrrshook_flag))-1
#dat_NA_removed$pnrrshook_flag = factor(dat_NA_removed$pnrrshook_flag)

dat_NA_removed$pnrrssor_flag = as.numeric(factor(dat_NA_removed$pnrrssor_flag))-1
#dat_NA_removed$pnrrssor_flag = factor(dat_NA_removed$pnrrssor_flag)

dat_NA_removed$dstworst_flag = as.numeric(factor(dat_NA_removed$dstworst_flag))-1
#dat_NA_removed$dstworst_flag = factor(dat_NA_removed$dstworst_flag)

dat_NA_removed$pnrllottm_flag = as.numeric(factor(dat_NA_removed$pnrllottm_flag))-1
#dat_NA_removed$pnrllottm_flag = factor(dat_NA_removed$pnrllottm_flag)

dat_NA_removed$pnrllimit_flag = as.numeric(factor(dat_NA_removed$pnrllimit_flag))-1
#dat_NA_removed$pnrllimit_flag = factor(dat_NA_removed$pnrllimit_flag)

dat_NA_removed$pnrlndmor_flag = as.numeric(factor(dat_NA_removed$pnrlndmor_flag))-1
#dat_NA_removed$pnrlndmor_flag = factor(dat_NA_removed$pnrlndmor_flag)

dat_NA_removed$pnrlcutdn_flag = as.numeric(factor(dat_NA_removed$pnrlcutdn_flag))-1
#dat_NA_removed$pnrlcutdn_flag = factor(dat_NA_removed$pnrlcutdn_flag)

################################
#Processing No-OUD data

# raw_data_no_oud = read.csv("C:\\Users\\sahil\\OneDrive - Northeastern University\\NEU\\OneDrive_1_5-19-2022\\2015-19_no_OUD_data1.csv")
raw_data_no_oud = read.csv(paste(workspace,"2015-19_no_OUD_data1.csv",sep=""))

raw_data_no_oud$year = as.numeric(factor(raw_data_no_oud$year))
raw_data_no_oud$year = factor(raw_data_no_oud$year)


dat_NA_removed_no_oud = subset(raw_data_no_oud, dstworst_flag != "Not Available")

#dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, c(pnrlndmor_flag, ANYHLTI2_flag,pnrwygamt_flag,pnrwyoftn_flag,pnrwylngr_flag,pnrrspain_flag,
#                                          pnrrsrelx_flag,pnrrsexpt_flag,pnrrshigh_flag,pnrrsslep_flag,pnrrsemot_flag,pnrrsdgfx_flag,
#                                          pnrrshook_flag,pnrrssor_flag,pnrllottm_flag,pnrllimit_flag,pnrlcutdn_flag)!= "Not Available")

dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrlndmor_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, ANYHLTI2_flag != "Not Available") 
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrwygamt_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrwyoftn_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrwylngr_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrspain_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsrelx_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsexpt_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrshigh_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsslep_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsemot_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsdgfx_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrshook_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrssor_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrllottm_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrllimit_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrlcutdn_flag!= "Not Available")

dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrlndmor_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, ANYHLTI2_flag != "Not Availaible") 
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrwygamt_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrwyoftn_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrwylngr_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrspain_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsrelx_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsexpt_flag!= "Not Available")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrshigh_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsslep_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsemot_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrsdgfx_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrshook_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrrssor_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrllottm_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrllimit_flag!= "Not Availaible")
dat_NA_removed_no_oud = subset(dat_NA_removed_no_oud, pnrlcutdn_flag!= "Not Availaible")




dat_NA_removed_no_oud$CATAG3 = factor(dat_NA_removed_no_oud$CATAG3)
dat_NA_removed_no_oud$NEWRACE2 = factor(dat_NA_removed_no_oud$NEWRACE2)
dat_NA_removed_no_oud$eduhighcat = factor(dat_NA_removed_no_oud$eduhighcat)
dat_NA_removed_no_oud$irwrkstat = factor(dat_NA_removed_no_oud$irwrkstat)
dat_NA_removed_no_oud$income = factor(dat_NA_removed_no_oud$income)
dat_NA_removed_no_oud$PDEN10 = factor(dat_NA_removed_no_oud$PDEN10)
dat_NA_removed_no_oud$HEALTH2 = factor(dat_NA_removed_no_oud$HEALTH2)


#raw_data_no_oud = subset(raw_data_no_oud, select = -c(PNRLWD3SX_flag) )

dat_NA_removed_no_oud$ANYHLTI2_flag = as.numeric(factor(dat_NA_removed_no_oud$ANYHLTI2_flag))-1
#dat_NA_removed_no_oud$ANYHLTI2_flag = factor(dat_NA_removed_no_oud$ANYHLTI2_flag)

dat_NA_removed_no_oud$pnrwygamt_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrwygamt_flag))-1
#dat_NA_removed_no_oud$pnrwygamt_flag = factor(dat_NA_removed_no_oud$pnrwygamt_flag)

dat_NA_removed_no_oud$pnrwyoftn_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrwyoftn_flag))-1
#dat_NA_removed_no_oud$pnrwyoftn_flag = factor(dat_NA_removed_no_oud$pnrwyoftn_flag)

dat_NA_removed_no_oud$pnrwylngr_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrwylngr_flag))-1
#dat_NA_removed_no_oud$pnrwylngr_flag = factor(dat_NA_removed_no_oud$pnrwylngr_flag)

dat_NA_removed_no_oud$pnrrspain_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrspain_flag))-1
#dat_NA_removed_no_oud$pnrrspain_flag = factor(dat_NA_removed_no_oud$pnrrspain_flag)

dat_NA_removed_no_oud$pnrrsrelx_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrsrelx_flag))-1
#dat_NA_removed_no_oud$pnrrsexpt_flag = factor(dat_NA_removed_no_oud$pnrrsexpt_flag)

dat_NA_removed_no_oud$pnrrsexpt_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrsexpt_flag))-1
#dat_NA_removed_no_oud$pnrrsexpt_flag = factor(dat_NA_removed_no_oud$pnrrsexpt_flag)

dat_NA_removed_no_oud$pnrrshigh_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrshigh_flag))-1
#dat_NA_removed_no_oud$pnrrshigh_flag = factor(dat_NA_removed_no_oud$pnrrshigh_flag)

dat_NA_removed_no_oud$pnrrsslep_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrsslep_flag))-1
#dat_NA_removed_no_oud$pnrrsslep_flag = factor(dat_NA_removed_no_oud$pnrrsslep_flag)

dat_NA_removed_no_oud$pnrrsemot_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrsemot_flag))-1
#dat_NA_removed_no_oud$pnrrsemot_flag = factor(dat_NA_removed_no_oud$pnrrsemot_flag)

dat_NA_removed_no_oud$pnrrsdgfx_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrsdgfx_flag))-1
#dat_NA_removed_no_oud$pnrrsdgfx_flag = factor(dat_NA_removed_no_oud$pnrrsdgfx_flag)

dat_NA_removed_no_oud$pnrrshook_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrshook_flag))-1
#dat_NA_removed_no_oud$pnrrshook_flag = factor(dat_NA_removed_no_oud$pnrrshook_flag)

dat_NA_removed_no_oud$pnrrssor_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrssor_flag))-1
#dat_NA_removed_no_oud$pnrrssor_flag = factor(dat_NA_removed_no_oud$pnrrssor_flag)

dat_NA_removed_no_oud$dstworst_flag = as.numeric(factor(dat_NA_removed_no_oud$dstworst_flag))-1
#dat_NA_removed_no_oud$dstworst_flag = factor(dat_NA_removed_no_oud$dstworst_flag)

dat_NA_removed_no_oud$pnrllottm_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrllottm_flag))-1
#dat_NA_removed_no_oud$pnrllottm_flag = factor(dat_NA_removed_no_oud$pnrllottm_flag)

dat_NA_removed_no_oud$pnrllimit_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrllimit_flag))-1
#dat_NA_removed_no_oud$pnrllimit_flag = factor(dat_NA_removed_no_oud$pnrllimit_flag)

dat_NA_removed_no_oud$pnrlndmor_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrlndmor_flag))-1
#dat_NA_removed_no_oud$pnrlndmor_flag = factor(dat_NA_removed_no_oud$pnrlndmor_flag)

dat_NA_removed_no_oud$pnrlcutdn_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrlcutdn_flag))-1
#dat_NA_removed_no_oud$pnrlcutdn_flag = factor(dat_NA_removed_no_oud$pnrlcutdn_flag)

#dat_NA_removed_no_oud$pnrrsexpt_flag = as.numeric(factor(dat_NA_removed_no_oud$pnrrsexpt_flag))-1
#dat_NA_removed_no_oud$pnrrsexpt_flag = factor(dat_NA_removed_no_oud$pnrrsexpt_flag)

###### General Settings ######

itr_max = 500

max_var = 70

# declare a zero matrix to hold att
att.sample = zeros(itr_max)
att.expert = zeros(itr_max)

att.enh_svm_sig = zeros(itr_max)
att.enh_log_tanh = zeros(itr_max)
att.oaenet = zeros(itr_max)
att.oal = zeros(itr_max)
att.bacr = zeros(itr_max)
att.bcee = zeros(itr_max)
att.Boruta_T = zeros(itr_max)
att.Boruta_Y = zeros(itr_max)
att.DWR = zeros(itr_max)
att.CTMLE = zeros(itr_max)

#declare zero matrices to store the indexes of selected variables
select.var.list.enh_svm_sig = zeros(itr_max, max_var)
select.var.list.enh_log_tanh = zeros(itr_max, max_var)
select.var.list.oaenet = zeros(itr_max, max_var)
select.var.list.oal = zeros(itr_max, max_var)
select.var.list.bacr = zeros(itr_max, max_var)
select.var.list.bcee = zeros(itr_max, max_var)
select.var.list.bort = zeros(itr_max, max_var)
select.var.list.bory = zeros(itr_max, max_var)
select.var.list.dwr = zeros(itr_max, max_var)
select.var.list.ctmle = zeros(itr_max, max_var)

# #Declare zero matrices to store the total number of variables selected
total.num.var.enh_svm_sig = zeros(itr_max)
total.num.var.enh_log_tanh = zeros(itr_max)
total.num.var.oaenet = zeros(itr_max)
total.num.var.oal = zeros(itr_max)
total.num.var.bacr = zeros(itr_max)
total.num.var.bcee = zeros(itr_max)
total.num.var.bort = zeros(itr_max)
total.num.var.bory = zeros(itr_max)
total.num.var.dwr = zeros(itr_max)
total.num.var.ctmle = zeros(itr_max)

# #Declare zero matrices to store the total number of variables selected
time.enh_svm_sig = zeros(itr_max)
time.enh_log_tanh = zeros(itr_max)
time.oaenet = zeros(itr_max)
time.oal = zeros(itr_max)
time.bacr = zeros(itr_max)
time.bcee = zeros(itr_max)
time.bort = zeros(itr_max)
time.bory = zeros(itr_max)
time.dwr = zeros(itr_max)
time.ctmle = zeros(itr_max)

#one hot encoding 
df_full = rbind(dat_NA_removed,dat_NA_removed_no_oud)

df_full <- one_hot(as.data.table(df_full))
df_full

#########
#saving one-hot encoded treatment and control 
treatment_onehot = one_hot(as.data.table(dat_NA_removed))
control_onehot = one_hot(as.data.table(dat_NA_removed_no_oud))
write.csv(treatment_onehot, paste(workspace,"Processed/NSDUH_treatment_onehot.csv", sep=""))
write.csv(control_onehot, paste(workspace,"Processed/NSDUH_control_onehot.csv", sep=""))

####### Penalty Smooth Function ##########
vec_sigmoid <- function(pen) {
  pen <- 1/(1+exp(-pen))
}

vec_tanh <- function(pen){
  pen <- tanh(pen)
}


#########
for (itr in 1:1){
  
  print(paste("Iteration ", itr, sep=""))
  #combining both OUD and No-Oud data
  dat_NA_removed_no_oud_sample = sample_n(dat_NA_removed_no_oud, 5000)
  combined_data = rbind(dat_NA_removed,dat_NA_removed_no_oud_sample)


  ########################
  # Preparing the covariates, outcomes, and treat data frame
  
  df = combined_data
  Y = df$suicide_flag
  A = df$udpypnr
  
  ############################
  #one-hot encoding
  df <- one_hot(as.data.table(df))
  df
  ############################
  
  XA = subset(df, select = -c(suicide_flag)) # predict treatment
  XX = subset(XA, select = -c(udpypnr)) # predict outcome
  
  ##### Sample var no select ####
  sample.var.index = 1:(ncol(df)-2)
  sample.var.list = names(XX[, ..sample.var.index])
  treat.form = formula(paste("udpypnr~",paste(sample.var.list, collapse="+")))
  ## performing matching on df_full instead of df
  mm = matchit(treat.form, data = df, method = "nearest", distance = "glm", ratio = 1,
               caliper = .25, link = "logit", estimand = "ATT", replace = FALSE, discard = "control")
  effect.form = formula(paste("suicide_flag ~ udpypnr+",paste(sample.var.list,collapse="+")))
  fit.effect = lm(effect.form, data = match.data(mm))
  final.coef = fit.effect$coefficients
  #final.coef[2]
  print("Sample ATT:")
  print(final.coef[2])
  att.sample[itr] = final.coef[2]
  
  ##### Expert var select ####
  expert.var.index = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                 21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,37,38,39,42,43,44,45,46,65)
  expert.var.list = names(XX[, ..expert.var.index])
  treat.form = formula(paste("udpypnr~",paste(expert.var.list, collapse="+")))
  ## performing matching on df_full instead of df
  mm = matchit(treat.form, data = df, method = "nearest", distance = "glm", ratio = 1,
               caliper = .25, link = "logit", estimand = "ATT", replace = FALSE, discard = "control")
  effect.form = formula(paste("suicide_flag ~ udpypnr+",paste(expert.var.list, collapse="+")))
  fit.effect = lm(effect.form, data = match.data(mm))
  final.coef = fit.effect$coefficients
  #final.coef[2]
  print("Estimated true ATT:")
  print(final.coef[2])
  att.expert[itr] = final.coef[2]
  
  
  ########### Enh-ESVM-sigmoid ########
  tic()

  fit1 = svm(udpypnr ~.,
                  data = XA,
                  type = 'C-classification',
                  kernel = 'linear', scale = TRUE)
  # fit1 = glm(udpypnr ~.,family=binomial(link='logit'),data=XA)

  weight_logit = coef(fit1)
  weight_logit = weight_logit[-c(1)] # exclude the intercept

  factor = 1
  penalty = abs(weight_logit)
  penalty[is.na(penalty)] <- 10000000000

  penalty <- vec_sigmoid(penalty)^factor

  # elastic net estimator
  fit2 = cv.glmnet(data.matrix(XX), data.matrix(Y), alpha = 0.5,
                   nfolds = 5L, family = "binomial", gamma = c(0, 0.25, 0.5, 0.75, 1),
                   relax = TRUE, parallel=TRUE, standardize=TRUE,
                   penalty.factor=penalty)

  beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
  beta_el_allvar = as.matrix(beta_el[2:(70)])

  gamma = 1
  penalty_adel = (1/abs(beta_el_allvar))^gamma

  # # penalty_adel <- vec_L2_norm(penalty_adel)
  # vec_L2_norm <- function(vec){
  #   vec <- vec*vec # square element in vec
  #   vec / sqrt(sum(vec^2))
  # }
  # penalty_adel <- vec_L2_norm(penalty_adel)
  penalty_adel[is.na(penalty_adel)] <- 10000000000

  # adaptive elastic net estimator
  fit3 = aenet(data.matrix(XX), data.matrix(Y), family = "binomial", alphas = seq(0.05, 0.95, 0.05),
               nfolds = 5L, rule = "lambda.1se",
               penalty.factor = penalty_adel, parallel = TRUE)

  beta_adnet_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
  beta_adnet_allvar = as.matrix(beta_adnet_raw[1:(69)]) #remove intercept to the last item

  beta_adnet_non_zero = row(beta_adnet_allvar)[which(abs(beta_adnet_allvar) >= 0.1)] # check
  sel.var.index  = beta_adnet_non_zero

  t=toc()
  t_elasped_enhesvmsig = t$toc - t$tic

  print("enh-svm-sigmoid_var_index:")

  print(sel.var.index)
  xx_true = XX[, ..beta_adnet_non_zero]

  enhesvmsig.true.var.list = names(xx_true)
  treat.form = formula(paste("udpypnr~",paste(enhesvmsig.true.var.list ,collapse="+")))
  ## performing matching on df_full instead of df
  mm = matchit(treat.form, data = df, method = "nearest", distance = "glm", ratio = 1,
               caliper = .25, link = "logit", estimand = "ATT", replace = FALSE, discard = "control")

  effect.form = formula(paste("suicide_flag ~ udpypnr+",paste(enhesvmsig.true.var.list ,collapse="+")))
  fit.effect = lm(effect.form, data = match.data(mm))
  final.coef = fit.effect$coefficients
  #final.coef[2]
  enhesvmsig_coef = final.coef[2]
  print("enh-svm-sigmoid estimated ATT:")
  print(enhesvmsig_coef)

  #### store vars
  select.var.list.enh_svm_sig[itr, 1:length(sel.var.index)] = sel.var.index
  total.num.var.enh_svm_sig[itr] = length(sel.var.index)
  att.enh_svm_sig[itr] = final.coef[2]
  time.enh_svm_sig[itr] = t_elasped_enhesvmsig

  ######## Enh-ELogR-tanh #######
  tic()

  fit1 = glm(udpypnr ~.,family=binomial(link='logit'),data=XA)

  weight_logit = coef(fit1)
  weight_logit = weight_logit[-c(1)]

  factor = 0.5
  penalty = abs(weight_logit)
  penalty[is.na(penalty)] <- 10000000000

  penalty <- vec_tanh(penalty)^factor
  # penalty <- vec_sigmoid(penalty)^factor

  # elastic net estimator
  fit2 = cv.glmnet(data.matrix(XX), data.matrix(Y), alpha = 0.5,
                   nfolds = 5L, gamma = c(0, 0.25, 0.5, 0.75, 1),
                   relax = TRUE, parallel=TRUE, standardize=TRUE,
                   family = "binomial", penalty.factor = penalty)

  # K.cv.fit = cv.glmnet(data.matrix(XX), data.matrix(A), type.measure = "mse",
  # nfolds = 5, gamma = c(0, 0.25, 0.5, 0.75, 1), relax = TRUE, family="binomial")

  beta_el = as.matrix(coef(fit2, s=fit2$lambda.1se))
  beta_el_allvar = as.matrix(beta_el[2:70])

  gamma = 1
  penalty_adel = (1/abs(beta_el_allvar))^gamma

  penalty_adel[is.na(penalty_adel)] <- 10000000000

  fit3 = aenet(data.matrix(XX), data.matrix(Y), family = "binomial", alphas = seq(0.05, 0.95, 0.05),
               nfolds = 5L, rule = "lambda.1se", penalty.factor = penalty_adel, parallel = TRUE)

  beta_aden_raw = as.matrix(coef(fit3, s=fit3$lambda.1se))
  beta_aden_allvar = as.matrix(beta_aden_raw[1:69]) #remove intercept to the last item

  beta_aden_non_zero = row(beta_aden_allvar)[which(abs(beta_aden_allvar) >= 0.1)] # check
  sel.var.index = beta_aden_non_zero

  t=toc()
  t_elapsed_enhlogtanh = t$toc - t$tic

  print("enh-eLogR-tanh selected variables")

  print(sel.var.index)
  xx_true = XX[, ..beta_aden_non_zero]

  enhelrtanh.true.var.list = names(xx_true)
  treat.form = formula(paste("udpypnr~",paste(enhelrtanh.true.var.list ,collapse="+")))

  # fit_coef = lm(formula(paste("suicide_flag~udpypnr+",paste(enhelrtanh.true.var.list, collapse="+"))), data = df)
  # ATT_avg_selected_coef = fit_coef$coefficients[2]

  ## performing matching on df_full instead of df
  mm = matchit(treat.form, data = df, method = "nearest", distance = "glm", ratio = 1,
               caliper = .25, link = "logit", estimand = "ATT", replace = FALSE, discard = "control")

  effect.form = formula(paste("suicide_flag ~ udpypnr+",paste(enhelrtanh.true.var.list ,collapse="+")))
  fit.effect = lm(effect.form, data = match.data(mm))
  final.coef = fit.effect$coefficients
  #final.coef[2]
  enhelrtanh_coef = final.coef[2]
  print("enh-eLogR-tanh selected estimated ATT")
  print(enhelrtanh_coef)

  #### store results
  select.var.list.enh_log_tanh[itr, 1:length(sel.var.index)] = sel.var.index
  total.num.var.enh_log_tanh[itr] = length(sel.var.index)
  att.enh_log_tanh[itr] = final.coef[2]
  time.enh_log_tanh[itr] = t_elapsed_enhlogtanh

  ##### Outcome adaptive elastic net #######
  tic()

  my_logit_mod = glm(suicide_flag ~ ., data = df, family = "binomial")
  odds_ratio = (coef(my_logit_mod))
  odds_ratio = odds_ratio[-c(1,43)] # remove odds ratio for intercept and udpypnr  #df <- subset(df, select = -c(a, c))

  # Setting the penalty
  gamma = 5
  penalty = 1/(abs(odds_ratio))^gamma
  penalty[is.na(penalty)] <- 10000000000 # large penalty

  #Perform cross validation
  K.cv.fit = cv.glmnet(data.matrix(XX), data.matrix(A),
                       nfolds = 5, gamma = c(0, 0.25, 0.5, 0.75, 1), relax = TRUE,
                       family="binomial") #check

  cv.alpha = K.cv.fit$relaxed
  alpha.opt = cv.alpha$gamma.1se
  fit2=glmnet(as.matrix(XX), as.matrix(A), alpha =alpha.opt,
              lambda = K.cv.fit$lambda.1se, penalty.factor = penalty , family="binomial")

  beta_ps_raw = as.matrix(coef(fit2))
  beta_ps_raw = as.matrix(coef(fit2, s = 0.02)) # check this
  beta_ps_allvar = as.matrix(beta_ps_raw[2:(70)]) #remove intercept to the last item

  beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.1)] # check
  sel.var.index  = beta_ps_non_zero
  sel.var.index
  print("Onet vars")
  print(sel.var.index)

  t=toc()
  t_elapsed_oaenet = t$toc - t$tic

  xx_true = XX[, ..beta_ps_non_zero]

  onet.true.var.list = names(xx_true)
  treat.form = formula(paste("udpypnr~",paste(onet.true.var.list,collapse="+")))
  ## performing matching on df_full instead of df
  mm = matchit(treat.form, data = df, method = "nearest", distance = "glm",
               ratio = 1, caliper = .25, link = "logit", estimand = "ATT",
               replace = FALSE, discard = "control")

  effect.form = formula(paste("suicide_flag ~ udpypnr+",paste(onet.true.var.list,collapse="+")))
  fit.effect = lm(effect.form, data = match.data(mm))
  final.coef = fit.effect$coefficients
  #final.coef[2]
  oaenetcoef = final.coef[2]
  print("OAENet")
  print(oaenetcoef)

  ### store results
  select.var.list.oaenet[itr, 1:length(sel.var.index)] = sel.var.index
  total.num.var.oaenet[itr] = length(sel.var.index)
  att.oaenet[itr] = final.coef[2]
  time.oaenet[itr] = t_elapsed_oaenet


  ######### Outcome Adaptive LASSO #####
  tic()

  K.cv.fit = cv.glmnet(data.matrix(XX), data.matrix(A), penalty.factor = penalty,
                       nfolds = 5, gamma =c(0, 0.25, 0.5, 0.75, 1),
                       relax = FALSE, family="binomial")
  fit2=glmnet(as.matrix(XX),as.matrix(A), alpha = 1, penalty.factor = penalty, lambda = K.cv.fit$lambda.1se, family="binomial")

  beta_ps_raw = as.matrix(coef(fit2))
  beta_ps_allvar = as.matrix(beta_ps_raw[2:70]) #check

  beta_ps_non_zero = row(beta_ps_allvar)[which(abs(beta_ps_allvar) >= 0.1)]
  #beta_ps_non_zero = row(beta_ps_allvar)[which(!beta_ps_allvar == 0)]
  sel.var.index  = beta_ps_non_zero

  t=toc()
  t_elapsed_oal = t$toc - t$tic

  print("adaptive lasso")
  print(sel.var.index)

  xx_true = XX[,..beta_ps_non_zero]

  oal.true.var.list = names(xx_true)
  treat.form = formula(paste("udpypnr~",paste(oal.true.var.list,collapse="+")))
  #used df_full to match over entire dataset
  mm = matchit(treat.form, data = df, method = "nearest", distance = "glm",
               ratio = 1, link = "logit", estimand = "ATT", replace = FALSE,
               discard = "control") #, caliper = .25

  effect.form = formula(paste("suicide_flag ~ udpypnr+",paste(oal.true.var.list,collapse="+")))
  fit.effect = lm(effect.form, data = match.data(mm))
  final.coef = fit.effect$coefficients
  oalcoef = final.coef[2]
  print("Adap LASSO")
  print(oalcoef)

  ### store results
  select.var.list.oal[itr, 1:length(sel.var.index)] = sel.var.index
  total.num.var.oal[itr] = length(sel.var.index)
  att.oal[itr] = final.coef[2]
  time.oal[itr] = t_elapsed_oal

  # ###### BACR ######
  # tic()
  # bacr = bac(data=((df)), exposure="udpypnr", outcome="suicide_flag", confounders=paste(names(XX), sep = ""),  #changed combined data to DF might have ti change this
  #            interactors=NULL, familyX="binomial", familyY="binomial", omega=Inf,
  #            num_its=200, burnM=1, burnB=1, thin=1)
  # t=toc()
  # time.bacr[itr] = t$toc - t$tic
  # 
  # bacr_df = data.frame(bacr$models)
  # bacr_df = bacr_df[1,1:dim(XX)[2]]
  # sel.var.index = which(bacr_df>0)
  # 
  # print("Bacr selected variables:")
  # print(sel.var.index)
  # 
  # bacr.cov = df[ , ..sel.var.index]
  # 
  # n_bacr = length(sel.var.index)
  # bacr.true.var.list = names(bacr.cov)
  # bacr.data = data.frame(as.matrix(bacr.cov), A, Y)
  # 
  # treat.form = formula(paste("udpypnr~",paste(bacr.true.var.list,collapse="+")))
  # mm = matchit(treat.form, data = df, method = "nearest", distance = "glm",
  #              ratio = 1, link = "logit", estimand = "ATT", replace = FALSE,
  #              discard = "control") #, caliper = .25
  # 
  # effect.form = formula(paste("suicide_flag ~ udpypnr+",paste(bacr.true.var.list,collapse="+")))
  # fit.effect = lm(effect.form, data = match.data(mm))
  # final.coef = fit.effect$coefficients
  # 
  # # df_b = matched_data[,1:(n_bacr + 2)]
  # # fit.effect = lm(Y ~ A + ., data = df_b)
  # # final.coef = fit.effect$coefficients
  # att.bacr[itr] = final.coef[2]
  # 
  # #Storing the Indexes of selected variables
  # select.var.list.bacr[itr, 1:length(sel.var.index)] = sel.var.index
  # 
  # #Storing the total number of variables selected
  # total.num.var.bacr [itr] = n_bacr
  # 
  # #att.BACR[itr] = mean(bacr$ACE)
  # 
  # print("BACR")
  # print(att.bacr[itr])
  # 
  # 
  # 
  # 

  # #########Boruta_T ########
  # tic()
  # bort = Boruta(XX,A,pValue = 0.10)
  # t=toc()
  # time.bort[itr] = t$toc - t$tic
  # #
  # bort_df = data.frame(bort$finalDecision)
  # bort_df$Index = 1:dim.data.frame(bort_df)[1]
  # bort_df = bort_df[bort_df$bort.finalDecision=="Confirmed",]
  # sel.var.index = bort_df$Index
  # bort.cov = XX[ , sel.var.index]
  # #
  # n_bort = length(sel.var.index)
  # bort.cov = data.matrix(bort.cov)
  # bor.data = data.frame(bort.cov, A, Y)
  # bor.data_1 = data.frame(bort.cov, A)
  # mm = matchit(A ~ bort.cov, data = bor.data, method = "nearest", distance = "glm", ratio = 1, caliper = .25, link = "logit", estimand = "ATT", replace = FALSE, discard = "control")
  # #
  # matched_data = (match.data(mm))
  # #
  # df_b = matched_data[,1:(n_bort+2)]
  # fit.effect = lm(Y ~ A + . , data = df_b)
  # final.coef = fit.effect$coefficients
  # att.Boruta_T[itr] =  final.coef[2]
  # #
  # # #Storing the Indexes of selected variables
  # select.var.list.bort[itr, 1:length(sel.var.index)] = sel.var.index
  #
  # #Storing the total number of variables selected
  # total.num.var.bort [itr] = n_bort
  # #
  # print("Boruta T")
  # print(att.Boruta_T[itr])
  #
  #
  # #########Boruta_Y ########
  # tic()
  # bory = Boruta(XX,Y,pValue = 0.10)
  # t=toc()
  # time.bory[itr] = t$toc - t$tic
  #
  # bory_df = data.frame(bory$finalDecision)
  # bory_df$Index = 1:dim.data.frame(bory_df)[1]
  # bory_df = bory_df[bory_df$bory.finalDecision=="Confirmed",]
  # sel.var.index = bory_df$Index
  # bory.cov = XX[ , sel.var.index]
  # #
  # n_bory = length(sel.var.index)
  # bory.cov = data.matrix(bory.cov)
  # bor.y.data = data.frame(bory.cov, A, Y)
  # mm = matchit(A ~ bory.cov, data = bor.y.data, method = "nearest", distance = "glm", ratio = 1, caliper = .25, link = "logit", estimand = "ATT", replace = FALSE, discard = "control")
  #
  # matched_data = (match.data(mm))
  #
  # df_by = matched_data[,1:(n_bory + 2)]
  # fit.effect = lm(Y ~ A + . , data = df_by)
  # final.coef = fit.effect$coefficients
  # att.Boruta_Y[itr] = final.coef[2]
  # #
  # # #Storing the Indexes of selected variables
  # select.var.list.bory[itr, 1:length(sel.var.index)] = sel.var.index
  # #
  # #Storing the total number of variables selected
  # total.num.var.bory[itr] = n_bory
  # #
  # print("Boruta Y")
  # print(att.Boruta_Y[itr])
  #
  # print(i)
  #
  #
  #
  #

  # ######### Deluna(DWR) #######
  # tic()
  #
  # sink("file")
  # fit1 <- try(cov.sel(as.numeric(A), as.numeric(Y), type = "np", XX))
  # sink()
  # if ('try-error' %in% class(fit1)){
  #   att.DWR[itr] = NA
  #   select.var.list.dwr[itr, 1:length(sel.var.index)] = 0
  #   total.num.var.dwr [itr] = NA
  # }
  # else{
  #   dwr_decision = fit1$covar
  #   sel.var.index = which(colnames(XX)==dwr_decision)
  #
  #   t=toc()
  #   time.dwr[itr] = t$toc - t$tic
  #
  #   dwr.cov = Data[ , sel.var.index]
  #   n_dwr<-dim(dwr.cov)[2]
  #   dwr.cov = as.matrix(dwr.cov)
  #   dwr.data = data.frame(dwr.cov, A, Y)
  #   mm <- try(matchit(aa ~ dwr.cov, data = dwr.data, method = "nearest",
  #                     distance = "glm", ratio = 1, caliper = .25, link = "logit",
  #                     estimand = "ATT", replace = FALSE, discard = "both"))
  #
  #   if ('try-error' %in% class(mm)){
  #     att.DWR[itr] = NA
  #     select.var.list.dwr[itr, 1:length(sel.var.index)] = 0
  #     total.num.var.dwr [itr] = NA
  #   }else{
  #     matched_data = (match.data(mm))
  #
  #     df_dwr = matched_data[,1:(n_dwr + 2)]
  #     fit.effect = lm(yy ~ aa + . , data = df_dwr)
  #
  #     ### post record
  #     final.coef = fit.effect$coefficients
  #     # att
  #     att.DWR[itr] = final.coef[2]
  #     #Storing the Indexes of selected variables
  #     select.var.list.dwr[itr, 1:length(sel.var.index)] = sel.var.index
  #     #Storing the total number of variables selected
  #     total.num.var.dwr[itr] = n_dwr
  #     # Percent var select
  #     for (v in sel.var.index){
  #       Var_Select[label, v] <- Var_Select[label, v]+1
  #     }
  #
  #     print("select.var.list.dwr")
  #     print(sel.var.index)
  #     print("att.dwr")
  #     print(att.DWR[itr])
  #   }
  # }
  #
  # ######### CTMLE(Collaborative Targeted Maximum Likelihood Estimation) #######
  # tic()
  # n = length(Y)
  # Q <- cbind(rep(mean(Y[A == 0]), n), rep(mean(Y[A == 1]), n))
  # fit1 <- try(ctmleDiscrete(Y = Y, A = A, W = XX, Q = Q, preOrder = TRUE, detailed = TRUE))
  # if ('try-error' %in% class(fit1)){
  #   att.CTMLE[itr] = NA
  #   select.var.list.ctmle[itr, 1:length(sel.var.index)] = 0
  #   total.num.var.ctmle[itr] = NA
  #   fail_times[label] = fail_times[label]+1
  # }else{
  #   ctmle.cov = fit1$candidates$terms[-c(1)]
  #   sel.var.index = which(colnames(XX)==ctmle.cov)
  # 
  #   t=toc()
  #   time.ctmle[itr] = t$toc - t$tic
  # 
  #   ctmle.cov = Data[ , sel.var.index]
  #   n_ctmle<-dim(ctmle.cov)[2]
  #   ctmle.cov = as.matrix(ctmle.cov)
  #   ctmle.data = data.frame(ctmle.cov, aa, yy)
  #   mm <- try(matchit(aa ~ ctmle.cov, data = ctmle.data, method = "nearest",
  #                     distance = "glm", ratio = 1, caliper = .25, link = "logit",
  #                     estimand = "ATT", replace = FALSE, discard = "both"))
  # 
  #   if ('try-error' %in% class(mm)){
  #     att.CTMLE[itr] = NA
  #     select.var.list.ctmle[itr, 1:length(sel.var.index)] = 0
  #     total.num.var.ctmle[itr] = NA
  #   }else{
  #     matched_data = (match.data(mm))
  # 
  #     if (is.null(n_ctmle)){
  #       att.CTMLE[itr] = NA
  #       select.var.list.ctmle[itr, 1:length(sel.var.index)] = 0
  #       total.num.var.ctmle[itr] = NA
  #     }
  #     else{
  #       df_ctmle <- matched_data[,1:(n_ctmle + 2)]
  #       fit.effect = lm(yy ~ aa + . , data = df_ctmle)
  # 
  #       ### post record
  #       final.coef = fit.effect$coefficients
  #       # att
  #       att.CTMLE[itr] = final.coef[2]
  #       #Storing the Indexes of selected variables
  #       select.var.list.ctmle[itr, 1:length(sel.var.index)] = sel.var.index
  #       #Storing the total number of variables selected
  #       total.num.var.ctmle[itr] = n_ctmle
  #       # Percent var select
  #       for (v in sel.var.index){
  #         Var_Select[label, v] <- Var_Select[label, v]+1
  #       }
  # 
  #       print("select.var.list.ctmle")
  #       print(sel.var.index)
  #       print("att.ctmle")
  #       print(att.CTMLE[itr])
  #     }
  #   }
  # }
  # 
  # 
  # print(t_elapse[,k])

  }

# ####### Writing att data frame ########
# df_att = data.frame(att.sample, att.expert, att.enh_svm_sig, att.enh_log_tanh, att.oaenet, att.oal, att.bacr)
# #df_att = rbind(df_att,ATT)
# write.csv(df_att, paste(workspace, "Processed/NSDUH_att.boot_new.csv", sep=""))
# 
# 
# ###################### Writing att data frame ##########
# df_time = data.frame(time.enh_svm_sig, time.enh_log_tanh, time.oaenet, time.oal, time.bacr)
# #df_att = rbind(df_att,ATT)
# write.csv(df_time, paste(workspace,"Processed/NSDUH_time.boot_new.csv", sep=""))
# 
# 
# ##################### Writing total number of variable selected ##########
# df_sel.var = data.frame(total.num.var.enh_svm_sig, total.num.var.enh_log_tanh, 
#                         total.num.var.oaenet, total.num.var.oal, total.num.var.bacr)
# write.csv(df_sel.var, paste(workspace,"Processed/NSDUH_total.var.boot_new.csv", sep=""))
# 
# 
# #################### Writing the selected variable indexes for each iteration
# write.csv(select.var.list.enh_svm_sig, paste(workspace,
#                                              "Processed/NSDUH_enh_svm_sig.sel.var.list.boot_new.csv"
#                                              ,sep=""))
# write.csv(select.var.list.enh_log_tanh, paste(workspace,
#                                               "Processed/NSDUH_enh_log_tanh.sel.var.list.boot_new.csv"
#                                               , sep=""))
# write.csv(select.var.list.oaenet, paste(workspace,"Processed/NSDUH_oaenet.sel.var.list.boot_new.csv",
#                                         sep=""))
# write.csv(select.var.list.oal, paste(workspace,"Processed/NSDUH_oal.var.list.boot_new.csv", sep=""))
# write.csv(select.var.list.bacr, paste(workspace,"Processed/NSDUH_bacr.var.list.boot_new.csv", sep=""))
# write.csv(select.var.list.bort, paste(workspace,"Processed/NSDUH_bort.var.list.boot_new.csv", sep=""))
# write.csv(select.var.list.bory, paste(workspace,"Processed/NSDUH_bory.var.list.boot_new.csv", sep=""))
# write.csv(select.var.list.dwr, paste(workspace,"Processed/NSDUH_dwr.var.list.boot_new.csv", sep=""))
# write.csv(select.var.list.ctmle, paste(workspace,"Processed/NSDUH_ctmle.var.list.boot_new.csv", sep=""))
# 
# 
# 
# ##############################
# #comb_data = rbind(dat_NA_removed,dat_NA_removed_no_oud)
# write.csv(dat_NA_removed, paste(workspace,"Processed/Final_OUD.csv", sep=""))
# write.csv(dat_NA_removed_no_oud, paste(workspace,"Processed/Final_NO_OUD.csv", sep=""))
# 
# 

# ############## BCEE ##########################
# xxx = as.matrix(unname(as.matrix(XX)))
# tic()
# bcee = GBCEE((as.matrix(A)), as.matrix(Y), xxx, omega = 300*sqrt(6659), niter = 5000, family.X = "binomial", family.Y = "binomial",
#              X1 = 1, X0 = 0, priorX = NA, priorY = NA, maxsize = NA, OR = 20, truncation = c(0.01, 0.99), var.comp = "asymptotic", B = 200)
# t=toc()
# time.bcee[itr] = t$toc - t$tic
# 
# #Considering the Posterior distribution of the outcome model
# bcee_df = data.frame(bcee$models.Y)
# 
# #Selecting the top 5 models as they includes 80% of the posterior probability
# freq_var_sel = colSums(bcee_df[1:5,])
# 
# #Select the variables who appear in at least 50% of the models
# sel.var.index = which((freq_var_sel[1:100] > 2))
# bcee.cov = df[, ..sel.var.index]
# 
# n_bcee = length(sel.var.index)
# bcee.cov = as.matrix(bcee.cov)
# df_bcee = data.frame(bcee.cov, A, Y)
# mm = matchit(A ~ bcee.cov, data = df_bcee, method = "nearest", distance = "glm", ratio = 1, caliper = .25, link = "logit", estimand = "ATT", replace = FALSE, discard = "both")
# 
# matched_data = (match.data(mm))
# 
# df_b = matched_data[,1:(n_bcee + 2)]
# fit.effect = lm(Y ~ A + ., data = df_b)
# final.coef = fit.effect$coefficients
# att.BCEE[itr] = final.coef[2]
# 
# #Storing the Indexes of selected variables
# select.var.list.bcee[itr, 1:length(sel.var.index)] = sel.var.index
# 
# #Storing the total number of variables selected
# total.num.var.bcee [itr] = n_bcee
# 
# print("BCEE")
# print(att.BCEE[itr])

##########################################



# # one average how many variables are selected?
# dummy=0
# control_dummy = 0
# minimim = 1
# maximum = 0
# difference = 0
# for(i in 1:1000){
#   control_dummy = dummy
#   for(j in 1:70){
#     if (select.var.list.oal[i,j]==0){
#       next
#     }
#     dummy=dummy+1
#   }
#   difference = dummy - control_dummy
#   if(difference>maximum){
#     maximum = difference
#   }
#   if(difference<minimim){
#     minimim = difference
#   }
# }
# 
# 
# 
# 
# 
# 
# ###########
# count = zeros(70)
# #find if the variables are selected at least 80% of the times
# for (k in 1:70){
#   for(i in 1:70){
#     for (j in 1:1000){
#       if(select.var.list.oal[j,i]==k){
#         count[k]=count[k]+1
#       }
#     }
#   }
# 
# }
# 
# count = count/1000
# 
# 
# ###########
# #which method on average selected similar variables
# # you could do something like, if the percentage difference is within 10 points then the methods are similar
# 
# 
# ##########
# #select the variables instead of numbers
# 
# 
# 
# #########
# #some kind of visualization for display the selected variables
# 

