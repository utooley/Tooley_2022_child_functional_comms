library(dplyr)
library(summarytools)
library(PerformanceAnalytics)
library(stringr)
library(car)
library(R.matlab)
library(mgcv)
library(RLRsim)
library(lm.beta)
library(visreg)
library(ggplot2)
library(data.table)
library(vioplot)
library(performance)

# Load data ---------------------------------------------------------------
data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
output_data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/'
net_data_dir='/cbica/projects/spatial_topography/data/imageData/net_stats/'
raw_data_dir='/data/jux/mackey_group/public_data/ABCD/'
sites <- read.delim(paste0(data_dir,"abcd_lt01.txt"), stringsAsFactors = F)# change strings as factors = FALSE next time
sites = sites[-1,]
sites = droplevels(sites)
sites <- select(sites, -c(collection_id:dataset_id))
sites$ID <- sites$subjectkey
socio <- read.delim(paste0(data_dir,"acspsw03.txt"), stringsAsFactors = F)
socio = socio[-1,]
socio = droplevels(socio)
socio <- select(socio, -c(collection_id:dataset_id)) 
socio$ID <- socio$subjectkey
income <- read.delim(paste0(data_dir,"pdem02.txt")) 
income = income[-1,]
income = droplevels(income)
income <- select(income, -c(collection_id:dataset_id))
income$ID <- income$subjectkey
#get NIH Toolbox cognition data
nihtbx <- read.delim(paste0(data_dir,"abcd_tbss01.txt"), stringsAsFactors = F) 
nihtbx = nihtbx[-1,]
nihtbx = droplevels(nihtbx)
nihtbx$ID <- nihtbx$subjectkey
nihtbx <- select(nihtbx, -c(collection_id:dataset_id)) 
nihtbx <- nihtbx %>% mutate_at(vars(matches("uncorrected|agecorrected|theta|itmcnt|fc|cs")), as.numeric)
#get Pearson scores from WISC-V cognition data
wiscv <- read.delim(paste0(data_dir,"abcd_ps01.txt"), stringsAsFactors = F) 
wiscv = wiscv[-1,]
wiscv = droplevels(wiscv)
wiscv$ID <- wiscv$subjectkey
wiscv <- select(wiscv, -c(collection_id:dataset_id)) 
wiscv <- wiscv %>% mutate_at(vars(matches("ravlt|wiscv")), as.numeric)
#get task fMRI task scores
taskfmri <- read.delim(paste0(data_dir,"abcd_mrinback02.txt"), stringsAsFactors = F) 
taskfmri = taskfmri[-1,]
taskfmri = droplevels(taskfmri)
taskfmri$ID <- taskfmri$subjectkey
taskfmri <- select(taskfmri, -c(collection_id:dataset_id)) 
taskfmri <- taskfmri %>% mutate_at(vars(tfmri_nback_beh_switchflag:tfmri_nb_r2_beh_c0bpnl_stdrt), as.numeric)
#get recognition memory post-scanner task scores
recmem <- read.delim(paste0(data_dir,"mribrec02.txt"), stringsAsFactors = F) 
recmem = recmem[-1,]
recmem = droplevels(recmem)
recmem$ID <- recmem$subjectkey
recmem <- select(recmem, -c(collection_id:dataset_id)) 
recmem <- recmem %>% mutate_at(vars(tfmri_rec_all_beh_newnf_hr:tfmri_rec_all_beh_place_dp), as.numeric)

main<- left_join(sites,socio, by=c("ID", "eventname"))
main <- left_join(main, income, by=c("ID", "eventname"))
main <- left_join(main, nihtbx, by=c("ID", "eventname"))
main <- left_join(main, wiscv, by=c("ID", "eventname"))
main <- left_join(main, taskfmri, by=c("ID", "eventname"))

#####################################
########### TRAINING SAMPLE #########
#####################################
#Get network connectivity data
#From Schaefer400-Yeo7
net_stats_schaeferyeo7 <- read.csv(paste0(net_data_dir,"n670_training_sample_schaefer400_yeo7_network_stats.csv"))
#From Schaefer400-WSBM
net_stats_schaeferwsbm <- read.csv(paste0(net_data_dir,"n670_site16_training_sample_schaefer400_wsbm_network_stats.csv"))
#From fsaverage6-Yeo dev

##need to get XCP mean FD and # of outliers and control for that
runs <- read.csv(paste0(subjlist_dir,'parcellation/n670_filtered_runs_site16_postprocess.csv'))
runs <- select(runs, id, var1:var2) #take only the first 2 runs that were used
runs<- melt(runs, measure=c("var1", "var2")) %>% arrange(., id) %>% select(., -variable) %>% rename(.,ID=id,run=value)#reshape them and take out extra
#read in the list of runs that were used and merge it with xcp qa vars
qa_vars <- read.csv(paste0(raw_data_dir, "bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/XCP_QAVARS_with_nVols.csv"))
qa_vars <- rename(qa_vars, ID=id0, run=id1)
qa <- left_join(runs, qa_vars, by=c("ID", "run"))

# Data cleaning -----------------------------------------------------------
#summarize motion and merge in
#make summary variables of volumes and censored volumes
qa <- qa %>% group_by(ID) %>% mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(nVols)) %>% ungroup()

qa <- qa %>% mutate(perc_vols=nVols/totalSizet, relMeanRMSMotion_weight=relMeanRMSMotion*perc_vols, pctSpikesFD_weight=pctSpikesFD*perc_vols) %>% group_by(ID) %>% 
  mutate(fd_mean_avg=sum(relMeanRMSMotion_weight), pctVolsCensored=(totalnVolCensored/totalSizet), pctSpikesFD_avg=sum(pctSpikesFD_weight)) %>% select(fd_mean_avg, pctVolsCensored, pctSpikesFD_avg, totalSizet)
#average motion across the two runs, weighted by the length of each run as a percentage of the total, same for percent spikes FD. select only the averaged variables
qa <- qa[!duplicated(qa$ID), ]

#Only baseline visits from cognition
main <- filter(main, eventname=="baseline_year_1_arm_1") %>% filter(.,site_id_l=="site16") #filter out only the baseline visits
#for site 16 only
#subjlist1 <- read.table(paste0(subjlist_dir, "/n670_subjects_only_filtered_runs_site16_postprocess.txt"), col.names = c("subjectid"))
#Keep only n670 subject list
main$ID <- str_remove(main$ID, "_") #subjectid has a _ in it after NDAR in the NDAR data, doesn't match subject list
main$ID <- paste0("sub-",main$ID)
#make age numeric
main$age <- as.numeric(main$interview_age.x)/12
#make sex a factor
main$gender <- as.factor(main$gender.x)

main_schaeferyeo7 <- left_join(net_stats_schaeferyeo7, main, by="ID")
main_schaeferyeo7 <- left_join(main_schaeferyeo7, qa, by="ID")
main_schaeferwsbm <- left_join(net_stats_schaeferwsbm, main, by="ID")
main_schaeferwsbm <- left_join(main_schaeferwsbm, qa, by="ID")

#recode income
main$demo_comb_income_numeric <- recode(as.numeric(as.character(main$demo_comb_income_v2)), "1 = 2500; 2 = 8500; 3 = 14000; 4 = 20500; 5 = 30000; 6 = 42500; 7 = 62500; 8 = 87500; 9 = 150000; 10 = 200000; 999 = NA ; 777 = NA")  

# Plot data descriptives --------------------------------------------------
measures=select(main_schaeferyeo7,one_of("nihtbx_list_uncorrected","pea_wiscv_tss", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate"))
par(mfrow=c(4,3))
i=1
for (meas in measures){
  name=colnames(measures)[i]
  hist(meas, main=name, col = "lightblue") #hist of measure
  i=i+1
  scatter.smooth(main_schaeferyeo7$age,meas,  col = "blue", main="versus age")
  p <- vioplot(meas~main_schaeferyeo7$gender, main="versus gender") #gender
}
#take out outliers!

# Cognitive measures-------------------------------------------------------
#add checks for model fit and performance with check_model() and model_performance()

measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate")
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7")
#Schaefer400-Yeo7
for (meas in measures){
  print(meas)
for (net in nets){
  print(net)
  name<-paste0("lm_", meas,"_",net)
  formula<-formula(paste0(net,'~age+gender+', meas))
  assign(name, lm(formula, data=main_schaeferyeo7))
  print(summary(get(name)))
  print(lm.beta(get(name)))
  name<-paste0("gam_wm_ls_",net)
  gamformula<-formula(paste0(meas,'~s(age)+gender+fd_mean_avg+avgweight+', net))
  assign(name, gam(gamformula, data=main_schaeferyeo7, REML=TRUE))
  print(exactRLRT(gamm(gamformula, data=main_schaeferyeo7, REML=TRUE)$lme))
}
}
#
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys4to7), data=main_schaeferyeo7, REML=TRUE)$lme)
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys6to7), data=main_schaeferyeo7, REML=TRUE)$lme)
visreg(lm_nihtbx_list_uncorrected_sys6to6)#this does not look believably non-linear

#FPN-FPN predicts WM, with or without avgweight, but AIC and BIC are lower without avg weight
#Model selection
BIC(lm_wm_ls_sys6to6)
AIC(lm_wm_ls_sys6to6)

#Schaefer400-WSBM
measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate")
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7")
for (meas in measures){
print(meas)
for (net in nets){
  print(net)
  name<-paste0("lm_", meas,"_",net)
  formula<-formula(paste0(net,'~age+gender+', meas))
  assign(name, lm(formula, data=main_schaeferwsbm))
  print(summary(get(name)))
  print(lm.beta(get(name)))
  name<-paste0("gam_wm_ls_",net)
  gamformula<-formula(paste0(meas,'~s(age)+gender+fd_mean_avg+avgweight+', net))
  assign(name, gam(gamformula, data=main_schaeferwsbm, REML=TRUE))
  print(exactRLRT(gamm(gamformula, data=main_schaeferwsbm, REML=TRUE)$lme))
}
}
#
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys4to7), data=main_schaeferwsbm, REML=TRUE)$lme)
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys6to7), data=main_schaeferwsbm, REML=TRUE)$lme)
visreg(lm_nihtbx_list_uncorrected_sys6to6)#this does not look believably non-linear

#No prediction

#Model selection, worse than Schaefer400-Yeo7 for FPN-FPN
BIC(lm_wm_ls_sys6to6)
AIC(lm_wm_ls_sys6to6)

#From fsaverage6-Yeo dev
WAITING!
  
  
#####################################
########### TEST SAMPLE #########
#####################################
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/'

#MAKE SURE TO TAKE OUT THE TWO SUBJECTS WHO NEED TO BE EXCLUDED.

#Get network connectivity data
#From Schaefer400-Yeo7
net_stats_schaeferyeo7 <- read.csv(paste0(net_data_dir,"n544_test_sample_schaefer400_yeo7_network_stats.csv"))
#From Schaefer400-WSBM
net_stats_schaeferwsbm <- read.csv(paste0(net_data_dir,"n544_site14site20_test_sample_schaefer400_wsbm_network_stats.csv"))
#From fsaverage6-Yeo dev

#remake main
main<- left_join(sites,socio, by=c("ID", "eventname"))
main <- left_join(main, income, by=c("ID", "eventname"))
main <- left_join(main, nihtbx, by=c("ID", "eventname"))
main <- left_join(main, wiscv, by=c("ID", "eventname"))
main <- left_join(main, taskfmri, by=c("ID", "eventname"))

##need to get XCP mean FD and # of outliers and control for that
runs <- read.csv(paste0(subjlist_dir,'n544_filtered_runs_site14site20_postprocess.csv'))
runs <- select(runs, id, var1:var2) #take only the first 2 runs that were used
runs<- melt(runs, measure=c("var1", "var2")) %>% arrange(., id) %>% select(., -variable) %>% rename(.,ID=id,run=value)#reshape them and take out extra
#read in the list of runs that were used and merge it with xcp qa vars
qa_vars <- read.csv(paste0(raw_data_dir, "bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/XCP_QAVARS_with_nVols.csv"))
qa_vars <- rename(qa_vars, ID=id0, run=id1)
qa <- left_join(runs, qa_vars, by=c("ID", "run"))

# Data cleaning -----------------------------------------------------------
#summarize motion and merge in
#make summary variables of volumes and censored volumes
qa <- qa %>% group_by(ID) %>% mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(nVols)) %>% ungroup()

qa <- qa %>% mutate(perc_vols=nVols/totalSizet, relMeanRMSMotion_weight=relMeanRMSMotion*perc_vols, pctSpikesFD_weight=pctSpikesFD*perc_vols) %>% group_by(ID) %>% 
  mutate(fd_mean_avg=sum(relMeanRMSMotion_weight), pctVolsCensored=(totalnVolCensored/totalSizet), pctSpikesFD_avg=sum(pctSpikesFD_weight)) %>% select(fd_mean_avg, pctVolsCensored, pctSpikesFD_avg, totalSizet)
#average motion across the two runs, weighted by the length of each run as a percentage of the total, same for percent spikes FD. select only the averaged variables
qa <- qa[!duplicated(qa$ID), ]

#Only baseline visits from sites 14 and 20
main <- filter(main, eventname=="baseline_year_1_arm_1") %>% filter(.,site_id_l=="site20"|site_id_l=="site14") #filter out only the baseline visits
main$ID <- str_remove(main$ID, "_") #subjectid has a _ in it after NDAR in the NDAR data, doesn't match subject list
main$ID <- paste0("sub-",main$ID)
#make age numeric
main$age <- as.numeric(main$interview_age.x)/12
#make sex a factor
main$gender <- as.factor(main$gender.x)

main_schaeferyeo7 <- left_join(net_stats_schaeferyeo7, main, by="ID")
main_schaeferyeo7 <- left_join(main_schaeferyeo7, qa, by="ID")
main_schaeferwsbm <- left_join(net_stats_schaeferwsbm, main, by="ID")
main_schaeferwsbm <- left_join(main_schaeferwsbm, qa, by="ID")

# Plot data descriptives --------------------------------------------------
measures=select(main_schaeferyeo7,one_of("nihtbx_list_uncorrected","pea_wiscv_tss", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate"))
par(mfrow=c(4,3))
i=1
for (meas in measures){
  name=colnames(measures)[i]
  hist(meas, main=name, col = "lightblue") #hist of measure
  i=i+1
  scatter.smooth(main_schaeferyeo7$age,meas,  col = "blue", main="versus age")
  p <- vioplot(meas~main_schaeferyeo7$gender, main="versus gender") #gender
}

# Cognitive measures-------------------------------------------------------
measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate")
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7")
#Schaefer400-Yeo7
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(net,'~age+gender+fd_mean_avg+avgweight+', meas))
    assign(name, lm(formula, data=main_schaeferyeo7))
    print(summary(get(name)))
    print(lm.beta(get(name)))
    name<-paste0("gam_wm_ls_",net)
    gamformula<-formula(paste0(meas,'~s(age)+gender+fd_mean_avg+avgweight+', net))
    assign(name, gam(gamformula, data=main_schaeferyeo7, REML=TRUE))
    print(exactRLRT(gamm(gamformula, data=main_schaeferyeo7, REML=TRUE)$lme))
  }
}

#Schaefer400-WSBM
measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate")
for (meas in measures){
  nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7")
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(net,'~age+gender+fd_mean_avg+avgweight+', meas))
    assign(name, lm(formula, data=main_schaeferwsbm))
    print(summary(get(name)))
    print(lm.beta(get(name)))
    name<-paste0("gam_wm_ls_",net)
    gamformula<-formula(paste0(meas,'~s(age)+gender+fd_mean_avg+avgweight+', net))
    assign(name, gam(gamformula, data=main_schaeferwsbm, REML=TRUE))
    print(exactRLRT(gamm(gamformula, data=main_schaeferwsbm, REML=TRUE)$lme))
  }
}

