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

# Load data ---------------------------------------------------------------
data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
output_data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/'
net_data_dir='/cbica/projects/spatial_topography/data/imageData/net_stats/'
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

# Data cleaning -----------------------------------------------------------

#Only baseline visits
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
main_schaeferwsbm <- left_join(net_stats_schaeferwsbm, main, by="ID")

#recode income
main$demo_comb_income_numeric <- recode(as.numeric(as.character(main$demo_comb_income_v2)), "1 = 2500; 2 = 8500; 3 = 14000; 4 = 20500; 5 = 30000; 6 = 42500; 7 = 62500; 8 = 87500; 9 = 150000; 10 = 200000; 999 = NA ; 777 = NA")  

# WM - List Sorting -------------------------------------------------------
#Schaefer400-Yeo7
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7")
for (net in nets){
  name<-paste0("lm_wm_ls_",net)
  #formula<-formula(paste0(net, '~age_scan*ses_composite_med+male+fd_mean+avgweight+pctSpikesFD+size_t'))
  #need to get mean fd and avgweight?
  formula<-formula(paste0('nihtbx_list_uncorrected~age+gender+', net))
  assign(name, lm(formula, data=main_schaeferyeo7))
  print(summary(get(name)))
  print(lm.beta(get(name)))
  name<-paste0("gam_wm_ls_",net)
  gamformula<-formula(paste0('nihtbx_list_uncorrected~age+gender+s(', net,')'))
  assign(name, gam(gamformula, data=main_schaeferyeo7, REML=TRUE))
  print(exactRLRT(gamm(gamformula, data=main_schaeferyeo7, REML=TRUE)$lme))
}
#
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys4to7), data=main_schaeferyeo7, REML=TRUE)$lme)
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys6to7), data=main_schaeferyeo7, REML=TRUE)$lme)
visreg(gam_wm_ls_sys6to7)#this does not look believably non-linear

#FPN-FPN predicts WM, with or without avgweight, but AIC and BIC are lower without avg weight
#Model selection
BIC(lm_wm_ls_sys6to6)
AIC(lm_wm_ls_sys6to6)

#Schaefer400-WSBM
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7")
for (net in nets){
  name<-paste0("lm_wm_ls_",net)
  #formula<-formula(paste0(net, '~age_scan*ses_composite_med+male+fd_mean+avgweight+pctSpikesFD+size_t'))
  #need to get mean fd and avgweight?
  formula<-formula(paste0('pnihtbx_list_uncorrected~age+gender+', net))
  assign(name, lm(formula, data=main_schaeferwsbm))
  print(summary(get(name)))
  print(lm.beta(get(name)))
  name<-paste0("gam_wm_ls_",net)
  gamformula<-formula(paste0('nihtbx_list_uncorrected~age+gender+s(', net,')'))
  assign(name, gam(gamformula, data=main_schaeferwsbm, REML=TRUE))
  print(exactRLRT(gamm(gamformula, data=main_schaeferwsbm, REML=TRUE)$lme))
}
#
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys4to7), data=main_schaeferwsbm, REML=TRUE)$lme)
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys6to7), data=main_schaeferwsbm, REML=TRUE)$lme)
visreg(gam_wm_ls_sys4to7)#this does not look believably non-linear

#No prediction

#Model selection, worse than Schaefer400-Yeo7 for FPN-FPN
BIC(lm_wm_ls_sys6to6)
AIC(lm_wm_ls_sys6to6)

#From fsaverage6-Yeo dev
WAITING!
  
# TEST SAMPLE -------------------------------------------------------------

#for site 20 and 14
#subjlist2 <- read.table(paste0(subjlist_dir, "site14site20/n611_release2_site14site20_0.2mm.txt"), col.names = c("subjectid"))
#subjlist2$subjectid <- str_remove(subjlist2$subjectid, "sub-")
#wholesubjectlist <- rbind(subjlist1, subjlist2)