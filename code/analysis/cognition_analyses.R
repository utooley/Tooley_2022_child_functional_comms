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
library(parallel)
library(PerformanceAnalytics)

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
main <- left_join(main, recmem, by=c("ID", "eventname"))

#####################################
########### TRAINING SAMPLE #########
#####################################
#Get network connectivity data
#From Schaefer400-Yeo7
net_stats_schaeferyeo7 <- read.csv(paste0(net_data_dir,"n670_training_sample_schaefer400_yeo7_network_stats.csv"))
#From Schaefer400-WSBM
net_stats_schaeferwsbm <- read.csv(paste0(net_data_dir,"n670_site16_training_sample_schaefer400_wsbm_network_stats.csv"))
#From fsaverage6-Yeo dev
net_stats_yeodev <- read.csv(paste0(net_data_dir,"n670_training_sample_fsaverage6_yeodev_network_stats.csv"))

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
#average nback face measures
main$nback_2back_place_rate_correct <- main$tfmri_nb_all_beh_c2bp_rate #is rate of correct for 2-back places
main$nback_2back_face_rate_correct <- (main$tfmri_nb_all_beh_c2bpf_rate+main$tfmri_nb_all_beh_c2bnf_rate+main$tfmri_nb_all_beh_c2bngf_rate)/3 #average together all the face rates
#average dprime measures
main$dprime_faces <- (main$tfmri_rec_all_beh_posf_dpr+main$tfmri_rec_all_beh_neutf_dp+main$tfmri_rec_all_beh_negf_dp)/3
main$dprime_places <- main$tfmri_rec_all_beh_place_dp

main_schaeferyeo7 <- left_join(net_stats_schaeferyeo7, main, by="ID")
main_schaeferyeo7 <- left_join(main_schaeferyeo7, qa, by="ID")
main_schaeferwsbm <- left_join(net_stats_schaeferwsbm, main, by="ID")
main_schaeferwsbm <- left_join(main_schaeferwsbm, qa, by="ID")
main_yeodev <- left_join(net_stats_yeodev, main, by="ID")
main_yeodev <- left_join(main_yeodev, qa, by="ID")

#average participation coefficient
net_stats_schaeferyeo7$partcoef <- (net_stats_schaeferyeo7$sub_partcoef_neg_yeo+net_stats_schaeferyeo7$sub_partcoef_pos_yeo)/2
net_stats_schaeferwsbm$partcoef <- (net_stats_schaeferwsbm$sub_partcoef_neg_wsbm+net_stats_schaeferwsbm$sub_partcoef_pos_wsbm)/2
#net_stats_yeodev$partcoef <- (net_stats_yeodev$+net_stats_schaeferyeo7$sub_partcoef_pos_yeo)/2 DIDN'T DO IN YEO-DEV

#take out the one outlier in within- and between- connectivity
main_schaeferyeo7 <- filter(main_schaeferyeo7, ID != "sub-NDARINV4H7G4RXD")
main_schaeferwsbm <- filter(main_schaeferwsbm, ID != "sub-NDARINV4H7G4RXD")
main_yeodev <- filter(main_yeodev, ID != "sub-NDARINV4H7G4RXD")

#recode income if I end up using it
main$demo_comb_income_numeric <- dplyr::recode(as.numeric(as.character(main$demo_comb_income_v2)), "1 = 2500; 2 = 8500; 3 = 14000; 4 = 20500; 5 = 30000; 6 = 42500; 7 = 62500; 8 = 87500; 9 = 150000; 10 = 200000; 999 = NA ; 777 = NA")  

# Plot data descriptives --------------------------------------------------
measures=select(main_schaeferyeo7,one_of("nihtbx_list_uncorrected","pea_wiscv_trs", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate","dprime_faces", "dprime_places"))
par(mfrow=c(5,3))
i=1
for (meas in measures){
  name=colnames(measures)[i]
  hist(meas, main=name, col = "lightblue") #hist of measure
  i=i+1
  scatter.smooth(main_schaeferyeo7$age,meas,  col = "blue", main="versus age")
  p <- vioplot(meas~main_schaeferyeo7$gender, main="versus gender") #gender
}
#plot each of them with each of the others
chart.Correlation(measures)

#take out those subjects who were flagged for poor performance on the nback, about 5% of them
main_schaeferyeo7_nback <- filter(main_schaeferyeo7, tfmri_nback_beh_performflag==1 | is.na(tfmri_nback_beh_performflag))# %>% filter(., nback_2back_place_rate_correct > 0.6)
main_schaeferwsbm_nback <- filter(main_schaeferwsbm, tfmri_nback_beh_performflag==1 | is.na(tfmri_nback_beh_performflag))# %>% filter(., nback_2back_place_rate_correct > 0.6)
main_yeodev_nback <- filter(main_yeodev, tfmri_nback_beh_performflag==1 | is.na(tfmri_nback_beh_performflag)) %>% filter(., nback_2back_place_rate_correct > 0.6)

# Cognitive measures-------------------------------------------------------
#add checks for model fit and performance with check_model() and model_performance()

measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "dprime_places", "dprime_faces")
nets=c("system_segreg_yeo","mean_within_sys_yeo", "mean_between_sys_yeo","modul_yeo")
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7", "sys5to7", "sys5to6")
measures=c("pea_wiscv_trs")
#Schaefer400-Yeo7
for (meas in measures){
  print(meas)
for (net in nets){
  print(net)
  name<-paste0("lm_", meas,"_",net)
  formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+', net))
  assign(name, lm(formula, data=main_schaeferyeo7))
  print(summary(get(name)))
  print(lm.beta(get(name)))
}
}

par(mfrow=c(3,4))
for (net in nets){
  #visreg(get(paste0("lm_dprime_faces_",net)), net)
  avPlot(get(paste0("lm_dprime_faces_",net)), net)
  #plot_model(get(paste0("lm_dprime_faces_",net)), terms = net, type="eff")
}
  
measures=c("tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate", "nback_2back_place_rate_correct","nback_2back_face_rate_correct")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(net,'~age+gender+fd_mean_avg+avgweight+', meas))
    assign(name, lm(formula, data=main_schaeferyeo7_nback))
    print(summary(get(name)))
  print(lm.beta(get(name)))
  #name<-paste0("gam_wm_ls_",net)
  #gamformula<-formula(paste0(meas,'~s(age)+gender+fd_mean_avg+avgweight+', net))
  #assign(name, gam(gamformula, data=main_schaeferyeo7_nback, REML=TRUE))
  #print(exactRLRT(gamm(gamformula, data=main_schaeferyeo7_nback, REML=TRUE)$lme))
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
measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "dprime_places", "dprime_faces")
nets=c("system_segreg_wsbm","mean_within_sys_wsbm", "mean_between_sys_wsbm","modul_wsbm")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+', net))
    assign(name, lm(formula, data=main_schaeferwsbm))
    print(summary(get(name)))
    print(lm.beta(get(name)))
  }
}

for (net in nets){
  crPlot(get(paste0("lm_dprime_faces_",net)), net)
}

measures=c("tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate", "nback_2back_place_rate_correct","nback_2back_face_rate_correct")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(net,'~age+gender+fd_mean_avg+avgweight+', meas))
    assign(name, lm(formula, data=main_schaeferwsbm_nback))
    print(summary(get(name)))
    print(lm.beta(get(name)))
    #name<-paste0("gam_wm_ls_",net)
    #gamformula<-formula(paste0(meas,'~s(age)+gender+fd_mean_avg+avgweight+', net))
    #assign(name, gam(gamformula, data=main_schaeferyeo7_nback, REML=TRUE))
    #print(exactRLRT(gamm(gamformula, data=main_schaeferyeo7_nback, REML=TRUE)$lme))
  }
}


for (net in nets){
  visreg(get(paste0("lm_dprime_faces_",net)), net)
}

exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys4to7), data=main_schaeferwsbm, REML=TRUE)$lme)
exactRLRT(gamm(nihtbx_list_uncorrected~age+gender+s(sys6to7), data=main_schaeferwsbm, REML=TRUE)$lme)
visreg(lm_nback_2back_place_rate_correct_sys6to6)#this does not look believably non-linear

#No prediction

#Model selection, worse than Schaefer400-Yeo7 for FPN-FPN
BIC(lm_wm_ls_sys6to6)
AIC(lm_wm_ls_sys6to6)

#From fsaverage6-Yeo dev
#measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "dprime_places", "dprime_faces")
nets=c("system_segreg_yeodev", "mean_within_sys_yeodev","mean_between_sys_yeodev")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+', net))
    assign(name, lm(formula, data=main_yeodev))
    print(summary(get(name)))
    print(lm.beta(get(name)))
  }
}

for (net in nets){
  visreg(get(paste0("lm_dprime_faces_",net)), net)
}

measures=c("tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate", "nback_2back_place_rate_correct","nback_2back_face_rate_correct")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(net,'~age+gender+fd_mean_avg+avgweight+', meas))
    assign(name, lm(formula, data=main_yeodev_nback))
    print(summary(get(name)))
    print(lm.beta(get(name)))
  }
}
visreg(lm_nback_2back_face_rate_correct_sys4to4)


# Age effects-descriptive -------------------------------------------------------------
yeo_measures=select(main_schaeferyeo7,one_of("age", "system_segreg_yeo","mean_within_sys_yeo", "mean_between_sys_yeo","modul_yeo"))
wsbm_measures=select(main_schaeferwsbm,one_of("age", "system_segreg_wsbm","mean_within_sys_wsbm", "mean_between_sys_wsbm","modul_wsbm"))
yeodev_measures=select(main_yeodev,one_of("age", "system_segreg_yeodev","mean_within_sys_yeodev", "mean_between_sys_yeodev"))

par(mfrow=c(5,3))
i=1
for (meas in yeo_measures){
  name=colnames(measures)[i]
  hist(meas, main=name, col = "lightblue") #hist of measure
  i=i+1
  scatter.smooth(main_schaeferyeo7$age,meas,  col = "blue", main="versus age")
  p <- vioplot(meas~main_schaeferyeo7$gender, main="versus gender") #gender
}
par(mfrow=c(3,5))
for (meas in names(yeodev_measures)){
  name=meas
  print(meas)
  name<-paste0("lm_", meas)
  formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight'))
  assign(name, lm(formula, data=main_yeodev))
  visreg(get(name), "age", main=meas)
  print(summary(get(name)))
}

#plot each of them with each of the others
chart.Correlation(yeo_measures)
chart.Correlation(wsbm_measures)
chart.Correlation(yeodev_measures)

# Age and sex effects-comparison --------------------------------------------------
# Make dataframes for each partitions with the strength of age effects and pvals!
covariates="~age+gender+fd_mean_avg+avgweight" #to look at age
covariates="~gender+age+fd_mean_avg+avgweight" #to look at sex

partitions=c("schaeferyeo7","schaeferwsbm", "yeodev")
ends=c("yeo", "wsbm", "yeodev")
i=1
for (partition in partitions){
  end=ends[i]
main_unique <- dplyr::select(get(paste0("main_",partition)), -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
#run and compare for multiple comparisons again
cols <- names((dplyr::select(ungroup(main_unique),matches("system_segreg| mean_within_sys_ | mean_between_sys_|sys"))))
m <- mclapply(cols, function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
#betas
networks_Age_betas <- mclapply(m, function(sys) {summary(lm(formula = sys,data=main_unique))$coefficients[2,1]},mc.cores=1)
networks_Age_betas <- as.data.frame(networks_Age_betas)
networks_Age_betas <- t(networks_Age_betas)
networks_Age_betas <- as.data.frame(networks_Age_betas)
#pvals
networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[2,4]},mc.cores=1)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
networks_Age_pvals <- t(networks_Age_pvals)
networks_Age_pvals <- as.data.frame(networks_Age_pvals)
#bonferroni correct
networks_Age_pvals_fdr <- p.adjust(networks_Age_pvals$V1, method="fdr")
networks_age_fdr <- data.frame(networks_Age_betas[,1],networks_Age_pvals,c("system_segreg","mean_within_sys","mean_between_sys",names(dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
colnames(networks_age_fdr) <- c("beta","pvalue", "network")
networks_age_fdr$partition <- rep(end, length(networks_age_fdr$beta))
#FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
print(networks_age_fdr)
assign(paste0("networks_age_fdr_",end), networks_age_fdr)
i=i+1
}

#Plot the comparison
full <- rbind(networks_age_fdr_yeo,networks_age_fdr_wsbm, networks_age_fdr_yeodev)
dim(full)
colnames(full)

#reorder them to have the right order
full$partition <- factor(full$partition,levels = ends)
full$network <- factor(full$network, levels=c("system_segreg","mean_within_sys","mean_between_sys",names(dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
ggplot(aes(x = partition, y = network, fill = beta), data=full)+geom_tile()+scale_fill_gradient2(low = ("lightblue"),
                                                                                                 mid="white",high = "red", midpoint=0, limits=c(-0.021, 0.021))+theme_ipsum()+geom_text(aes(label=ifelse(pvalue<0.05, "*", ""))) 
# Sex effects -------------------------------------------------------------
yeo_measures=select(main_schaeferyeo7,one_of("age", "system_segreg_yeo","mean_within_sys_yeo", "mean_between_sys_yeo","modul_yeo"))
wsbm_measures=select(main_schaeferwsbm,one_of("age", "system_segreg_wsbm","mean_within_sys_wsbm", "mean_between_sys_wsbm","modul_wsbm"))
yeodev_measures=select(main_yeodev,one_of("age", "system_segreg_yeodev","mean_within_sys_yeodev", "mean_between_sys_yeodev"))

par(mfrow=c(3,3))
i=1
for (meas in yeo_measures){
  name=colnames(yeo_measures)[i]
  i=i+1
  p <- vioplot(meas~main_schaeferyeo7$gender, main=name) #gender
}


#####################################
########### TEST SAMPLE #########
#####################################
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/'

#Get network connectivity data
#From Schaefer400-Yeo7
net_stats_schaeferyeo7 <- read.csv(paste0(net_data_dir,"n544_test_sample_schaefer400_yeo7_network_stats.csv"))
#From Schaefer400-WSBM
net_stats_schaeferwsbm <- read.csv(paste0(net_data_dir,"n544_site14site20_test_sample_schaefer400_wsbm_network_stats.csv"))
#From fsaverage6-Yeo dev
net_stats_yeodev <- read.csv(paste0(net_data_dir,"n544_test_sample_fsaverage6_yeodev_network_stats.csv"))

#remake main
main<- left_join(sites,socio, by=c("ID", "eventname"))
main <- left_join(main, income, by=c("ID", "eventname"))
main <- left_join(main, nihtbx, by=c("ID", "eventname"))
main <- left_join(main, wiscv, by=c("ID", "eventname"))
main <- left_join(main, taskfmri, by=c("ID", "eventname"))
main <- left_join(main, recmem, by=c("ID", "eventname"))

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
#make site a factor
main$site <- as.factor(main$site_id_l)
#average nback face measures
main$nback_2back_place_rate_correct <- main$tfmri_nb_all_beh_c2bp_rate #is rate of correct for 2-back places
main$nback_2back_face_rate_correct <- (main$tfmri_nb_all_beh_c2bpf_rate+main$tfmri_nb_all_beh_c2bnf_rate+main$tfmri_nb_all_beh_c2bngf_rate)/3 #average together all the face rates
#average dprime measures
main$dprime_faces <- (main$tfmri_rec_all_beh_posf_dpr+main$tfmri_rec_all_beh_neutf_dp+main$tfmri_rec_all_beh_negf_dp)/3
main$dprime_places <- main$tfmri_rec_all_beh_place_dp
#take out the one outlier in dprime_faces
main <- filter(main, dprime_faces > -3)

main_schaeferyeo7 <- left_join(net_stats_schaeferyeo7, main, by="ID")
main_schaeferyeo7 <- left_join(main_schaeferyeo7, qa, by="ID")
main_schaeferwsbm <- left_join(net_stats_schaeferwsbm, main, by="ID")
main_schaeferwsbm <- left_join(main_schaeferwsbm, qa, by="ID")
main_yeodev<- left_join(net_stats_yeodev, main, by="ID")
main_yeodev <- left_join(main_yeodev, qa, by="ID")

#take out the one outlier in mean within_sys_conn
main_schaeferyeo7 <- filter(main_schaeferyeo7, mean_within_sys_yeo < 0.6) 

main_schaeferyeo7 <- filter(main_schaeferyeo7, ID != "sub-NDARINVJV77KDEJ")
main_schaeferwsbm <- filter(main_schaeferwsbm, ID != "sub-NDARINVJV77KDEJ")
main_yeodev <- filter(main_yeodev, ID != "sub-NDARINVJV77KDEJ")

# Plot data descriptives --------------------------------------------------

#take out those subjects who were flagged for poor performance on the nback, about 5% of them
main_schaeferyeo7_nback <- filter(main_schaeferyeo7, tfmri_nback_beh_performflag==1 | is.na(tfmri_nback_beh_performflag))# %>% filter(., nback_2back_place_rate_correct > 0.6)
main_schaeferwsbm_nback <- filter(main_schaeferwsbm, tfmri_nback_beh_performflag==1 | is.na(tfmri_nback_beh_performflag))# %>% filter(., nback_2back_place_rate_correct > 0.6)
main_yeodev_nback <- filter(main_yeodev, tfmri_nback_beh_performflag==1 | is.na(tfmri_nback_beh_performflag))

measures=select(main_schaeferyeo7_nback,one_of("nihtbx_list_uncorrected","pea_wiscv_tss", "tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate"))
par(mfrow=c(4,3))
i=1
for (meas in measures){
  name=colnames(measures)[i]
  hist(meas, main=name, col = "lightblue") #hist of measure
  i=i+1
  scatter.smooth(main_schaeferyeo7_nback$age,meas,  col = "blue", main="versus age")
  p <- vioplot(meas~main_schaeferyeo7_nback$gender, main="versus gender") #gender
}

# Cognitive measures-------------------------------------------------------
measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "dprime_places", "dprime_faces")
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7", "sys5to6", "sys5to7")
nets=c("system_segreg_yeo","mean_within_sys_yeo", "mean_between_sys_yeo","modul_yeo")
#Schaefer400-Yeo7
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+site+', net))
    assign(name, lm(formula, data=main_schaeferyeo7))
    print(summary(get(name)))
    print(lm.beta(get(name)))
  }
}

par(mfrow=c(3,4))
for (net in nets){
  visreg(get(paste0("lm_dprime_faces_",net)), net)
  #avPlot(get(paste0("lm_dprime_faces_",net)), net)
  #plot_model(get(paste0("lm_dprime_faces_",net)), terms = net, type="eff")
}

measures=c("tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate", "nback_2back_place_rate_correct","nback_2back_face_rate_correct")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(net,'~age+gender+fd_mean_avg+avgweight+site+', meas))
    assign(name, lm(formula, data=main_schaeferyeo7_nback))
    print(summary(get(name)))
  }
}

#Schaefer400-WSBM
measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "dprime_places", "dprime_faces")
nets=c("system_segreg_wsbm","mean_within_sys_wsbm", "mean_between_sys_wsbm","modul_wsbm")
#nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7", "sys5to6")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+site+', net))
    assign(name, lm(formula, data=main_schaeferwsbm))
    print(summary(get(name)))
  }
}

for (net in nets){
  visreg(get(paste0("lm_dprime_faces_",net)), net)
  #avPlot(get(paste0("lm_dprime_faces_",net)), net)
  #plot_model(get(paste0("lm_dprime_faces_",net)), terms = net, type="eff")
}

measures=c("tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate", "nback_2back_place_rate_correct","nback_2back_face_rate_correct")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+site+', net))
    assign(name, lm(formula, data=main_schaeferwsbm_nback))
    print(summary(get(name)))
  }
}

visreg(lm_nback_2back_place_rate_correct_sys4to7)

#Yeo-dev
measures=c("nihtbx_list_uncorrected","pea_wiscv_tss", "dprime_places", "dprime_faces")
nets=c("system_segreg_yeodev","mean_within_sys_yeodev", "mean_between_sys_yeodev")
nets=c("sys4to4","sys6to6", "sys4to7", "sys6to7", "sys5to6")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+site+', net))
    assign(name, lm(formula, data=main_yeodev))
    print(summary(get(name)))
  }
}
for (net in nets){
  visreg(get(paste0("lm_dprime_faces_",net)), net)
  #avPlot(get(paste0("lm_dprime_faces_",net)), net)
  #plot_model(get(paste0("lm_dprime_faces_",net)), terms = net, type="eff")
}

measures=c("tfmri_nb_all_beh_ctotal_rate","tfmri_nb_all_beh_c2b_rate", "nback_2back_place_rate_correct","nback_2back_face_rate_correct")
for (meas in measures){
  print(meas)
  for (net in nets){
    print(net)
    name<-paste0("lm_", meas,"_",net)
    formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+avgweight+site+', net))
    assign(name, lm(formula, data=main_yeodev_nback))
    print(summary(get(name)))
  }
}

# Age effects-descriptive -------------------------------------------------------------
yeo_measures=select(main_schaeferyeo7,one_of("age", "system_segreg_yeo","mean_within_sys_yeo", "mean_between_sys_yeo","modul_yeo"))
wsbm_measures=select(main_schaeferwsbm,one_of("age", "system_segreg_wsbm","mean_within_sys_wsbm", "mean_between_sys_wsbm","modul_wsbm"))
yeodev_measures=select(main_yeodev,one_of("age", "system_segreg_yeodev","mean_within_sys_yeodev", "mean_between_sys_yeodev"))

par(mfrow=c(5,3))
i=1
for (meas in yeo_measures){
  name=colnames(measures)[i]
  hist(meas, main=name, col = "lightblue") #hist of measure
  i=i+1
  scatter.smooth(main_schaeferyeo7$age,meas,  col = "blue", main="versus age")
  p <- vioplot(meas~main_schaeferyeo7$gender, main="versus gender") #gender
}
par(mfrow=c(3,5))
for (meas in names(yeodev_measures)){
  name=meas
  print(meas)
  name<-paste0("lm_", meas)
  formula<-formula(paste0(meas,'~age+gender+fd_mean_avg+site+avgweight'))
  assign(name, lm(formula, data=main_yeodev))
  visreg(get(name), "age", main=meas)
  print(summary(get(name)))
}

#plot each of them with each of the others
chart.Correlation(yeo_measures)
chart.Correlation(wsbm_measures)
chart.Correlation(yeodev_measures)

# Age and sex effects-comparison --------------------------------------------------
# Make dataframes for each partitions with the strength of age effects and pvals!
covariates="~age+gender+fd_mean_avg+site+avgweight" #to look at age
covariates="~gender+age+fd_mean_avg+site+avgweight" #to look at sex


partitions=c("schaeferyeo7","schaeferwsbm", "yeodev")
ends=c("yeo", "wsbm", "yeodev")
i=1
for (partition in partitions){
  end=ends[i]
  main_unique <- dplyr::select(get(paste0("main_",partition)), -c(sys2to1,sys3to1,sys3to2,sys4to1,sys4to2,sys4to3,sys5to1,sys5to2,sys5to3,sys5to4,sys6to1,sys6to2,sys6to3,sys6to4,sys6to5,sys7to1,sys7to2,sys7to3,sys7to4,sys7to5,sys7to6))
  #run and compare for multiple comparisons again
  cols <- names((dplyr::select(ungroup(main_unique),matches("system_segreg| mean_within_sys_ | mean_between_sys_|sys"))))
  m <- mclapply(cols, function(sys) {as.formula(paste(sys, covariates, sep=""))},mc.cores=2)
  #betas
  networks_Age_betas <- mclapply(m, function(sys) {summary(lm(formula = sys,data=main_unique))$coefficients[2,1]},mc.cores=1)
  networks_Age_betas <- as.data.frame(networks_Age_betas)
  networks_Age_betas <- t(networks_Age_betas)
  networks_Age_betas <- as.data.frame(networks_Age_betas)
  #pvals
  networks_Age_pvals <- mclapply(m, function(sys) { summary(lm(formula = sys,data=main_unique))$coefficients[2,4]},mc.cores=1)
  networks_Age_pvals <- as.data.frame(networks_Age_pvals)
  networks_Age_pvals <- t(networks_Age_pvals)
  networks_Age_pvals <- as.data.frame(networks_Age_pvals)
  #bonferroni correct
  networks_Age_pvals_fdr <- p.adjust(networks_Age_pvals$V1, method="fdr")
  networks_age_fdr <- data.frame(networks_Age_betas[,1],networks_Age_pvals_fdr,c("system_segreg","mean_within_sys","mean_between_sys",names(dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
  colnames(networks_age_fdr) <- c("beta","pvalue", "network")
  networks_age_fdr$partition <- rep(end, length(networks_age_fdr$beta))
  #FDR correction shows DMN to attentional networks and visual to dorsal attention is marginal.
  print(networks_age_fdr)
  assign(paste0("networks_age_fdr_",end), networks_age_fdr)
  i=i+1
}

#Plot the comparison
full <- rbind(networks_age_fdr_yeo,networks_age_fdr_wsbm, networks_age_fdr_yeodev)
dim(full)
colnames(full)

#reorder them to have the right order
full$partition <- factor(full$partition,levels = ends)
full$network <- factor(full$network, levels=c("system_segreg","mean_within_sys","mean_between_sys",names(dplyr::select(ungroup(main_unique),sys1to1:sys7to7))))
ggplot(aes(x = partition, y = network, fill = beta), data=full)+geom_tile()+scale_fill_gradient2(low = ("lightblue"),
                                                                                                 mid="white",high = "red", midpoint=0, limits=c(-0.021, 0.021))+theme_ipsum()+geom_text(aes(label=ifelse(pvalue<0.05, "*", ""))) 
