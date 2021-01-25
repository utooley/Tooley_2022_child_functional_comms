library(mclustcomp)
library(fsbrain)
library(R.matlab)
library(igraph)
library(aricode)
library(data.table)

# Load in data to make demo data ------------------------------------------
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

main<- left_join(sites,socio, by=c("ID", "eventname"))
main <- left_join(main, income, by=c("ID", "eventname"))

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
subjlist1 <- read.table(paste0(subjlist_dir, "/n670_subjects_only_filtered_runs_site16_postprocess.txt"), col.names = c("ID"))
#Keep only n670 subject list
main$ID <- str_remove(main$ID, "_") #subjectid has a _ in it after NDAR in the NDAR data, doesn't match subject list
main$ID <- paste0("sub-",main$ID)
#make age numeric
main$age <- as.numeric(main$interview_age.x)/12
#make sex a factor
main$gender <- as.factor(main$gender.x)

main_schaeferyeo7 <- left_join(subjlist1, main, by="ID")
main_schaeferyeo7 <- left_join(main_schaeferyeo7, qa, by="ID")
main_schaeferwsbm <- left_join(subjlist1, main, by="ID")
main_schaeferwsbm <- left_join(main_schaeferwsbm, qa, by="ID")
main_yeodev <- left_join(subjlist1, main, by="ID")
main_yeodev <- left_join(main_yeodev, qa, by="ID")

# Demographic statistics --------------------------------------------------
#length
dim(main_yeodev)[1]

#Race
sum(main_yeodev$race_ethnicity==1)/670
sum(main_yeodev$race_ethnicity==2)/670

#Education
main_yeodev$demo_prtnr_ed_v2
#parent edu variable
summary(main_yeodev$demo_prnt_ed_v2)
main_yeodev$demo_prnt_ed_num_v2 <- as.numeric(main_yeodev$demo_prnt_ed_v2)
main_yeodev$demo_prnt_ed_num_v2_recoded <- car::recode(as.character(main_yeodev$demo_prnt_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_yeodev$demo_prnt_ed_num_v2_recoded <- as.numeric(main_yeodev$demo_prnt_ed_num_v2_recoded)

#partner edu variable
summary(main_yeodev$demo_prtnr_ed_v2)
main_yeodev$demo_prtnr_ed_num_v2 <- as.numeric(main_yeodev$demo_prtnr_ed_v2)
main_yeodev$demo_prtnr_ed_num_v2_recoded <- car::recode(as.character(main_yeodev$demo_prtnr_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_yeodev$demo_prtnr_ed_num_v2_recoded <- as.numeric(main_yeodev$demo_prtnr_ed_num_v2_recoded)

#creating a combined parent education variable
main_yeodev$demo_comb_parent_edu <- rowMeans(main_yeodev[c('demo_prnt_ed_num_v2_recoded', 'demo_prtnr_ed_num_v2_recoded')], na.rm=TRUE)
View(main_yeodev$demo_comb_parent_edu)
hist(as.numeric(main_yeodev$demo_comb_parent_edu), xlab = "Combined averaged years of parent education", breaks = 15, label = TRUE, col = "red")
median(main_yeodev$demo_comb_parent_edu)
# 
# save(main_yeodev, file="~/Downloads/Demo_data.Rdata")
# load("Demo_data.Rdata")
# Replication dataset demographics ----------------------------------------
data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
output_data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/'
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

main<- left_join(sites,socio, by=c("ID", "eventname"))
main <- left_join(main, income, by=c("ID", "eventname"))

##need to get XCP mean FD and # of outliers and control for that
runs <- read.csv(paste0(subjlist_dir,'n544_filtered_runs_site14site20_postprocess.csv'))
runs <- select(runs, id, var1:var2) #take only the first 2 runs that were used
runs<- melt(runs, measure=c("var1", "var2")) %>% arrange(., id) %>% select(., -variable) %>% rename(.,ID=id,run=value)#reshape them and take out extra
#read in the list of runs that were used and merge it with xcp qa vars
qa_vars <- read.csv(paste0(raw_data_dir, "bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/XCP_QAVARS_with_nVols.csv"))
qa_vars <- rename(qa_vars, ID=id0, run=id1)
qa <- left_join(runs, qa_vars, by=c("ID", "run"))

#make summary variables of volumes and censored volumes
qa <- qa %>% group_by(ID) %>% mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(nVols)) %>% ungroup()

qa <- qa %>% mutate(perc_vols=nVols/totalSizet, relMeanRMSMotion_weight=relMeanRMSMotion*perc_vols, pctSpikesFD_weight=pctSpikesFD*perc_vols) %>% group_by(ID) %>% 
  mutate(fd_mean_avg=sum(relMeanRMSMotion_weight), pctVolsCensored=(totalnVolCensored/totalSizet), pctSpikesFD_avg=sum(pctSpikesFD_weight)) %>% select(fd_mean_avg, pctVolsCensored, pctSpikesFD_avg, totalSizet)
#average motion across the two runs, weighted by the length of each run as a percentage of the total, same for percent spikes FD. select only the averaged variables
qa <- qa[!duplicated(qa$ID), ]

#Only baseline visits from cognition
main <- filter(main, eventname=="baseline_year_1_arm_1") %>% filter(.,site_id_l=="site14" | site_id_l=="site20") #filter out only the baseline visits
#for site 16 only
subjlist1 <- read.table(paste0(subjlist_dir, "/n544_subjects_only_filtered_runs_site14site20_postprocess.txt"), col.names = c("ID"))
#Keep only n670 subject list
main$ID <- str_remove(main$ID, "_") #subjectid has a _ in it after NDAR in the NDAR data, doesn't match subject list
main$ID <- paste0("sub-",main$ID)
#make age numeric
main$age <- as.numeric(main$interview_age.x)/12
#make sex a factor
main$gender <- as.factor(main$gender.x)

main_rep_schaeferyeo7 <- left_join(subjlist1, main, by="ID")
main_rep_schaeferyeo7 <- left_join(main_rep_schaeferyeo7, qa, by="ID")
main_rep_schaeferwsbm <- left_join(subjlist1, main, by="ID")
main_rep_schaeferwsbm <- left_join(main_rep_schaeferwsbm, qa, by="ID")
main_rep_yeodev <- left_join(subjlist1, main, by="ID")
main_rep_yeodev <- left_join(main_rep_yeodev, qa, by="ID")

#Race
sum(main_rep_yeodev$race_ethnicity==1)/544
sum(main_rep_yeodev$race_ethnicity==2)/544
sum(main_rep_yeodev$race_ethnicity==3)/544
sum(main_rep_yeodev$race_ethnicity==5)/544

#Education
main_rep_yeodev$demo_prtnr_ed_v2
#parent edu variable
summary(main_rep_yeodev$demo_prnt_ed_v2)
main_rep_yeodev$demo_prnt_ed_num_v2 <- as.numeric(main_rep_yeodev$demo_prnt_ed_v2)
main_rep_yeodev$demo_prnt_ed_num_v2_recoded <- car::recode(as.character(main_rep_yeodev$demo_prnt_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_rep_yeodev$demo_prnt_ed_num_v2_recoded <- as.numeric(main_rep_yeodev$demo_prnt_ed_num_v2_recoded)

#partner edu variable
summary(main_rep_yeodev$demo_prtnr_ed_v2)
main_rep_yeodev$demo_prtnr_ed_num_v2 <- as.numeric(main_rep_yeodev$demo_prtnr_ed_v2)
main_rep_yeodev$demo_prtnr_ed_num_v2_recoded <- car::recode(as.character(main_rep_yeodev$demo_prtnr_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_rep_yeodev$demo_prtnr_ed_num_v2_recoded <- as.numeric(main_rep_yeodev$demo_prtnr_ed_num_v2_recoded)

#creating a combined parent education variable
main_rep_yeodev$demo_comb_parent_edu <- rowMeans(main_rep_yeodev[c('demo_prnt_ed_num_v2_recoded', 'demo_prtnr_ed_num_v2_recoded')], na.rm=TRUE)
View(main_rep_yeodev$demo_comb_parent_edu)
hist(as.numeric(main_rep_yeodev$demo_comb_parent_edu), xlab = "Combined averaged years of parent education", breaks = 15, label = TRUE, col = "red")
median(main_rep_yeodev$demo_comb_parent_edu)

save(main_yeodev, main_rep_yeodev, file="~/Downloads/Demo_data.Rdata")
load("Demo_data.Rdata")

# Partition comparisons ---------------------------------------------------
#Yeo7
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";
lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','Yeo2011_7Networks_N1000' )
yeo7_lh <- as.numeric(str_split(lh$label_names, "_", simplify = TRUE)[,2])
yeo7_lh[is.na(yeo7_lh)] <- 0
rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','Yeo2011_7Networks_N1000' )
yeo7_rh <- as.numeric(str_split(rh$label_names, "_", simplify = TRUE)[,2])
yeo7_rh[is.na(yeo7_rh)] <- 0
yeo_full <- c(yeo7_lh, yeo7_rh)

#Developmental clustering
yeo_dev_dir="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_mot_outliers/"
yeo_dev_partition <- readMat(paste0(yeo_dev_dir,"yeo7_n670_2runsonly_1000tries_mot_outliers.mat"), drop = )
yeo_dev_lh <- yeo_dev_partition$lh.labels
yeo_dev_rh <- yeo_dev_partition$rh.labels
yeo_dev_full <- c(yeo_dev_lh,yeo_dev_rh)

#copy WSBM annotation into local CBIG subjects dir
wsbm_rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.consensus.fsaverage6')
wsbm_rh<- as.numeric(as.character(wsbm_rh$label_names))
wsbm_lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.consensus.fsaverage6')
wsbm_lh<- as.numeric(as.character(wsbm_lh$label_names))
wsbm_full <- c(wsbm_lh,wsbm_rh)

#NMI
NMI_yeodev_yeo7<- compare(yeo_dev_full, yeo_full, method = "nmi")
NID_yeodev_yeo7<- clustComp(yeo_dev_full, yeo_full)$NID

NMI_wsbm_yeo7<- compare(wsbm_full, yeo_full, method = "nmi")
NID_wsbm_yeo7<- clustComp(wsbm_full, yeo_full)$NID

#Replication dataset yeodev
yeo_dev_replication="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n544_2runsonly_1000tries/"
yeo_dev_partition_rep <- readMat(paste0(yeo_dev_replication,"yeo7_n544_2runsonly_1000tries.mat"), drop = )
yeo_dev_rep_lh <- yeo_dev_partition_rep$lh.labels;yeo_dev_rep_rh <- yeo_dev_partition_rep$rh.labels
yeo_dev_rep <- c(yeo_dev_rep_lh,yeo_dev_rep_rh)
compare <- as.numeric(yeo_dev_full==yeo_dev_rep)

NMI_yeodev_yeodevrep<- compare(yeo_dev_full, yeo_dev_rep, method = "nmi")
NID_yeodev_yeodevrep<- clustComp(yeo_dev_full, yeo_dev_rep)$NID

#Original WSBM parcel-wise
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
partitions <- readMat(paste0(wsbm_datadir,"n670_training_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled
schaefer_atlas_region_names_lh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");
wsbm_consensus_lh=as.list(setNames(c(0, consensus_iterative_labels[1:200]), schaefer_atlas_region_names_lh));wsbm_consensus_rh=as.list(setNames(c(0, consensus_iterative_labels[201:400]), schaefer_atlas_region_names_rh))

#Replication dataset WSBM parcel-wise
wsbm_datadir_rep="/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample/brains/"
partitions <- readMat(paste0(wsbm_datadir_rep,"n544_test_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_iterative_labels_replication <- partitions$consensus.iter.mode.yeorelabeled
wsbm_consensus_lh_rep=as.list(setNames(c(0, consensus_iterative_labels_replication[1:200]), schaefer_atlas_region_names_lh));wsbm_consensus_rh_rep=as.list(setNames(c(0, consensus_iterative_labels_replication[201:400]), schaefer_atlas_region_names_rh))
compar_lh <- as.numeric(as.numeric(wsbm_consensus_lh)==as.numeric(wsbm_consensus_lh_rep));compar_rh <- as.numeric(as.numeric(wsbm_consensus_rh)==as.numeric(wsbm_consensus_rh_rep))

NMI_wsbm_wsbmrep<- compare(as.numeric(c(wsbm_consensus_lh,wsbm_consensus_rh)),as.numeric(c(wsbm_consensus_lh_rep,wsbm_consensus_rh_rep)), method = "nmi")
NID_wsbm_wsbmrep<- clustComp(as.numeric(c(wsbm_consensus_lh,wsbm_consensus_rh)),as.numeric(c(wsbm_consensus_lh_rep,wsbm_consensus_rh_rep)))$NID

#save(yeo_dev_full, yeo_full, wsbm_full, file="~/Downloads/Partition_comparisons.RData")
save(yeo_dev_full, yeo_full, wsbm_full,yeo_dev_rep,NMI_yeodev_yeo7, NMI_wsbm_yeo7, NID_wsbm_yeo7, NID_yeodev_yeo7,NMI_yeodev_yeodevrep, NID_yeodev_yeodevrep,NMI_wsbm_wsbmrep, NID_wsbm_wsbmrep, file="~/Downloads/Partition_comparison_stats.RData")

# Confidence data ---------------------------------------------------------

#silhouettes 
silhouette_lh_dev <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh_dev <-  yeo_dev_partition$rh.s

subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
silhouette_lh_yeo7 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/lh.silhouette.fsaverage6.curv"))
silhouette_rh_yeo7 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/rh.silhouette.fsaverage6.curv"))

#Yeo7 adult confidence stats
silhouette_yeo <- c(silhouette_lh_fs6,silhouette_rh_fs6)
yeo7 <- c(yeo7_lh, yeo7_rh)
data <- data.frame(silhouette_yeo, as.character(yeo7))
data <- data[data$as.character.yeo7. != "0",]#remove the medial wall
colnames(data) <- c("silhouette", "yeo7")
data <- data %>% filter(silhouette!=1) #remove the 1's for silhouette around the medial wall

kruskal.test(silhouette~yeo7, data = data)
pairwise.wilcox.test(data$silhouette, data$yeo7,
                     p.adjust.method = "bonferroni")

#Yeodev child confidence stats
silhouette_lh_dev <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh_dev <-  yeo_dev_partition$rh.s
silhouette_dev <- rbind(silhouette_lh_dev,silhouette_rh_dev)
yeodev <- c(yeo_dev_lh, yeo_dev_rh)
data <- data.frame(silhouette_dev, as.character(yeodev))
data <- data[data$as.character.yeodev. != "0",]#remove the medial wall
colnames(data) <- c("silhouette", "yeodev")
data <- data %>% filter(silhouette!=1) #remove the 1's for silhouette around the medial wall

kruskal.test(silhouette~yeodev, data = data)
pairwise.wilcox.test(data$silhouette, data$yeodev,
                     p.adjust.method = "bonferroni")

#Yeo7 child confidence stats
silhouette_lh_dev <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh_dev <-  yeo_dev_partition$rh.s
# silhouette_lh_dev <- silhouette_lh_yeo7
# silhouette_rh_dev <- silhouette_rh_yeo7
silhouette_dev <- c(silhouette_lh_dev,silhouette_rh_dev) 
yeo7 <- c(yeo7_lh, yeo7_rh)
data <- data.frame(silhouette_dev, as.character(yeo7))
data <- data[data$as.character.yeo7 != "0",]#remove the medial wall
colnames(data) <- c("silhouette", "yeo7")
data <- data %>% filter(silhouette!=1) #remove the 1's for silhouette around the medial wall

kruskal.test(silhouette~yeo7, data = data)
pairwise.wilcox.test(data$silhouette, data$yeo7,
                     p.adjust.method = "bonferroni")

#Variability in assignment by WSBM partition
#Read in the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
freq<- t(z$freq)
freq_continuous <- abs(freq-670) #change into the number of non-modal assignments out of 670
partitions <- readMat(paste0(wsbm_datadir,"n670_training_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled
data <- data.frame(as.character(consensus_iterative_labels),freq_continuous)
colnames(data) <- c('wsbm', 'freq')

#Stats
kruskal.test(freq~wsbm, data = data)
pairwise.wilcox.test(data$freq, data$wsbm,
                     p.adjust.method = "bonferroni")

#Variability in assignment by Yeo7 partition
yeo_nodes=read.table('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt', col.names = "yeo_nodes")
yeo_nodes <- as.character(yeo_nodes$yeo_nodes)
#Read in the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
freq<- t(z$freq)
freq_continuous <- abs(freq-670) #change into the number of non-modal assignments out of 670
data <- data.frame(yeo_nodes,freq_continuous)
colnames(data) <- c('yeo7', 'freq')

#Stats
kruskal.test(freq~yeo7, data = data)
pairwise.wilcox.test(data$freq, data$yeo7,
                     p.adjust.method = "bonferroni")


#Task betas for FP
data <- data.frame(yeo7_6_betas, yeo_dev_6_betas, wsbm_6_betas)
colnames(data) <- c('yeo7', 'yeo_dev', 'wsbm')
longdata <- melt(data)
longdata <- longdata[longdata$value != 0.0000000,]#remove the medial wall
longdata <- longdata[longdata$value > 0,]#remove the medial wall
colnames(longdata) <- c('partition', 'betas')
longdata$partition <- factor(longdata$partition,levels = c("wsbm","yeo_dev","yeo7"))
longdata = longdata %>% 
  dplyr::group_by(partition) %>% 
  mutate(med = median(betas))


longdata %>% group_by(partition) %>% summarise_all(mean)
longdata %>% group_by(partition) %>% summarise_all(sd)
#Stats
kruskal.test(betas~partition, data = longdata)
pairwise.wilcox.test(longdata$betas, longdata$partition,
                     p.adjust.method = "bonferroni")

#Task betas for DM
data <- data.frame(yeo7_7_betas, yeo_dev_7_betas, wsbm_7_betas)
colnames(data) <- c('yeo7', 'yeo_dev', 'wsbm')
longdata <- melt(data)
longdata <- longdata[longdata$value != 0.0000000,]#remove the medial wall
longdata <- longdata[longdata$value < 0,]#look at what is most negative and not least positive
colnames(longdata) <- c('partition', 'betas')
longdata$partition <- factor(longdata$partition,levels = c("wsbm","yeo_dev","yeo7"))
longdata = longdata %>% 
  dplyr::group_by(partition) %>% 
  mutate(med = median(betas))

longdata %>% group_by(partition) %>% summarise_all(mean)
longdata %>% group_by(partition) %>% summarise_all(sd)