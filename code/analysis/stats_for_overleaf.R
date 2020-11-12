


# Demographic statistics --------------------------------------------------
save(main_yeodev, file="~/Downloads/Demo_data.Rdata")
load("Demo_data.Rdata")

#length
dim(main_yeodev)[1]

#Race
sum(main_yeodev$race_ethnicity==1)/669
sum(main_yeodev$race_ethnicity==2)/669

#Education
main_yeodev$demo_prtnr_ed_v2
#parent edu variable
summary(main_yeodev$demo_prnt_ed_v2)
main_yeodev$demo_prnt_ed_num_v2 <- as.numeric(main_yeodev$demo_prnt_ed_v2)
main_yeodev$demo_prnt_ed_num_v2_recoded <- recode(as.character(main_yeodev$demo_prnt_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_yeodev$demo_prnt_ed_num_v2_recoded <- as.numeric(main_yeodev$demo_prnt_ed_num_v2_recoded)

#partner edu variable
summary(main_yeodev$demo_prtnr_ed_v2)
main_yeodev$demo_prtnr_ed_num_v2 <- as.numeric(main_yeodev$demo_prtnr_ed_v2)
main_yeodev$demo_prtnr_ed_num_v2_recoded <- recode(as.character(main_yeodev$demo_prtnr_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_yeodev$demo_prtnr_ed_num_v2_recoded <- as.numeric(main_yeodev$demo_prtnr_ed_num_v2_recoded)

#creating a combined parent education variable
main_yeodev$demo_comb_parent_edu <- rowMeans(main_yeodev[c('demo_prnt_ed_num_v2_recoded', 'demo_prtnr_ed_num_v2_recoded')], na.rm=TRUE)
View(main_yeodev$demo_comb_parent_edu)
hist(as.numeric(main_yeodev$demo_comb_parent_edu), xlab = "Combined averaged years of parent education", breaks = 15, label = TRUE, col = "red")
median(main_yeodev$demo_comb_parent_edu)


# Partition comparisons ---------------------------------------------------
wsbm_full <- c(wsbm_lh,wsbm_rh)
#save(yeo_dev_full, yeo_full, wsbm_full, file="~/Downloads/Partition_comparisons.RData")
save(yeo_dev_full, yeo_full, wsbm_full,NMI_yeodev_yeo7, NMI_wsbm_yeo7, NID_wsbm_yeo7, NID_yeodev_yeo7, file="~/Downloads/Partition_comparison_stats.RData")

#NMI
NMI_yeodev_yeo7<- compare(yeo_dev_full, yeo_full, method = "nmi")
NID_yeodev_yeo7<- clustComp(yeo_dev_full, yeo_full)$NID

NMI_wsbm_yeo7<- compare(wsbm_full, yeo_full, method = "nmi")
NID_wsbm_yeo7<- clustComp(wsbm_full, yeo_full)$NID


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

