library(dplyr)
library(psych)
library(mgcv)
library(stringi)
library(stringr)
library(summarytools)
library(lm.beta)
library(ggplot2)
require(reshape2)
library(lme4)
# SETUP -------------------------------------------------------------------
#Cluster mounted locally on personal computer
sublistdir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20/"
netdatadir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/"
wsbmdatadir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site14site20_test_sample/"
qadir="/Users/utooley/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/"
demo_dir="~/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjData/Release2_fixed/"

# Read in files -----------------------------------------------------------
wsbm_consensus<-read.csv(paste0(wsbmdatadir, "n546_test_sample_schaefer400_wsbm_consensus_iter_network_stats.csv"))
wsbm_similarity<-read.csv(paste0(wsbmdatadir, "n546_test_sample_schaefer400_wsbm_consensus_network_stats.csv"))
yeo_nets <- read.csv(paste0(netdatadir, "n546_test_sample_schaefer400_yeo_network_stats.csv"))
qaData <- read.csv(paste0(qadir,"XCP_QAVARS.csv")) #I only took the first two runs here so make sure to only keep those variables
training_subjlist <- read.csv("~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/parcellation/n670_filtered_runs_site16_postprocess.csv")
subjlist <- read.csv(paste0(sublistdir, "n546_filtered_runs_site14site20_postprocess.csv")) #I only took the first two runs here so make sure to only keep those variables and then mutate them
subjlist <- subjlist %>% select(., id:var2)
subjlist <- reshape(subjlist, direction="long", varying=list(2:3)) %>% arrange(id) %>% rename(run=var1, index=time, ID=id)

## READ IN DEMO DATA ##
# right now just getting age and sex
demoData <- read.csv(paste0(demo_dir, "pdem02.txt"), sep = "\t", stringsAsFactors = F)
demoData <- demoData[-1,] 
demoData <- droplevels(demoData)
demoData$age <- as.numeric(demoData$interview_age)/12
demoData$gender <- as.factor(demoData$gender)
demoData <- select(demoData,src_subject_id:demo_prim, age, gender)
demoData$ID <- demoData$src_subject_id
demoData$ID <- paste0("sub-",gsub("_", "", demoData$ID))
# Summarise any run-wise statistics ---------------------------------------

## FIND TOTALSIZET FOR PARTICIPANTS WITH NO CENSORED VOLS BY HAND! ##
qaData$nVolCensored[is.na(qaData$nVolCensored)]<- 0
qaData <- rename(qaData, ID=id0, run=id1)
#filter out extraneous QA variables and make summary variables of volumes and censored volumes
qaData <- qaData %>% group_by(ID) %>% 
  mutate(totalnVolCensored=sum(nVolCensored), totalSizet=sum(nVols)) %>% ungroup()

qaData <- qaData %>% mutate(perc_vols=nVols/totalSizet, relMeanRMS_weight=relMeanRMSMotion*perc_vols, pctSpikesFD_weight=pctSpikesFD*perc_vols) %>% group_by(ID) %>% 
  mutate(relMeanRMS_avg=sum(relMeanRMS_weight), pctVolsCensored=(totalnVolCensored/totalSizet), pctSpikesFD_avg=sum(pctSpikesFD_weight))
#average motion across the two runs, weighted by the length of each run as a percentage of the total, same for percent spikes FD.

# Data Cleaning -----------------------------------------------------------
#In training sample sub-NDARINVK0U3PRRR
#In test sample, it's  sub-NDARINV2RD4CZ7T

qaData <- left_join(subjlist,qaData,by = c("ID", "run")) 
#make ID a character vector
subjlist$ID <- as.character(subjlist$ID)
wsbm_consensus$ID <- as.character(wsbm_consensus$subjlist)
wsbm_similarity$ID <- as.character(wsbm_similarity$subjlist)
yeo_nets$ID <- as.character(yeo_nets$subjlist)

master <- merge(wsbm_consensus, yeo_nets, by=c("ID"), suffixes = c("_wsbm_consensus", "_yeo"))
master_similarity<- merge(wsbm_similarity, yeo_nets,by="ID", suffixes=c("_similarity", "_yeo")) #make this better

master <- left_join(master,demoData,by="ID")
master <- left_join(master,qaData,by="ID")  %>% group_by(ID) %>% filter(row_number() == 1)

# T.tests of net metrics between 3 partitions ------------------------------------------
#SD of edge weights
t.test(master$deviation_edge_weights_consensus, master$deviation_edge_weights_yeo,paired=TRUE)#higher in yeo
t.test(master$deviation_edge_weights_yeo, master_similarity$deviation_edge_weights_consensus, paired=TRUE)
t.test(master$deviation_edge_weights_consensus, master_similarity$deviation_edge_weights_consensus, paired=TRUE)

#Within and between
t.test(master$mean_within_sys_yeo, master$mean_within_sys_consensus,paired=TRUE) 
t.test(master$mean_between_sys_consensus, master$mean_between_sys_yeo, paired=TRUE)

t.test(master$sub_partcoef_pos_yeo, master$sub_partcoef_pos_consensus,paired=TRUE)
t.test(master$sub_partcoef_neg_yeo, master$sub_partcoef_neg_consensus,paired=TRUE)

#System segregation
t.test(master$system_segreg_consensus,master$system_segreg_yeo,paired=TRUE)
t.test(master$system_segreg_consensus,master_similarity$system_segreg_consensus,paired=TRUE)

#Modularity quality
t.test(master$modul_consensus, master$modul_yeo, paired=TRUE) 

# Make violin plots of net metrics in the 3 partitions -------------------------------------------------------------
## EDGE WEIGHT STANDARD DEVIATION
library(lme4)
library(lmerTest)
df <- melt(data.frame(master$deviation_edge_weights_consensus,master$ID, master_similarity$ID, master_similarity$deviation_edge_weights_consensus, master$deviation_edge_weights_yeo))
colnames(df) <- c("ID","ID2","which", "value")
print(df)
#plot
anova(lmer(value~as.factor(which)+(1|ID), data=df)) #try mixed effects model instead of t-test for comparisons, still sig
ggplot(df, aes(which, value))+ geom_violin(mapping = aes(which, value), color="darkblue", fill="lightblue")+ggtitle(df$which[1])+stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
                                                                                                                                              geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_classic()+scale_color_gradient()
##PARTICIPATION COEF
df <- melt(data.frame(master$sub_partcoef_pos_consensus, master$ID, master_similarity$sub_partcoef_pos_consensus, master_similarity$ID, master$sub_partcoef_pos_yeo))
colnames(df) <- c("ID","ID2","which", "value")
print(df)
#plot
anova(lmer(value~as.factor(which)+(1|ID), data=df))
ggplot(df, aes(which, value))+ geom_violin(binaxis="y", mapping = aes(which, value),stackdir="center", color="darkblue", fill="darkblue")+ggtitle(df$which[1])+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                                            geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_classic()+scale_color_gradient()

df <- melt(data.frame(master$sub_partcoef_neg_consensus, master$ID, master_similarity$sub_partcoef_neg_consensus, master_similarity$ID, master$sub_partcoef_neg_yeo))
colnames(df) <- c("ID","ID2","which", "value")
print(df)
anova(lmer(value~as.factor(which)+(1|ID), data=df))
ggplot(df, aes(which, value))+ geom_violin(binaxis="y", mapping = aes(which, value),stackdir="center", color="darkblue", fill="darkblue")+ggtitle(df$which[1])+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                                            geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_classic()+scale_color_gradient()
## BETWEEN AND WITHIN SYSTEM CONNECTIVITY
df <- melt(data.frame( master$mean_within_sys_consensus, master$ID,  master_similarity$mean_within_sys_consensus, master_similarity$ID, master$mean_within_sys_yeo))
colnames(df) <- c("ID","ID2","which", "value")
print(df)
#plot
anova(lmer(value~as.factor(which)+(1|ID), data=df))
ggplot(df, aes(which, value))+ geom_violin(binaxis="y", mapping = aes(which, value),stackdir="center", color="darkblue", fill="lightyellow")+ggtitle(df$which[1])+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                                               geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_bw()+scale_color_gradient()
df <- melt(data.frame( master$mean_between_sys_consensus, master$ID,  master_similarity$mean_between_sys_consensus, master_similarity$ID, master$mean_between_sys_yeo))
colnames(df) <- c("ID","ID2","which", "value")
print(df)
#plot
anova(lmer(value~as.factor(which)+(1|ID), data=df))
ggplot(df, aes(which, value))+ geom_violin(mapping = aes(which, value),color="darkblue", fill="lightyellow")+ggtitle(df$which[1])+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                               geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_classic()+scale_color_gradient()
## SEGREGATION
df <- melt(data.frame(master$system_segreg_consensus, master$ID, master_similarity$system_segreg_consensus, master_similarity$ID, master$system_segreg_yeo))
colnames(df) <- c("ID","ID2","which", "value")
print(df)
#plot
anova(lmer(value~as.factor(which)+(1|ID), data=df))
ggplot(df, aes(which, value))+ geom_violin(binaxis="y", mapping = aes(which, value),stackdir="center", color="blue", fill="lightyellow")+ggtitle(df$which[1])+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                                           geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_classic()+scale_color_gradient()
## MODULARITY QUALITY
df <- melt(data.frame(master$modul_consensus, master$ID, master_similarity$modul_consensus, master_similarity$ID, master$modul_yeo))
colnames(df) <- c("ID","ID2","which", "value")
print(df)
#plot
anova(lmer(value~as.factor(which)+(1|ID), data=df))
ggplot(df, aes(which, value))+ geom_violin(binaxis="y", mapping = aes(which, value),stackdir="center", color="blue", fill="lightgreen")+ggtitle(df$which[1])+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                                          geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_classic()+scale_color_gradient()
# Descriptives on WSBM Data ----------------------------------------------------
view(dfSummary(master))
describe(master$deviation_edge_weights_yeo)

l2 <- lm(sub_log_evidence ~ age+gender+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(avgweight_yeo ~ age+gender+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(mean_within_sys_consensus ~ age+gender+avgweight_yeo+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(mean_between_sys_consensus ~ age+gender+avgweight_yeo+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(system_segreg_wsbm ~ age+gender+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_pos_wsbm ~ age+gender+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_neg_wsbm ~ age+gender+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(zRandstat ~ age+gender+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

# Descriptives on Yeo Data ------------------------------------------
l2 <- lm(mean_within_sys_yeo ~ age+gender+avgweight_yeo+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(mean_between_sys_yeo ~ age+gender+avgweight_yeo+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(system_segreg_yeo ~ age+gender+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_pos_yeo ~ age+gender+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_neg_yeo ~ age+gender+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

# Plot variance in community assignment by Yeo system ------------------------
library(R.matlab)
yeo_nodes=read.table('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt', col.names = "yeo_nodes")
yeo_nodes <- as.character(yeo_nodes$yeo_nodes)
wsbmvectorsdatadir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling/"
#Read in the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbmvectorsdatadir,"consensus_iter_mode.mat"), drop = )
freq<- t(z$freq)
freq_continuous <- abs(freq-670) #change into the number of non-modal assignments out of 670
data <- data.frame(yeo_nodes,freq_continuous)
se <- function(x){sd(x)/sqrt(length(x))}
my_dat <- data %>% group_by(yeo_nodes) %>% summarise(my_mean = mean(freq_continuous),my_se = se(freq_continuous))

p<-ggplot(data, aes(x=yeo_nodes, y=freq_continuous)) + geom_jitter(position=position_jitter(0.1), cex=1, alpha=0.2) +geom_bar(data=my_dat, aes(y=my_mean,fill=yeo_nodes, x=yeo_nodes,ymin=my_mean-my_se,ymax=my_mean+my_se), stat="identity", width = 0.75) + 
  geom_errorbar(data=my_dat, aes(y=my_mean,x=yeo_nodes,ymin=my_mean-my_se,ymax=my_mean+my_se), width = 0.5) + theme_classic()+
  scale_fill_manual(values = c("purple","blue","aquamarine","khaki1","magenta","orange", "red"))
p

# Plotting partitions on brains ---------------------------------------------------------
library(ggseg3d)
library(ggsegExtra)
library(tidyr)
library(R.matlab)
library(dplyr)
#Read in the vector of WSBM community assignments
wsbmvectorsdatadir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling/"
partitions <- readMat(paste0(wsbmvectorsdatadir,"n670_training_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_representative_labels <-  partitions$consensus.represent.yeorelabeled
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled

#Read in entropy of co-occurence matrix and the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbmvectorsdatadir,"consensus_iter_mode.mat"), drop = )
entropy<- z$node.entropy
freq<- t(z$freq)

#Change the output directory for figures to be on the cluster brains
setwd("~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/figures/brains/")

#Make a network vector to use for coloring the Yeo brain
someData_Yeonets = schaefer7_3d %>% 
  filter(surf == "inflated") %>% 
  unnest(ggseg_3d) %>% 
  select(annot,network) %>% 
  na.omit() %>% 
  mutate(community = ifelse(network=="Vis", 1, 
                            ifelse(network=="SomMot", 2, 
                                   ifelse(network=="DorsAttn", 3, 
                                          ifelse(network=="SalVentAttn", 4, 
                                                 ifelse(network=="Limbic", 5,
                                                        ifelse(network=="Cont", 6,
                                                               ifelse(network=="Default", 7, NA)))))))) #needed to make this a numeric vector not a label

someData = schaefer7_3d %>% 
  filter(surf == "inflated") %>% 
  unnest(cols=ggseg_3d) %>% 
  select(annot, network) %>% 
  na.omit() %>% 
  mutate(consensus_iter=rbind(0,consensus_iterative_labels),consensus_represent=rbind(0, consensus_representative_labels), freq=rbind(0,freq), entropy=rbind(0,entropy))

# YEO PARTITION
hemispheres=c("left", "right")
pans=c("medial", "lateral")
for (hemi in hemispheres){
  for (pan in pans){
    yeo<- ggseg3d(.data = someData_Yeonets, 
                  atlas = schaefer7_3d, hemisphere = hemi,
                  colour = "community", text = "network", palette = c("gray"=0,"purple"=1,"blue"=2,"aquamarine"=3,"khaki1"=4,"magenta"=5,"orange"=6, "red"=7)) %>% 
      pan_camera(paste0(hemi, " ", pan)) %>% remove_axes()
    orca(yeo, paste0("yeo", hemi, pan,".png"))
  }
}

#CONSENSUS ITER PARTITION
hemispheres=c("left", "right")
pans=c("medial", "lateral")
for (hemi in hemispheres){
  for (pan in pans){
    consensus_iter <- ggseg3d(.data = someData, 
                              atlas = schaefer7_3d, hemisphere = hemi,
                              colour = "consensus_iter", text = "network", palette = c("gray"=0,"purple"=1,"blue"=2,"aquamarine"=3,"khaki1"=4,"magenta"=5,"orange"=6, "red"=7)) %>% 
      pan_camera(paste0(hemi, " ", pan)) %>% remove_axes()
    orca(consensus_iter, paste0("consensus_iter_", hemi, pan,".png"))
  }
}

#CONSENSUS REPRESENTATIVE/SIMILARITY PARTITION
hemispheres=c("left", "right")
pans=c("medial", "lateral")
for (hemi in hemispheres){
  for (pan in pans){
    consensus_represent <- ggseg3d(.data = someData, 
                                   atlas = schaefer7_3d, hemisphere = hemi,
                                   colour = "consensus_represent", text = "network", palette = c("gray"=0,"purple"=1,"blue"=2,"aquamarine"=3,"khaki1"=4,"magenta"=5,"orange"=6, "red"=7)) %>% 
      pan_camera(paste0(hemi, " ", pan)) %>% remove_axes()
    orca(consensus_represent, paste0("consensus_representative", hemi, pan,".png"))
  }
}

#VARIANCE, TIES IN CONSENSUS COMMUNITY ASSIGNMENT
hemispheres=c("left", "right")
pans=c("medial", "lateral")
for (hemi in hemispheres){
  for (pan in pans){
    ties <- ggseg3d(.data = someData, 
                    atlas = schaefer7_3d, hemisphere = hemi,
                    colour = "freq", text = "network", palette = c("blue"=300,"gray"=670)) %>% 
      pan_camera(paste0(hemi, " ", pan)) %>% remove_axes()
    orca(ties, paste0("ties_comm_assignment_", hemi, pan,".png"))
  }
}

#RGB of colors for MATLAB
colors <- t(col2rgb(c("black","gray","purple","blue","aquamarine","khaki1","magenta","orange", "red")))
# Plot the search over k ------------------------------------------------------
k_dir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/search_over_k/"
logevidence=matrix(NA,670,18)
for (i in 1:18){
  k=i+3
  try(mine <- read.csv(paste0(k_dir, "wsbm_search_over_k",k,"_n670_site16_30trials.csv")))
  logevidence[,i] <- with(mine, if("subjlist_.4" %in% colnames(mine)) subjlist_.4 else Var4)
}
#remove the outlier subject sub-NDARINVK0U3PRRR
colnames(logevidence)=lapply(4:21, function(x) paste0("k",x)) #give them column names
#k10 is full of 0's, filter out those who haven't finished k9
logevidence2 <- filter(data.frame(logevidence), !is.na(logevidence[,6]))
logevidence[logevidence==0] <- NA
longdata <- melt(logevidence, varnames = c("x","k"))
p<-ggplot(longdata, aes(x=k, y=value)) + geom_jitter(position=position_jitter(0.2), cex=1, alpha=0.5)+stat_summary(aes(group=1),
                                                                                                                          fun.y = median, geom="line", col=rgb(101, 58, 150, maxColorValue = 255)) + theme_classic()
p

# Plot system connectivity matrices ---------------------------------------
library(R.matlab)
library(mdpeer)
avg_matrix <- readMat(paste0(wsbmdatadir,"n546_mean_wsbm_consensus_connectivity.mat"))[1]
avgmatrix <- as.matrix(avg_matrix$mean.system.conn.mat)
# palf <- colorRampPalette(c("red","white", "blue")) 
# heatmap(avgmatrix, Rowv = NA, Colv = NA, col = palf(20))
plot(vizu.mat(avgmatrix, fill.scale.limits = c(-1, 1), x.lab="WSBM consensus partition"))

# Archive: Bar graph of nodal variance by Yeo system -------------------------------------------------------------
data=data.frame(t(wsbm_assign_variance))
colnames(data) <- "avgvariance"
#RGB colors of Yeo brain
yeo_in_order_cols <- c(rgb(120, 18, 134, maxColorValue = 255), rgb(70, 130, 180, maxColorValue = 255), rgb(0, 118, 14, maxColorValue = 255), rgb(196, 58, 250, maxColorValue = 255), 
                       rgb(220,248,164, maxColorValue = 255),rgb(230, 148, 34, maxColorValue = 255),  rgb(205, 62, 78, maxColorValue = 255))
Fig_variance <- ggplot(data=data, aes(1:7,avgvariance,)) +geom_bar(fill=yeo_in_order_cols, col="black",stat="identity", size=0)

# Training Sample Symmetry and Similarity to Yeo ----------------------------------------------
wsbmtrainingdatadir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/"

training_sample=read.csv(paste0(wsbmtrainingdatadir,"wsbm_k7_training_sample_n670_site16_50trials_norelabel.csv"))

# Look at how symmetric the hemispheres are in terms of labeling within a subject

#Look at subject similarity to the Yeo partition
