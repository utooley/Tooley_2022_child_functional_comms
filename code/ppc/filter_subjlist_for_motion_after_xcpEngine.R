# Original sample: Make subject lists ------------------------------------------------------
library(dplyr)
xcp_qa_vars="/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/"
subjlist_vars ="/cbica/projects/spatial_topography/data/subjLists/release2/site16/"

#main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS_with_nVols.csv")) # "XCP_QAVARS.csv")) 
main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS.csv"))
summary(main) #17 NAs, n=696, 2589 runs
length(unique(main$id0))

#Filtering out subjects based on:
#any outliers in coverage metrics (Dice, Jaccard), n=2 runs
main <- main %>% filter(., coregJaccard > 0.95) 
length(unique(main$id0)) #n=0 people, n=2 runs for bad coverage, but also 2589->2578, n=9 runs that had NAs because didn't finish removed

#total number of vols in ABCD to check for runs that are too short
main <- main %>% filter(.,  nVols  > 360 | is.na(nVols)) #I've replaced all the missing nVols from xcpEngine bug so that it reflects the correct number, people with NAs didn't finish.
summary(main)
length(unique(main$id0)) # n= 0 people
dim(main) #2571->2564 n=7 runs removed

#and any individual run with motion over 0.5 mm
hist(main$relMeanRMSMotion)
main <- main %>% filter(.,  relMeanRMSMotion  < 0.5)
dim(main) #2373->2367, n=6 runs removed, n=0 people removed
length(unique(main$id0))

#remove runs with over 50% of volumes censored
main$percVolCensored=main$nVolCensored/main$nVols
main <- main %>% filter(.,  percVolCensored  < 0.5 | is.na(percVolCensored))
#main <- main %>% filter(.,  percVolCensored  < 0.2 | is.na(percVolCensored)) #if we use 0.2 instead what do we get? 1323 runs, 540 people
dim(main) #2564->2373, n=191 runs removed, n=8 people removed entirely
summary(main)
length(unique(main$id0))
hist(main$percVolCensored)

#how much data do most people have left?
main$nVolsRemain <- main$nVols-main$nVolCensored
summary(main$nVolsRemain)

#and with less than two runs left now?
mainfilt <- main %>% group_by(id0) %>%
  filter(n() >= 2)
examinemeanvols <- mainfilt %>% group_by(id0) %>% filter(row_number() == 1 | row_number() ==2) %>% mutate(nVols_uncensored_total=sum(nVolsRemain))
summary(examinemeanvols$nVols_uncensored_total)
dim(mainfilt) #2367->2349, n=18 runs removed, n=18 people removed
length(unique(mainfilt$id0)) #if we use 0.2 instead we end up with only 429 people in the original sample

#save out filtered subject list
subjlist <- select(mainfilt, id0:id1) %>% ungroup() 
subjlist$id0 <- as.character(subjlist$id0)
write.csv(subjlist, paste0(subjlist_vars,"n",length(unique(subjlist$id0)),"_filtered_runs_site16_postprocess_0.2perc_outlier_cutoff.csv"))

# Save a list of motion files for each of the two runs used ---------------
library(dplyr)
subjlist <- read.csv("/cbica/projects/spatial_topography/data/subjLists/release2/site16/n670_filtered_runs_site16_postprocess.csv")
subjlist <-  select(subjlist,id:var4)

#truncate motion files to 370 lines using bash
#Run in Bash
"cd /cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols/;
for file in `find . -iname '*_tmask.1D' | cut -d. -f2`;
do
echo ${file}
head -n370 ${file:1}.1D > ${file:1}_truncate.1D
done"

subjlist <- subjlist %>% mutate(., name1=paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols/",id,"/",var1,"/confound2/mc/",id,"_",var1, "_tmask_truncate.1D"))
subjlist <- subjlist %>% mutate(., name2=paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols/",id,"/",var2,"/confound2/mc/",id,"_",var2, "_tmask_truncate.1D"))
subjlist <- subjlist %>% mutate(., name3=ifelse(is.na(var3), NA,paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols/",id,"/",var3,"/confound2/mc/",id,"_",var3, "_tmask_truncate.1D")))
subjlist <- subjlist %>% mutate(., name4=ifelse(is.na(var4), NA, paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols/",id,"/",var4,"/confound2/mc/",id,"_",var4, "_tmask_truncate.1D")))

#write it out somewhere 
write.csv(subjlist, "/cbica/projects/spatial_topography/data/subjLists/release2/site16/n670_site16_motion_files.csv")

# Examine inclusion-exclusion criteria for FD, reviewers ------------------
library(ggplot2)
main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS.csv"))
summary(main) #17 NAs, n=696, 2589 runs
length(unique(main$id0))

#Filtering out subjects based on:
#any outliers in coverage metrics (Dice, Jaccard), n=2 runs
main <- main %>% filter(., coregJaccard > 0.95) 
length(unique(main$id0)) #n=0 people, n=2 runs for bad coverage, but also 2589->2578, n=9 runs that had NAs because didn't finish removed

#total number of vols in ABCD to check for runs that are too short
main <- main %>% filter(.,  nVols  > 360 | is.na(nVols)) #I've replaced all the missing nVols from xcpEngine bug so that it reflects the correct number, people with NAs didn't finish.
summary(main)
length(unique(main$id0)) # n= 0 people
dim(main) #2571->2564 n=7 runs removed

#and any individual run with motion over 0.5 mm
hist(main$relMeanRMSMotion)
main <- main %>% filter(.,  relMeanRMSMotion  < 0.5)
dim(main) #2373->2367, n=6 runs removed, n=0 people removed
length(unique(main$id0))

#For % volumes censored in 0-0.5 in 0.01 increments, save out the number of unique subjects we'd have that meet criteria at each
main$percVolCensored=main$nVolCensored/main$nVols
data <- list()
percVolCensored <- vector()
number_of_part <- vector()
for (i in seq(0, 0.5, 0.01))
{
  percVolCensored <- c(percVolCensored,i)
  mainfilt <- main %>% filter(.,  percVolCensored  < i | is.na(percVolCensored))
  mainfilt <- mainfilt %>% group_by(id0) %>% 
    filter(n() >= 2)
  number_of_part <- c(number_of_part,length(unique(mainfilt$id0)))
}
data <-data.frame(percVolCensored, number_of_part)

ggplot(data=data, aes(x=percVolCensored, y=number_of_part))+geom_bar(stat="identity", fill="steelblue")+
  scale_x_reverse()+theme_classic()+ggtitle("# of participants included at different FD thresholds") +
  xlab("Number of participants in sample")+ylab("% threshold for excluding a run \n based on % of volumes flagged for motion > 0.2 mm FD")#+
  geom_hline(yintercept = 442)+geom_vline(xintercept = 0.2)

# Replication sample: Make subject lists ------------------------------------------------------
library(dplyr)
xcp_qa_vars="/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/"
subjlist_vars ="/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/"

main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS_with_nVols.csv")) 
#main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS.csv"))  
summary(main) #9 NAs, n=594, 2589 runs
length(unique(main$id0))

#Filtering out subjects based on:
#any outliers in coverage metrics (Dice, Jaccard), n=2 runs
main <- main %>% filter(., coregJaccard > 0.95) 
length(unique(main$id0)) 

#total number of vols in ABCD to check for runs that are too short
main <- main %>% filter(.,  nVols  > 360 )
length(unique(main$id0)) 
dim(main)

#remove subjects with mincost > 0.6 or the subject with the missing parcel for bad coverage (sub-NDARINV8GP2PFEE)
main <- main %>% filter(., !id0 %in% c("sub-NDARINVPUGJUJZ2", "sub-NDARINVW2JAJB5C","sub-NDARINV8GP2PFEE"))
length(unique(main$id0)) # n= 3 people
dim(main) 

#remove runs with over 50% of volumes censored
main$percVolCensored=main$nVolCensored/main$nVols
main <- main %>% filter(.,  percVolCensored  < 0.5 | is.na(percVolCensored))
#main <- main %>% filter(.,  percVolCensored  < 0.2 | is.na(percVolCensored)) #if we use 0.2 instead what do we get? 1323 runs, 540 people
dim(main)
summary(main)
length(unique(main$id0))

#and any individual run with motion over 0.5 mm
hist(main$relMeanRMSMotion)
main <- main %>% filter(.,  relMeanRMSMotion  < 0.5)
dim(main)
length(unique(main$id0))

#how much data do most people have left?
main$nVolsRemain <- main$nVols-main$nVolCensored
summary(main$nVolsRemain)

#and with less than two runs left now?
mainfilt <- main %>% group_by(id0) %>% 
  filter(n() >= 2)
dim(mainfilt) #n=544
length(unique(mainfilt$id0)) #if we use 0.2 instead we end up with only 429 people in the original sample

#save out filtered subject list
subjlist <- select(mainfilt, id0:id1) %>% ungroup() 
subjlist$id0 <- as.character(subjlist$id0)
write.csv(subjlist, paste0(subjlist_vars,"n",length(unique(subjlist$id0)),"_filtered_runs_site16_postprocess_0.2perc_outlier_cutoff.csv"))

# Replication sample: Examine inclusion-exclusion criteria for FD, reviewers ------------------
xcp_qa_vars="/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/"
subjlist_vars ="/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/"

main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS_with_nVols.csv")) 
summary(main) #9 NAs, n=594, 2589 runs
length(unique(main$id0))

#Filtering out subjects based on:
#any outliers in coverage metrics (Dice, Jaccard), n=2 runs
main <- main %>% filter(., coregJaccard > 0.95) 
length(unique(main$id0)) 

#total number of vols in ABCD to check for runs that are too short
main <- main %>% filter(.,  nVols  > 360 )
length(unique(main$id0)) 
dim(main)

#remove subjects with mincost > 0.6 or the subject with the missing parcel for bad coverage (sub-NDARINV8GP2PFEE)
main <- main %>% filter(., !id0 %in% c("sub-NDARINVPUGJUJZ2", "sub-NDARINVW2JAJB5C","sub-NDARINV8GP2PFEE"))
length(unique(main$id0)) # n= 3 people
dim(main)

#For % volumes censored in 0-0.5 in 0.01 increments, save out the number of unique subjects we'd have that meet criteria at each
main$percVolCensored=main$nVolCensored/main$nVols
data <- list()
percVolCensored <- vector()
number_of_part <- vector()
for (i in seq(0, 0.5, 0.01))
{
  percVolCensored <- c(percVolCensored,i)
  mainfilt <- main %>% filter(.,  percVolCensored  < i | is.na(percVolCensored))
  mainfilt <- mainfilt %>% group_by(id0) %>% 
    filter(n() >= 2)
  number_of_part <- c(number_of_part,length(unique(mainfilt$id0)))
}
data <-data.frame(percVolCensored, number_of_part)

ggplot(data=data, aes(x=percVolCensored, y=number_of_part))+geom_bar(stat="identity", fill="steelblue")+
  scale_x_reverse()+theme_classic()+ggtitle("# of participants included at different FD thresholds") +
  xlab("Number of participants in sample")+ylab("% threshold for excluding a run \n based on % of volumes flagged for motion > 0.2 mm FD")#+
geom_hline(yintercept = 313)+geom_vline(xintercept = 0.2)
