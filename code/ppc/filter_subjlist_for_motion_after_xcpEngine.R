# Make subject lists ------------------------------------------------------
library(dplyr)
xcp_qa_vars="~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols/"
subjlist_vars ="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/"

main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS_with_nVols.csv")) 
summary(main) #17 NAs, n=696, 2589 runs
length(unique(main$id0))

#Filtering out subjects based on:
#any outliers in coverage metrics (Dice, Jaccard), n=2 runs
main <- main %>% filter(., coregJaccard > 0.95) 
length(unique(main$id0)) #n=0 people, n=2 runs for bad coverage, but also 2589->2578, n=9 runs that had NAs because didn't finish removed

#total number of vols in ABCD to check for runs that are too short
main <- main %>% filter(.,  nVols  > 360 ) #I've replaced all the missing nVols from xcpEngine bug so that it reflects the correct number, people with NAs didn't finish.
summary(main)
length(unique(main$id0)) # n= 0 people
dim(main) #2571->2564 n=7 runs removed

#remove runs with over 50% of volumes censored
main$percVolCensored=main$nVolCensored/main$nVols
main <- main %>% filter(.,  percVolCensored  < 0.5 | is.na(percVolCensored))
dim(main) #2564->2373, n=191 runs removed, n=8 people removed entirely
summary(main)
length(unique(main$id0))

#and any individual run with motion over 0.5 mm
hist(main$relMeanRMSMotion)
main <- main %>% filter(.,  relMeanRMSMotion  < 0.5)
dim(main) #2373->2367, n=6 runs removed, n=0 people removed
length(unique(main$id0))

#how much data do most people have left?
main$nVolsRemain <- main$nVols-main$nVolCensored
summary(main$nVolsRemain)

#and with less than two runs left now?
mainfilt <- main %>% group_by(id0) %>% 
  filter(n() >= 2)
dim(mainfilt) #2367->2349, n=18 runs removed, n=18 people removed
length(unique(mainfilt$id0))

#save out filtered subject list
subjlist <- select(mainfilt, id0:id1) %>% ungroup() 
subjlist$id0 <- as.character(subjlist$id0)
write.csv(subjlist, paste0(subjlist_vars,"n",length(unique(subjlist$id0)),"_filtered_runs_site16_postprocess.csv"))

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

