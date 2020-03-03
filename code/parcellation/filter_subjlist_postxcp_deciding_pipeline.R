qa_dir_one="~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_36p_gsr_multruns_scrub_dropvols/"
qa_dir_two="~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.5mm_dropvols/"

allparams=read.csv(paste0(qa_dir_one, "XCP_QAVARS.csv"))
marek=read.csv(paste0(qa_dir_two, "XCP_QAVARS.csv"))

allparams$diff=allparams$motionDVCorrInit-allparams$motionDVCorrFinal
marek$diff=marek$motionDVCorrInit-marek$motionDVCorrFinal

hist(allparams$diff, col="aquamarine",main=paste0("Mean diff b/w init FD-DVARs corr is ", round(mean(allparams$diff), digits = 2)))
hist(marek$diff, col="aquamarine3",main=paste0("Mean diff b/w init FD-DVARs corr is ", round(mean(marek$diff), digits = 2)))


# Make subject lists ------------------------------------------------------
library(dplyr)
xcp_qa_vars="~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/"
subjlist_vars ="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/"

main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS_with_nVols.csv")) ## NEED TO RECREATE SUBJECT LISTS NOW WITH FILTERING OUT SHORT SCANS--ANYONE WHOSE SHORT
#SCAN WAS INCORPORATED NEEDS A NEW AVERAGED SUBJECT MATRIX
summary(main) #17 NAs, n=696, 2589 runs
length(unique(main$id0))

#Filtering out subjects based on:
#any outliers in coverage metrics (Dice, Jaccard), n=2 runs
main <- main %>% filter(., coregJaccard > 0.95)
length(unique(main$id0)) #n=0 people, n=2 runs for bad coverage, but also 2589->2578, n=9 runs that had NAs because didn't finish removed

#total number of vols in ABCD to check for runs that are too short
main <- main %>% filter(.,  nVols  > 360 ) #I've replaced all the nVols so that it reflects the correct number, people with NAs didn't finish.
summary(main)
length(unique(main$id0)) # n= 0 people
dim(main) #2571->2564 n=xx runs removed

#remove people with over 50% of volumes censored
main$percVolCensored=main$nVolCensored/main$nVols
main <- main %>% filter(.,  percVolCensored  < 0.5 | is.na(percVolCensored))
dim(main) #2564->2373, n=xx runs removed, n=8 people removed
summary(main)
length(unique(main$id0))

#and any individual run motion over 0.5 mm
hist(main$relMeanRMSMotion)
main <- main %>% filter(.,  relMeanRMSMotion  < 0.5)
dim(main) #2373->2367, n=xx runs removed, n=0 people removed
length(unique(main$id0))

#how much data do most people have left?
main$nVolsRemain <- main$nVols-main$nVolCensored
summary(main$nVolsRemain)

#and with less than two runs left now?
mainfilt <- main %>% group_by(id0) %>% 
  filter(n() >= 2)
dim(mainfilt) #2367->2349, n=xx runs removed, n=18 people removed
length(unique(mainfilt$id0))

#save out filtered subject list
subjlist <- select(mainfilt, id0:id1) %>% ungroup() 
subjlist$id0 <- as.character(subjlist$id0)
write.csv(subjlist, paste0(subjlist_vars,"n",length(unique(subjlist$id0)),"_filtered_runs_site16_postprocess.csv"))

####
#save out the subject list for gwMRF function with all runs on one line
#add an index for how many runs there might be w/in a subject
subjlist <- subjlist %>%
  group_by(id0) %>%
  mutate(number = 1:n())
subjlist <- as.data.frame(subjlist)
subjlist <- dcast(subjlist, id0 ~ number, value.var="id1") #reshape it to wide
names(subjlist) <- c("id", "var1","var2","var3", "var4") #rename columns so not numbers
subjlist <- subjlist %>% mutate(., name1=paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var1, "_task-rest_residualised_reshape.nii.gz"))
subjlist <- subjlist %>% mutate(., name2=paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var2, "_task-rest_residualised_reshape.nii.gz"))
subjlist <- subjlist %>% mutate(., name3=ifelse(is.na(var3), NA,paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var3, "_task-rest_residualised_reshape.nii.gz")))
subjlist <- subjlist %>% mutate(., name4=ifelse(is.na(var4), NA, paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var4, "_task-rest_residualised_reshape.nii.gz")))

#write it out somewhere 
write.csv(subjlist, paste0(subjlist_vars,"parcellation/n",length(unique(subjlist$id)),"_filtered_runs_site16_postprocess.csv"))


# Save a list of motion files for each of the two runs used ---------------
library(dplyr)
subjlist <- read.csv("/cbica/projects/spatial_topography/data/subjLists/release2/site16/parcellation/n670_filtered_runs_site16_postprocess.csv")
subjlist <-  select(subjlist,id:var4)
#truncate motion files to 370 lines using bash
#Run in BASH
"cd /cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/;
for file in `find . -iname '*_tmask.1D' | cut -d. -f2`;
do
echo ${file}
head -n370 ${file:1}.1D > ${file:1}_truncate.1D
done"

subjlist <- subjlist %>% mutate(., name1=paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/",id,"/",var1,"/confound2/mc/",id,"_",var1, "_tmask_truncate.1D"))
subjlist <- subjlist %>% mutate(., name2=paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/",id,"/",var2,"/confound2/mc/",id,"_",var2, "_tmask_truncate.1D"))
subjlist <- subjlist %>% mutate(., name3=ifelse(is.na(var3), NA,paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/",id,"/",var3,"/confound2/mc/",id,"_",var3, "_tmask_truncate.1D")))
subjlist <- subjlist %>% mutate(., name4=ifelse(is.na(var4), NA, paste0("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/",id,"/",var4,"/confound2/mc/",id,"_",var4, "_tmask_truncate.1D")))

#write it out somewhere 
write.csv(subjlist, "/cbica/projects/spatial_topography/data/subjLists/release2/site16/n670_site16_motion_files.csv")

# Make histograms of registration coverage ---------------------------------------
xcp_qa_vars="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/"

main <- read.csv(paste0(xcp_qa_vars,"n696_site16_filt_XCP_QAVARS.csv"))
main2 <- main %>% filter(., coregJaccard > 0.95)

hist(main2$coregJaccard,  col="blue", main="Coreg Jaccard")

hist(main2$coregDice,  col="lightblue", main="Coreg Dice")

hist(main2$coregCrossCorr,  col="darkblue", main="Coreg CrossCorr")

