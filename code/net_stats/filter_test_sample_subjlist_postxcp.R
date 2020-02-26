# Make subject lists ------------------------------------------------------
library(dplyr)
xcp_qa_vars="~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/"
subjlist_vars ="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20/"

main <- read.csv(paste0(xcp_qa_vars,"XCP_QAVARS_with_nVols.csv"))
summary(main) #10 NAs, n=595, 2348 runs
length(unique(main$id0))

#Filtering out subjects based on:
#any outliers in coverage metrics (Dice, Jaccard), n=2 runs
main <- main %>% filter(., coregJaccard > 0.95)
length(unique(main$id0)) #n=0 people, n=2 people total for bad coverage, but also 2348->2332, n=9 runs for bad coverage
# n=2330 runs, n=1 person not finishing, n=592 after those that had NAs because didn't finish removed

#total number of vols in ABCD to check for runs that are too short
main <- main %>% filter(.,  nVols  > 360) #NA people had no volumes censored
summary(main)
length(unique(main$id0)) # n= 0 people
dim(main) #2330->2298 n=32 runs removed

#remove people with over 50% of volumes censored
main$percVolCensored=main$nVolCensored/main$nVols
main <- main %>% filter(.,  percVolCensored  < 0.5 | is.na(percVolCensored))
dim(main) #2298->2001, n=297 runs removed, n=27 people removed
summary(main)
length(unique(main$id0)) #n=565

#and any individual run motion over 0.5 mm
hist(main$relMeanRMSMotion)
main <- main %>% filter(.,  relMeanRMSMotion  < 0.5)
dim(main) #2001>2001, n=0 runs removed, n=0 people removed
length(unique(main$id0))

#how much data do most people have left?
main$nVolsRemain <- main$nVols-main$nVolCensored
summary(main$nVolsRemain)

#and with less than two runs left now?
mainfilt <- main %>% group_by(id0) %>% 
  filter(n() >= 2)
dim(mainfilt) #2001->1983, n=xx runs removed, n=18 people removed
length(unique(mainfilt$id0))

#save out filtered subject list
subjlist <- select(mainfilt, id0:id1) %>% ungroup() 
subjlist$id0 <- as.character(subjlist$id0)
write.csv(subjlist, paste0(subjlist_vars,"n",length(unique(subjlist$id0)),"_filtered_runs_site14site20_postprocess.csv"))

####
library(data.table)
#save out the subject list for gwMRF function with all runs on one line
#add an index for how many runs there might be w/in a subject
subjlist <- subjlist %>%
  group_by(id0) %>%
  mutate(number = 1:n())
subjlist <- as.data.frame(subjlist)
subjlist <- dcast(subjlist, id0 ~ number, value.var="id1") #reshape it to wide
names(subjlist) <- c("id", "var1","var2","var3", "var4") #rename columns so not numbers
subjlist <- subjlist %>% mutate(., name1=paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var1, "_task-rest_residualised_reshape.nii.gz"))
subjlist <- subjlist %>% mutate(., name2=paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var2, "_task-rest_residualised_reshape.nii.gz"))
subjlist <- subjlist %>% mutate(., name3=ifelse(is.na(var3), NA,paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var3, "_task-rest_residualised_reshape.nii.gz")))
subjlist <- subjlist %>% mutate(., name4=ifelse(is.na(var4), NA, paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var4, "_task-rest_residualised_reshape.nii.gz")))
subjlist <- subjlist %>% mutate(., name5=ifelse(is.na(var5), NA, paste0("/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/surfaces/",id,"/surf/lh.fs6_sm4.8_",id,"_",var5, "_task-rest_residualised_reshape.nii.gz")))

#write it out somewhere 
write.csv(subjlist, paste0(subjlist_vars,"n",length(unique(subjlist$id)),"_filtered_runs_site14site20_postprocess.csv"))


# Make histograms of registration coverage ---------------------------------------
xcp_qa_vars="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/"

main <- read.csv(paste0(xcp_qa_vars,"n696_site16_filt_XCP_QAVARS.csv"))
main2 <- main %>% filter(., coregJaccard > 0.95)

hist(main2$coregJaccard,  col="blue", main="Coreg Jaccard")

hist(main2$coregDice,  col="lightblue", main="Coreg Dice")

hist(main2$coregCrossCorr,  col="darkblue", main="Coreg CrossCorr")

