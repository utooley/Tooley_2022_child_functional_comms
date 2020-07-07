# Setup -------------------------------------------------------------------
options(rgl.useNULL=TRUE) #if fsbrain not working!
.rs.restartR()

options(scipen = 9)

library(fsbrain) #this may not work if you're editing the script directly on the cluster...
library(freesurferformats)
library(dplyr)
library(R.matlab)
library(stringr)
library(ggplot2)
library(tidyr)
library(data.table)

#rearrange the order of the brains in the T9 view of fsbrain
source("~/Documents/tools/fsbrain_fix_t9.R")
environment(brainview.t9) <- asNamespace('fsbrain')
assignInNamespace("brainview.t9", brainview.t9, ns = "fsbrain")

# Paths -------------------------------------------------------------------
#CUBIC cluster filepaths, mounted locally
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
yeo_dev_dir="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_mot_outliers/"
yeo7_ref_dir="/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/"
task_maps_dir="/cbica/projects/spatial_topography/data/imageData/task_maps/"

# Read in files -----------------------------------------------------------
#Vector of WSBM community assignments
partitions <- readMat(paste0(wsbm_datadir,"n670_training_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_representative_labels <-  partitions$consensus.represent.yeorelabeled
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled

#Read in entropy of co-occurence matrix for representative partition and the frequency of ties
consensus <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
entropy<- consensus$node.entropy #entropy 
freq<- t(consensus$freq) #frequency of ties in the consensus representative partition

#Read in the Yeo developmental vector of assignments
yeo_dev_partition <- readMat(paste0(yeo_dev_dir,"yeo7_n670_2runsonly_1000tries_mot_outliers.mat"), drop = )
#make yeo colors
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53")) #these are corr

# Set defaults ------------------------------------------------------------
subject_id = 'fsaverage';       # for function which use one subject only
atlas='Schaefer2018_400Parcels_7Networks_order'

#Get the Schaefer atlas region names
schaefer_atlas_region_names_lh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");

#Change the output directory for figures to be on the cluster brains
setwd("/cbica/projects/spatial_topography/output/images/brains/")

#Set RGL defaults, larger window and how to make snapshots!
rgloptions=list("windowRect"=c(50,50,1000,1000));

#FIGURE OUT HOW TO ROTATE ON A DIFFERENT AXIS FOR RGL
#This edits the hidden function vis.coloredmeshes.rotating, figure out a way to edit this permanently. 
#Change x=0, y=0, z=1, maybe rpm
#trace(fsbrain:::brainview.sr, edit=TRUE)

# Read in task data -------------------------------------------------------
#Read in the task data
#FIRST NEED TO CONVERT TO FSAVERAGE OUTSIDE OF R WITH THE SCRIPT transform_task_maps_to_fsaverage.sh, then read in
lh_nback <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.nBack_2_back_vs_0_back_performance_2back_.lh.41k_fsavg_L.func.gii"), element_index = 1L)
rh_nback <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.nBack_2_back_vs_0_back_performance_2back_.rh.41k_fsavg_R.func.gii"), element_index = 1L)
nback_2_vs_0_perf <- c(lh_nback,rh_nback)

lh_nback_faces_vs_places <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.nBack_face_vs_place_performance_2back_.lh.41k_fsavg_L.func.gii"), element_index = 1L)
rh_nback_faces_vs_places <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.nBack_face_vs_place_performance_2back_.rh.41k_fsavg_R.func.gii"), element_index = 1L)
nback_faces_vs_places <- c(lh_nback_faces_vs_places,rh_nback_faces_vs_places)

#Sanity check of plotting
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list('shift_hemis_apart'=list('min_dist'=100))
colFn_diverging = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, name="RdBu")));
makecmap_options=list('colFn'=colFn_diverging, 'n'=100, 'symm'= TRUE)

vis.data.on.subject(subjects_dir, 'fsaverage6',lh_nback_faces_vs_places, rh_nback_faces_vs_places, "inflated", views="t9", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = TRUE)

lh_mid_largeloss <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_large_loss_vs_neutral_.lh.41k_fsavg_L.func.gii"))
rh_mid_largeloss <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_large_loss_vs_neutral_.rh.41k_fsavg_R.func.gii"))
mid_largeloss_vs_neut <- c(lh_mid_largeloss,rh_mid_largeloss)

lh_mid_largewin <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_large_reward_vs_neutral_.lh.41k_fsavg_L.func.gii"))
rh_mid_largewin <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_large_reward_vs_neutral_.rh.41k_fsavg_R.func.gii"))
mid_largewin_vs_neut <- c(lh_mid_largewin,rh_mid_largewin)

lh_mid_ant_win <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_reward_vs_neutral_.lh.41k_fsavg_L.func.gii"))
rh_mid_ant_win <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_reward_vs_neutral_.rh.41k_fsavg_R.func.gii"))
mid_ant_win <- c(lh_mid_ant_win,rh_mid_ant_win)

lh_mid_ant_loss <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_loss_vs_neutral_.lh.41k_fsavg_L.func.gii"))
rh_mid_ant_loss <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_antic_loss_vs_neutral_.rh.41k_fsavg_R.func.gii"))
mid_ant_loss <- c(lh_mid_ant_loss,rh_mid_ant_loss)

lh_mid_feedback_win_posneg <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_reward_pos_vs_neg_feedback_.lh.41k_fsavg_L.func.gii"))
rh_mid_feedback_win_posneg <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_reward_pos_vs_neg_feedback_.rh.41k_fsavg_R.func.gii"))
mid_feedback_win_pos_vs_neg <- c(lh_mid_feedback_win_posneg,rh_mid_feedback_win_posneg)

lh_mid_feedback_loss_posneg <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_loss_pos_vs_neg_feedback_.lh.41k_fsavg_L.func.gii"))
rh_mid_feedback_loss_posneg <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.MID_loss_pos_vs_neg_feedback_.rh.41k_fsavg_R.func.gii"))
#mid_feedback_loss_pos_vs_neg <- c(lh_mid_feedback_win_posneg,rh_mid_feedback_win_posneg)

lh_sst_corr_stop_vs_failed <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.SST_correct_stop_vs_incorrect_stop_performance_.lh.41k_fsavg_L.func.gii"))
rh_sst_corr_stop_vs_failed <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.SST_correct_stop_vs_incorrect_stop_performance_.rh.41k_fsavg_R.func.gii"))
sst_crct_stop_vs_failed <- c(lh_sst_corr_stop_vs_failed,rh_sst_corr_stop_vs_failed)

lh_sst_corr_stop_vs_corr_go <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.SST_correct_stop_vs_correct_go_performance_.lh.41k_fsavg_L.func.gii"))
rh_sst_corr_stop_vs_corr_go <- read.fs.morph.gii(paste0(task_maps_dir,"ABCD.SST_correct_stop_vs_correct_go_performance_.rh.41k_fsavg_R.func.gii"))
sst_crct_stop_vs_corr_go <- c(lh_sst_corr_stop_vs_corr_go,rh_sst_corr_stop_vs_corr_go)

#Sanity check of plotting
rgloptions=list("windowRect"=c(50,50,1000,1000));
colFn_diverging = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, name="RdBu")));
makecmap_options=list('colFn'=colFn_diverging, 'n'=100, 'symm'=TRUE)
vis.data.on.subject(subjects_dir, 'fsaverage6', lh_mid_feedback_loss_posneg, rh_mid_feedback_loss_posneg, "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)

# Make a default system conjunction map -----------------------------------
# MID_antic_loss_vs_neutral_dat_cohen_c1.dscalar.nii, MID_antic_reward_vs_neutral_dat_cohen_c1.dscalar.nii for MID and then overlap with 2 vs 0 nback.
#Take the bottom 20% of vertices for each map
nback_2_vs_0_bin <- rep(0, length(nback_2_vs_0_perf))
nback_2_vs_0_bin[nback_2_vs_0_perf<quantile(nback_2_vs_0_perf, probs = 0.20)] <- 1

mid_ant_loss_bin <- rep(0, length(mid_ant_loss))
mid_ant_loss_bin[mid_ant_loss<quantile(mid_ant_loss , probs = 0.20)] <- 1

mid_ant_win_bin <- rep(0, length(mid_ant_win))
mid_ant_win_bin[mid_ant_win<quantile(mid_ant_win , probs = 0.20)] <- 1

#Overlay MID two maps
mid_ant_loss_win_overlap=mid_ant_loss_bin+mid_ant_win_bin
hist(mid_ant_loss_win_overlap)
overlap_colors=colorRampPalette(c("lightgray", "lightpink", "red"))
makecmap_options=list('colFn'=overlap_colors, 'n'=3)
vis.data.on.subject(subjects_dir, 'fsaverage6',mid_ant_loss_win_overlap[0:40962], mid_ant_loss_win_overlap[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)

#Try with continuous MID data and then bottom 20% of the summed data
mid_ant_loss_win_cont <- scale(mid_ant_win)+scale(mid_ant_loss)
mid_ant_loss_win_cont_bin<- rep(0, length(mid_ant_loss_win_cont))
mid_ant_loss_win_cont_bin[mid_ant_loss_win_cont<quantile(mid_ant_loss_win_cont, probs = 0.20)] <- 1
dmn_min <- nback_2_vs_0_bin+mid_ant_loss_win_cont_bin

#where they overlap.
dmn_min <- nback_2_vs_0_bin+mid_ant_loss_win_overlap
hist(dmn_min)

#Look at this overlap
overlap_colors=colorRampPalette(c("lightgray", "lightgray", "lightgray", "darkred"))
makecmap_options=list('colFn'=overlap_colors, 'n'=4)
vis.data.on.subject(subjects_dir, 'fsaverage6',dmn_min[0:40962], dmn_min[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)

# Make a frontoparietal system conjunction map -----------------------------------
#Just take top of nback 2 vs 0
nback_2_vs_0_top <- rep(0, length(nback_2_vs_0_perf))
nback_2_vs_0_top[nback_2_vs_0_perf>quantile(nback_2_vs_0_perf, probs = 0.90)] <- 1
#Look at this
overlap_colors=colorRampPalette(c("lightgray", "darkred"))
makecmap_options=list('colFn'=overlap_colors)
vis.data.on.subject(subjects_dir, 'fsaverage6',nback_2_vs_0_top[0:40962], nback_2_vs_0_top[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)


# Make a limbic system conjunction map -----------------------------------
#Take MID feedback, reward pos vs. negative.
mid_feedback_win_pos_vs_neg_bin <- rep(0, length(mid_feedback_win_pos_vs_neg))
mid_feedback_win_pos_vs_neg_bin[mid_feedback_win_pos_vs_neg > quantile(mid_feedback_win_pos_vs_neg, probs = 0.80)] <- 1
#Look at this
overlap_colors=colorRampPalette(c("lightgray", "darkred"))
makecmap_options=list('colFn'=overlap_colors)
vis.data.on.subject(subjects_dir, 'fsaverage6',mid_feedback_win_pos_vs_neg_bin[0:40962], mid_feedback_win_pos_vs_neg_bin[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)

###################################
########## COMPARE TO PARTITIONS ######
###################################
#Yeo adult
#other alternative is read in the Freesurfer version of Yeo7 in  fsaverage6 space
#this is more theoretically motivated, because this is how Yeo-dev was generated
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";
lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','Yeo2011_7Networks_N1000' )
yeo7_lh <- as.numeric(str_split(lh$label_names, "_", simplify = TRUE)[,2])
yeo7_lh[is.na(yeo7_lh)] <- 0
rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','Yeo2011_7Networks_N1000' )
yeo7_rh <- as.numeric(str_split(rh$label_names, "_", simplify = TRUE)[,2])
yeo7_rh[is.na(yeo7_rh)] <- 0
yeo7 <- c(yeo7_lh,yeo7_rh)

#Yeo dev
yeo_dev_lh <- yeo_dev_partition$lh.labels
yeo_dev_rh <- yeo_dev_partition$rh.labels
yeo_dev <- c(yeo_dev_lh,yeo_dev_rh)

#WSBM but in fsaverage6 space
#copy WSBM annotation into local CBIG subjects dir
get.atlas.region.names("wsbm.consensus.fsaverage6", template_subjects_dir = subjects_dir,template_subject='fsaverage6', hemi="rh");
wsbm_rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.consensus.fsaverage6')
wsbm_rh<- as.numeric(as.character(wsbm_rh$label_names))
wsbm_lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.consensus.fsaverage6')
wsbm_lh<- as.numeric(as.character(wsbm_lh$label_names))
wsbm <- c(wsbm_lh,wsbm_rh)

# Tally partition assignments-FP ---------------------------------------------
library(igraph);library(aricode)
partitions=cbind(yeo7, yeo_dev, wsbm)
names=c("yeo7", "yeo_dev", "wsbm")
for (i in 1:3){
  partition=partitions[,i]
  v <- partition[nback_2_vs_0_top==1]
  assign(paste0(names[i], "_FP_tally_vertices"), as.data.frame(table(v)))
}

vertices <- merge(wsbm_FP_tally_vertices, yeo_dev_FP_tally_vertices,by="v",all.x = T, suffixes = c("wsbm", "yeo_dev"))
#surf_area <- merge(surf_area, yeo7_FP_tally_vertices,by="v",all.x = T) Add in Yeo7 or not
#colnames(surf_area) <- c("community", "yeodev", "wsbm","yeo7" )
colnames(vertices) <- c("community", "yeodev", "wsbm")
longdata <- vertices %>% melt()
#sum of all diff surf areas is the same
ggplot(data=longdata, aes(x=community, y= value))+geom_bar(stat="identity", aes(fill = variable), position = "dodge")+labs(y="Number of vertices")+theme(axis.text.x=element_text(angle=90,hjust=1))+scale_fill_manual(values=yeo_colors(3))

#Compare with statistics!
wsbm_6 =ifelse(wsbm==6,1,0)
yeo_dev_6=ifelse(yeo_dev==6,2,0)
yeo7_6=ifelse(yeo7==6,1,0)#about the same as yeo_dev
c1 <- mclustcomp(wsbm_6,nback_2_vs_0_top, types = c("jaccard", "sdc")) #for binary vectors
c2 <- mclustcomp(yeo_dev_6,nback_2_vs_0_top,  types = c("jaccard", "sdc")) #for binary vectors
cbind(c1,c2)
#these seem contradictory with vertices....it's because WSBM could just have more vertices to FPN overall, but the pattern of FPN doesn't match better.

#PLOT IT
nback_2_vs_0_top[nback_2_vs_0_top==3] <- 6
total <- wsbm_6+yeo_dev_6+nback_2_vs_0_top
hist(total)
#Each color has a different WSBM=1, Yeo=2, nback=6, nback+WSBM=7, nback+yeo=8, all=9
overlap_colors=colorRampPalette(c("lightgray", "aquamarine", "cadetblue1", "cyan", "lightgray", "lightgray", "orange", "gold4", "gold2", "darkslateblue"), )
makecmap_options=list('colFn'=overlap_colors, 'n'=10)
vis.data.on.subject(subjects_dir, 'fsaverage6',total[0:40962], total[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)

# Tally partition assignments-Limbic ---------------------------------------------
library(igraph);library(aricode)
partitions=cbind(yeo7, yeo_dev, wsbm)
names=c("yeo7", "yeo_dev", "wsbm")
mid_feedback_win_pos_vs_neg_bin

#Compare with statistics!
wsbm_5 =ifelse(wsbm==5,1,0)
yeo_dev_5=ifelse(yeo_dev==5,2,0)
yeo7_5=ifelse(yeo7==5,1,0)#about the same as yeo_dev
c1 <- mclustcomp(wsbm_5,mid_feedback_win_pos_vs_neg_bin, types = c("jaccard", "sdc")) #for binary vectors
c2 <- mclustcomp(yeo_dev_5,mid_feedback_win_pos_vs_neg_bin,  types = c("jaccard", "sdc")) #for binary vectors
cbind(c1,c2)
#these seem contradictory with vertices....it's because WSBM could just have more vertices to FPN overall, but the pattern of FPN doesn't match better.

#PLOT IT
mid_feedback_win_pos_vs_neg_bin[mid_feedback_win_pos_vs_neg_bin==1] <- 6
total <- wsbm_5+yeo_dev_5+mid_feedback_win_pos_vs_neg_bin
hist(total)
#Each color has a different WSBM=1, Yeo=2, nback=6, nback+WSBM=7, nback+yeo=8, all=9
overlap_colors=colorRampPalette(c("lightgray", "aquamarine", "cadetblue1", "cyan", "lightgray", "lightgray", "orange", "gold4", "gold2", "darkslateblue"), )
makecmap_options=list('colFn'=overlap_colors, 'n'=10)
vis.data.on.subject(subjects_dir, 'fsaverage6',total[0:40962], total[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)

# Tally partition assignments-DMN ---------------------------------------------
library(igraph);library(aricode)
partitions=cbind(yeo7, yeo_dev, wsbm)
names=c("yeo7", "yeo_dev", "wsbm")
# DMN binary
#MID 20% overlap with n-back, binarized
mid_ant_loss_win_overlap=mid_ant_loss_bin+mid_ant_win_bin
dmn_min <- nback_2_vs_0_bin+mid_ant_loss_win_overlap
dmn_bin <- ifelse(dmn_min==3,1,0)

#alt DMN 
# dmn_min <- nback_2_vs_0_bin+mid_ant_loss_win_cont_bin
# dmn_bin <- ifelse(dmn_min==2,1,0)

#Compare with statistics!
wsbm_7 =ifelse(wsbm==7,1,0)
yeo_dev_7=ifelse(yeo_dev==7,2,0)
yeo7_7=ifelse(yeo7==7,1,0) #about the same as yeo_dev
c1 <- mclustcomp(wsbm_7,dmn_bin, types = c("jaccard", "sdc")) #for binary vectors
c2 <- mclustcomp(yeo_dev_7,dmn_bin,  types = c("jaccard", "sdc")) #for binary vectors
cbind(c1,c2)
#these seem contradictory with vertices....it's because WSBM could just have more vertices to FPN overall, but the pattern of FPN doesn't match better.

#PLOT IT
dmn_bin[dmn_bin==1] <- 6
total <- wsbm_7+yeo_dev_7+dmn_bin
hist(total)
#Each color has a different WSBM=1, Yeo=2, nback=6, nback+WSBM=7, nback+yeo=8, all=9
overlap_colors=colorRampPalette(c("lightgray", "aquamarine", "cadetblue1", "cyan", "lightgray", "lightgray", "orange", "gold4", "gold2", "darkslateblue"), )
makecmap_options=list('colFn'=overlap_colors, 'n'=10)
vis.data.on.subject(subjects_dir, 'fsaverage6',total[0:40962], total[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
                    rgloptions = rgloptions, draw_colorbar = TRUE)


##########################
####### PLOTTING #########
##########################

# Unneeded gifti brainstorming --------------------------------------------

#read gifti
left_gii <- readgii("/Users/utooley/Dropbox/projects/in_progress/arun_fc_metrics_motion/docs/Motion_FunctionalConnectivity/manuscript/Figure4/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii")
right_gii <- readgii("/Users/utooley/Dropbox/projects/in_progress/arun_fc_metrics_motion/docs/Motion_FunctionalConnectivity/manuscript/Figure4/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii")
# c <- surf_triangles(l)
# rgl::rgl.open()
# rgl::rgl.triangles(c$pointset)
# rgl::play3d(rgl::spin3d(axis = c(0,1, 1)), duration = 10)
#read cifti
cifti='/Users/utooley/Documents/projects/in_progress/spatial_topography_CUBIC/data/ursula/nBack_2_back_vs_0_back_performance_2back_dat_cohen_c1.dscalar.nii'
left = cifti_coords_gifti(
  cifti,
  gii_file = left_gii,
  structure = "CIFTI_STRUCTURE_CORTEX_LEFT")
right = cifti_coords_gifti(
  cifti,
  gii_file = right_gii,
  structure = "CIFTI_STRUCTURE_CORTEX_RIGHT")
c <- readcii(cifti)
