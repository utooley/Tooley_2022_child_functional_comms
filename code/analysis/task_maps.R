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

#Take the bottom 20% of vertices for each map
nback_2_vs_0_bin <- rep(0, length(nback_2_vs_0_perf))
nback_2_vs_0_bin[nback_2_vs_0_perf<quantile(nback_2_vs_0_perf, probs = 0.20)] <- 1

mid_largewin_vs_neut_bin <- rep(0, length(mid_largewin_vs_neut))
mid_largewin_vs_neut_bin[mid_largewin_vs_neut<quantile(mid_largewin_vs_neut , probs = 0.20)] <- 1

sst_crct_stop_vs_failed_bin <- rep(0, length(sst_crct_stop_vs_failed))
sst_crct_stop_vs_failed_bin[sst_crct_stop_vs_corr_go>quantile(sst_crct_stop_vs_corr_go, probs = 0.80)] <- 1

#Overlay where they overlap.
dmn_min <- nback_2_vs_0_bin+mid_largewin_vs_neut_bin+sst_crct_stop_vs_failed_bin
hist(dmn_min)

#Look at this overlap
overlap_colors=colorRampPalette(c("lightgray", "lightpink", "red", "darkred"))
makecmap_options=list('colFn'=overlap_colors, 'n'=4)
vis.data.on.subject(subjects_dir, 'fsaverage6',dmn_min[0:40962], dmn_min[40963:81924], "inflated", views="t4", makecmap_options = makecmap_options,
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
