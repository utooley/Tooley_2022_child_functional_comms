# Packages ----------------------------------------------------------------
#turn off scientific notation
options(scipen = 9)

library(fsbrain) #this may not work if you're editing the script directly on the cluster...
library(freesurferformats)
library(dplyr)
library(R.matlab)
library(stringr)
library(ggplot2)
library(tidyr)
library(igraph)
library(aricode)
library(reshape)

#rearrange the order of the brains in the T9 view of fsbrain
source("~/Documents/tools/fsbrain_fix_t9.R")
environment(brainview.t9) <- asNamespace('fsbrain')
assignInNamespace("brainview.t9", brainview.t9, ns = "fsbrain")

# SETUP -------------------------------------------------------------------
#CUBIC cluster filepaths, mounted locally
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
yeo_dev_dir="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_mot_outliers/"
yeo7_ref_dir="/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/"

# Read in files -----------------------------------------------------------
#Vector of WSBM community assignments
partitions <- readMat(paste0(wsbm_datadir,"n670_training_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled

#Read in the frequency of ties for consensus WSBM partition
consensus <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
freq<- t(consensus$freq) #frequency of ties in the consensus representative partition

#Read in the Yeo developmental vector of assignments (same with and w/o motion outliers)
yeo_dev_partition <- readMat(paste0(yeo_dev_dir,"yeo7_n670_2runsonly_1000tries_mot_outliers.mat"), drop = )

#make yeo colors
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53")) #these are correct

# Set defaults ------------------------------------------------------------
subject_id = 'fsaverage';       # for functions which use one subject only
atlas='Schaefer2018_400Parcels_7Networks_order'

#Get the Schaefer atlas region names
schaefer_atlas_region_names_lh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");

#Change the output directory for figures to be on the cluster brains
setwd("/cbica/projects/spatial_topography/output/images/brains/")

#Set RGL defaults, larger window and how to make snapshots!
rgloptions=list("windowRect"=c(50,50,1000,1000));

# Versatility  in the WSBM-Supp Fig S6a -------------------------------------
#load silhouette data
z <- readMat(paste0(wsbm_datadir,"versatility_n670.mat"), drop = )
versatility_wsbm <- z$V

silhouette_wsbm_lh=as.list(setNames(c(NA, versatility_wsbm[1:200]), schaefer_atlas_region_names_lh))
silhouette_wsbm_rh=as.list(setNames(c(NA, versatility_wsbm[201:400]), schaefer_atlas_region_names_rh))
silhouette_wsbm_lh

rgloptions=list("windowRect"=c(50,50,1000,1000));
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
rglactions=list("snapshot_png"=paste0(output_image_directory,"versatility_wsbm.png"))
colFn_diverging = colorRampPalette(rev(c("burlywood3","white", "white")));
makecmap_options=list('colFn'=colFn_diverging)

vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  silhouette_wsbm_lh,
                             silhouette_wsbm_rh, makecmap_options = makecmap_options, "inflated", views="t4", rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

#Compare to uncertainty in WSBM using z-score of Rand 
freq=abs(freq-670)
freq=freq/max(freq)
clustComp(c(versatility_wsbm), c(freq))

b <- readMat(paste0(wsbm_datadir,"test_freq.mat"), drop = )
freq2 <- b$freq

# Partitions with k= 17 ---------------------------------------------------
wsbm_dir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains"
#Vector of WSBM community assignments
partitions <- readMat(paste0(wsbm_datadir,"n670_k17training_sample_consensus_partition_wsbmrelabeled.mat"), drop = )
consensus_k17_labels <- partitions$consensus.iter.mode

unique(consensus_k17_labels) #there are only 8 communities

wsbm_consensus_lh=as.list(setNames(c(0, consensus_k17_labels[1:200]), schaefer_atlas_region_names_lh))
wsbm_consensus_rh=as.list(setNames(c(0, consensus_k17_labels[201:400]), schaefer_atlas_region_names_rh))

#make colormap of Yeo colors
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53","#55EBD1")) #added on end
makecmap_options=list('colFn'=yeo_colors)
rgloptions=list("windowRect"=c(50,50,1000,1000));
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus//"
rglactions=list("snapshot_png"=paste0(output_image_directory,"17_communities_wsbm_relabeled.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, makecmap_options = makecmap_options, "inflated", views="t4",rgloptions = rgloptions, rglactions = rglactions)

#Read in the Yeo k=17 partition
yeo_dev_partition <- readMat("/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/search_over_k/Cluster017.s00.tries100.rand100.znorm1.dim1..mat", drop = )
yeo_dev_labels_17 <- yeo_dev_partition$con.struct[[1]][[1]][[1]][[2]]
#somehow figure out how to get left and right out of that?
yeo_dev_lh <- yeo_dev_labels_17
yeo_dev_rh <- yeo_dev_partition$rh.labels

#make yeo colors
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53")) #these are correct
rgloptions=list("windowRect"=c(50,50,1000,1000));
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"
rglactions=list("snapshot_png"=paste0(output_image_directory,"17_communities.png"))
makecmap_options=list('colFn'=yeo_colors)
vis.data.on.subject(subjects_dir, 'fsaverage6',yeo_dev_lh, NULL, "inflated", makecmap_options = makecmap_options, views="t4", rgloptions = rgloptions)

