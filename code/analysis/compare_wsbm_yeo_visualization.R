library(ggseg3d)
library(ggsegExtra)
library(tidyr)
library(R.matlab)
library(dplyr)
library(fsbrain)

# SETUP -------------------------------------------------------------------
#CUBIC cluster filepaths, mounted locally
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";

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

# Set defaults ------------------------------------------------------------
subject_id = 'fsaverage';       # for function which use one subject only
atlas='Schaefer2018_400Parcels_7Networks_order'
#Change the output directory for figures to be on the cluster brains
setwd("/cbica/projects/spatial_topography/output/images/brains/")

#Set RGL defaults, larger window and how to make snapshots!
rgloptions=list("windowRect"=c(50,50,1000,1000));

# Plot Yeo7 atlas on brain ---------------------------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo7/"

#Get the Schaefer atlas region names
schaefer_atlas_region_names_lh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");

rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.subject.annot(subjects_dir, 'fsaverage', 'Schaefer2018_400Parcels_7Networks_order', 'both',  'inflated', views=c('t4'), rgloptions = rgloptions, rglactions = rglactions);

rgloptions=list("windowRect"=c(50,50,200,200), );
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"))
vis.subject.annot(subjects_dir, 'fsaverage', 'Schaefer2018_400Parcels_7Networks_order', 'rh',  'inflated', views=c('sr'), rgloptions = rgloptions, rglactions = rglactions);

# Plot WSBM consensus partition on brain ----------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"), "snapshot_png"=paste0(output_image_directory,"communities.png"))

#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
wsbm_consensus_lh=as.list(setNames(c(0, consensus_iterative_labels[1:200]), schaefer_atlas_region_names_lh))
wsbm_consensus_rh=as.list(setNames(c(0, consensus_iterative_labels[201:400]), schaefer_atlas_region_names_rh))
#make colormap of Yeo colors
yeo_colors=colorRampPalette(c("#78797A", "#7B287E", "#5CA1C8", "#9F5AA4", "#0A9045", "#B6C988", "#EF9C23", "#E34A53"))

rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage5', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, colormap=yeo_colors, "inflated", views="t4", rgloptions = rgloptions, rglactions = rglactions)

rgloptions=list("windowRect"=c(50,50,200,200));
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"))
vis.region.values.on.subject(subjects_dir, 'fsaverage5', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, colormap=yeo_colors, "inflated", views="sr", rgloptions = rgloptions, rglactions = rglactions)

# Plot the Yeo developmental partition on brain ---------------------------

