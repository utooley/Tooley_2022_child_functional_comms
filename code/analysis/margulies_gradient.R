library(ggseg3d)
library(ggsegExtra)
library(tidyr)
library(R.matlab)
library(dplyr)
library(fsbrain)
library(stringr)
library(data.table)
library(freesurferformats)
library(readr)
library(tidyr)
library(ggplot2)
library(aricode)
#turn off scientific notation
options(scipen = 9)

# SETUP -------------------------------------------------------------------
#CUBIC cluster filepaths, mounted locally
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
yeo_dev_dir="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_mot_outliers/"
yeo7_ref_dir="/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/"

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

# Margulies gradients -----------------------------------------------------
hcp_gradients_dir="/cbica/projects/spatial_topography/data/imageData/marguliesgradient/"
output_image_directory="/cbica/projects/spatial_topography/data/imageData/marguliesgradient/"
gradients <- readMat(paste0(hcp_gradients_dir,"func_G1_fsa5_sjhparc.mat"), drop = )
mask <- readMat(paste0(hcp_gradients_dir,"mask.mat"), drop = )$mask
gradients[mask==0] <- NA

rglactions=list("snapshot_png"=paste0(output_image_directory,"gradients.png"))
vis.data.on.subject(subjects_dir, 'fsaverage5',gradients[1:10242] , gradients[10243:20484], "inflated", colormap = colorRampPalette(c("darkblue","lightgreen","green", "yellow","red")),  views="t4", draw_colorbar = T,rgloptions = rgloptions, rglactions = rglactions)

#compare to the other 3 partitions
library(igraph)
yeo7 <- readMat(paste0(yeo7_ref_dir,"1000subjects_clusters007_ref.mat"), drop = )
yeo7_fsaverage5 <- c(yeo7$lh.labels, yeo7$rh.labels)

#average gradients values within communities of the 3 different partitions
        
        
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53")) #these are corr
