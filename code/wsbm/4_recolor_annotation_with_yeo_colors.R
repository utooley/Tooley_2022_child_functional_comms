library(fsbrain) #this may not work if you're editing the script directly on the cluster...
library(freesurferformats)
library(dplyr)
library(R.matlab)
library(stringr)
library(ggplot2)
library(tidyr)
library(igraph)
library(aricode)

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

#RH
wsbm_rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.consensus.fsaverage6')
#load colors
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";
lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','Yeo2011_7Networks_N1000' )

#modify annot parameters to rewrite
colortable=data.frame(lh$colortable_df)
colortable$struct_name=0:7

#might need to add 1 to the labels to start at 1 rather than 0
labels_as_indices_into_colortable=as.numeric(wsbm_rh$label_names)
labels_as_indices_into_colortable=labels_as_indices_into_colortable+1

write.fs.annot(paste0(wsbm_datadir,"rh.wsbm.fsaverage6.annot"), num_vertices=length(wsbm_rh$vertices), 
               colortable=colortable, labels_as_colorcodes=NULL, labels_as_indices_into_colortable=labels_as_indices_into_colortable)

#LH
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
wsbm_lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.consensus.fsaverage6')

#modify annot parameters to rewrite
colortable=data.frame(lh$colortable_df)
colortable$struct_name=0:7

#might need to add 1 to the labels to start at 1 rather than 0
labels_as_indices_into_colortable=as.numeric(wsbm_lh$label_names)
labels_as_indices_into_colortable=labels_as_indices_into_colortable+1

write.fs.annot(paste0(wsbm_datadir,"lh.wsbm.fsaverage6.annot"), num_vertices=length(wsbm_lh$vertices), 
               colortable=colortable, labels_as_colorcodes=NULL, labels_as_indices_into_colortable=labels_as_indices_into_colortable)
