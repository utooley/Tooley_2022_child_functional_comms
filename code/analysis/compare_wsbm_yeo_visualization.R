library(ggseg3d)
library(ggsegExtra)
library(tidyr)
library(R.matlab)
library(dplyr)
library(fsbrain)
library(stringr)
library(data.table)

# SETUP -------------------------------------------------------------------
#CUBIC cluster filepaths, mounted locally
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
yeo_dev_dir="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_mot_outliers/"

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
#Change the output directory for figures to be on the cluster brains
setwd("/cbica/projects/spatial_topography/output/images/brains/")

#Set RGL defaults, larger window and how to make snapshots!
rgloptions=list("windowRect"=c(50,50,1000,1000));

#FIGURE OUT HOW TO ROTATE ON A DIFFERENT AXIS FOR RGL
#This edits the hidden function vis.coloredmeshes.rotating
trace(fsbrain:::brainview.sr, edit=TRUE)

# Plot Yeo7 atlas on brain ---------------------------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo7/"

#Get the Schaefer atlas region names
schaefer_atlas_region_names_lh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");

rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.subject.annot(subjects_dir, 'fsaverage', 'Schaefer2018_400Parcels_7Networks_order', 'both',  'inflated', views=c('t4'), rgloptions = rgloptions, rglactions = rglactions);

rgloptions=list("windowRect"=c(50,50,200,200));
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

rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, colormap=yeo_colors, "inflated", views="t4", rgloptions = rgloptions, rglactions = rglactions)

rgloptions=list("windowRect"=c(50,50,200,200));
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"))
#need only one hemisphere here
vis.region.values.on.subject(subjects_dir, 'fsaverage5', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, colormap=yeo_colors,"inflated", views="sr", rgloptions = rgloptions, rglactions = rglactions)

# Plot the Yeo developmental partition on brain ---------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"

yeo_dev_lh <- yeo_dev_partition$lh.labels
yeo_dev_rh <- yeo_dev_partition$rh.labels

rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',yeo_dev_lh, yeo_dev_rh, "inflated", colormap = yeo_colors,  views="t4", rgloptions = rgloptions, rglactions = rglactions)

#This is still not working because it tries to display both instead of one! Could just write this to an annotation file, if I wanted.
rgloptions=list("windowRect"=c(50,50,200,200));
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"))
vis.data.on.subject(subjects_dir, 'fsaverage6',yeo_dev_lh, yeo_dev_rh, colormap=yeo_colors,"inflated", views="sr", rgloptions = rgloptions, rglactions = rglactions)

# Compare assignments in Yeo7 to Yeo-dev ----------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"

#other alternative is read in the Freesurfer version of Yeo7 in  fsaverage6 space
#this is more theoretically motivated, because this is how Yeo-dev was generated
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";
lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','Yeo2011_7Networks_N1000' )
yeo7_lh <- as.numeric(str_split(lh$label_names, "_", simplify = TRUE)[,2])
yeo7_lh[is.na(yeo7_lh)] <- 0
rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','Yeo2011_7Networks_N1000' )
yeo7_rh <- as.numeric(str_split(rh$label_names, "_", simplify = TRUE)[,2])
yeo7_rh[is.na(yeo7_rh)] <- 0

#LH
#Where do they differ?
comparison_lh <- yeo_dev_lh==yeo7_lh
compare_yeo7_dev_1_is_same_lh <- as.integer(as.logical(comparison_lh))
#RH
comparison_rh <- yeo_dev_rh==yeo7_rh
compare_yeo7_dev_1_is_same_rh <- as.integer(as.logical(comparison_rh))

#Make a gray and purple color palette
gray_colors=colorRampPalette(c("#9900FF", "#D3D3D3")) #purple for 0, gray for 1

#Plot the vector of differences on fsaverage6 surface
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"diff_in_assignment.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',compare_yeo7_dev_1_is_same_lh, compare_yeo7_dev_1_is_same_rh, 
                    "inflated", colormap= gray_colors, views="t4", rgloptions = rgloptions, rglactions = rglactions)

perc_vertices_unmatched=(sum(compare_yeo7_dev_1_is_same_lh==0)+sum(compare_yeo7_dev_1_is_same_rh==0))/(40962*2)
#actually, ignore medial wall in calculating percentage of vertices, which is 3544 + 3534 = 7078 vertices
perc_vertices_unmatched=(sum(compare_yeo7_dev_1_is_same_lh==0)+sum(compare_yeo7_dev_1_is_same_rh==0))/((40962*2)-7078) #39.5 of vertices unmatched

# Number of vertices in Yeo7 vs Yeo-dev by community ----------------------------
#Sum vertices in each hemisphere, these two have the same total # of vertices, in fsaverage6 space
yeo7 <- table(yeo7_lh) + table(yeo7_rh)
yeo_dev <- table(yeo_dev_lh) + table(yeo_dev_rh)

#Make them into a long datafram
num_vertices <- data.frame(yeo7, yeo_dev)
colnames(num_vertices) <- c("comm",    "yeo7"    ,   "community", "yeo_dev" )
num_vertices$community <- factor(num_vertices$community, ordered = TRUE, levels = c("Medial wall", "VIS", "SM", "DAN", "VAN", "LIM", "FPN", "DMN"))
longdata <- num_vertices %>% select(-comm) %>% melt()

ggplot(data=longdata, aes(x=community, y= value))+geom_bar(stat="identity", aes(fill = variable), position = "dodge")+labs(y="Surface area (# of vertices)")

# Surface area in Yeo7 vs Yeo-dev by community ----------------------------
#YEO7
#LH
yeo7_surfarea_lh <- read.table("/cbica/projects/spatial_topography/tools/parcellations/yeo7_from_freesurfer/fsaverage6/lh.7nets.surfarea.stats")
yeo7_surfarea_lh <- yeo7_surfarea_lh %>% select(V4:V5)
#RH
yeo7_surfarea_rh <- read.table("/cbica/projects/spatial_topography/tools/parcellations/yeo7_from_freesurfer/fsaverage6/rh.7nets.surfarea.stats")
yeo7_surfarea_rh <- yeo7_surfarea_rh %>% select(V4:V5)
#combine them
yeo7_surf_area <- data.frame((yeo7_surfarea_lh$V4+yeo7_surfarea_rh$V4), yeo7_surfarea_lh$V5)
colnames(yeo7_surf_area) <- c("surf_area_mm","community")

#YEO-DEV
#LH
yeodev_surfarea_lh <- read.table("/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/freesurfer/lh.surfarea.stats")
yeodev_surfarea_lh <- yeodev_surfarea_lh %>% select(V4:V5)
#RH
yeodev_surfarea_rh <- read.table("/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/freesurfer/rh.surfarea.stats")
yeodev_surfarea_rh <- yeodev_surfarea_rh %>% select(V4:V5)
#combine them
yeodev_surf_area <- data.frame((yeodev_surfarea_lh$V4+yeodev_surfarea_rh$V4), yeodev_surfarea_lh$V5)
colnames(yeodev_surf_area) <- c("surf_area_mm","community")

#Make them into a long datafram
surf_area <- data.frame(yeo7_surf_area, yeodev_surf_area )
colnames(surf_area) <- c("yeo7", "community", "yeo_dev","comm" )
longdata <- surf_area %>% select(-comm) %>% melt()

ggplot(data=longdata, aes(x=community, y= value))+geom_bar(stat="identity", aes(fill = variable), position = "dodge")+labs(y="Surface area (mm^2)")+theme(axis.text.x=element_text(angle=90,hjust=1))

# Compare assignments in Yeo7 to WSBM -------------------------------------
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";

#read in the Schaefer400 Yeo7 communities in fsaverage6 surface space
x <- subject.annot(subjects_dir, 'fsaverage6', 'rh','Schaefer2018_400Parcels_7Networks_order' )

