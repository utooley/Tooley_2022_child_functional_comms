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

#FIGURE OUT HOW TO ROTATE ON A DIFFERENT AXIS FOR RGL
#This edits the hidden function vis.coloredmeshes.rotating, figure out a way to edit this permanently. 
#Change x=0, y=0, z=1, maybe rpm
trace(fsbrain:::brainview.sr, edit=TRUE)

##########################
####### PLOTTING #########
##########################

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

# Plot Yeo17 atlas on brain  -------------------------------------------
#try to find if blue blob is head/face--plot yeo17
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities_17.png"))
vis.subject.annot(subjects_dir, 'fsaverage', 'Schaefer2018_400Parcels_17Networks_order', 'both',  'inflated', views=c('t4'), rgloptions = rgloptions, rglactions = rglactions);

#try to find if blue blob is head/face?
subjects_dir="~/Applications/freesurfer/subjects"
  
vis.subject.annot(subjects_dir, 'fsaverage', 'PALS_B12_Visuotopic', 'lh',  'inflated', views=c('t4'), rgloptions = rgloptions, rglactions = rglactions);

# Plot the Yeo developmental partition on brain ---------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"

yeo_dev_lh <- yeo_dev_partition$lh.labels
yeo_dev_rh <- yeo_dev_partition$rh.labels

rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',yeo_dev_lh, yeo_dev_rh, "inflated", colormap = yeo_colors,  views="t4", rgloptions = rgloptions, rglactions = rglactions)

#try visualizing the annotation file instead! This enables rotation. Need to copy the annotation into the CBIG/Schaefer/Freesurfer5.3 subjects directory first
#subjects_dir = "/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/lh.yeodev.fsaverage6.annot";
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.subject.annot(subjects_dir, 'fsaverage6', 'yeonets.fsaverage6', 'both',  'inflated', views=c('t4'), rgloptions = rgloptions, rglactions = rglactions);
 
rgloptions=list("windowRect"=c(50,50,200,200));
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"))
vis.subject.annot(subjects_dir, 'fsaverage6', 'yeonets.fsaverage6', 'rh',  'inflated', views=c('sr'), rgloptions = rgloptions, rglactions = rglactions);

####################################
######## YEO DEV ###################
####################################

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

# Partition overlap statistics Adjusted-Rand and NMI ------------------------------
library(igraph)
library(aricode)

## ZRAND or NMI
yeo_dev_full <- c(yeo_dev_lh,yeo_dev_rh)
yeo_full <- c(yeo7_lh, yeo7_rh)
compare(yeo_dev_full, yeo_full, method = "adjusted.rand") #these give the same values as zrand in matlab, but there's no zrand.
compare(yeo_dev_full, yeo_full, method = "nmi")
#use aricode 
NMI_value <- clustComp(yeo_dev_full, yeo_full)$NMI

#permutation testing
distr <- vector()
for (i in 1:100){
  shuffled <- sample(yeo_dev_full, replace = FALSE)
  distr[i] <- clustComp(shuffled, yeo_full)$NID
}
hist(distr)
abline(v =NMI_value)
sum(abs(distr) > NMI_value)/1000  # no values this high

# Examine the switches in community assignment from Yeo7 to Yeo dev ----------------------------
#Where do they differ?
compare_vector <- paste0(as.character(yeo7_lh),as.character(yeo_dev_lh))
table(compare_vector)
remove=c("00","11","22","33","44","55","66","77") #remove the same assignments and look at what is most common
table(compare_vector[!compare_vector %in% remove])
#75 24 67 46 are the most common reassignments from Yeo to Yeo dev on lh and rh, so assign these their own categories
#RH
compare_vector2 <- paste0(as.character(yeo7_rh),as.character(yeo_dev_rh))
table(compare_vector2)
remove=c("00","11","22","33","44","55","66","77") #remove the same assignments and look at what is most common
table(compare_vector2[!compare_vector %in% remove])
brain <- c(compare_vector,compare_vector2) #whole brain
perc <- table(brain[!brain %in% remove])/29588 #perc of total switches, total vertices that switched is 29588
switches <- rep(5,length(brain)) #create a new vector to compare to
switches[brain=="75"] <- 1 #7 to 5
switches[brain=="24"] <- 2 #2 to 4
switches[brain=="67"] <- 3 #6 to 7
switches[brain=="46"] <- 4 #4 to 6
switches[brain %in% remove ] <- 0
switches_lh <- switches[1:40962]
switches_rh <- switches[40963:81924]

#Make a palette for this
switch_colors=colorRampPalette(c("#D3D3D3", "#FF999A", "#C381FF", "#df613a","#d5668f", "#cec2b7")) #gray for 0(same),  ,tan for other switch
barplot(1:6,col=switch_colors(6))

#Plot this vector of switches on fsaverage6 surface
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"switches_in_assignment.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',switches_lh, switches_rh, 
                    "inflated", colormap= switch_colors, views="t4", rgloptions = rgloptions, rglactions = rglactions)

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

ggplot(data=longdata, aes(x=community, y= value))+geom_bar(stat="identity", aes(fill = variable), position = "dodge")+
  labs(y="Surface area (mm^2)")+ ylim(0, 47000) + theme(axis.text.x=element_text(angle=90,hjust=1))

# Viz confidence maps for Yeo7 -------------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo7/"

#Read in the 1000subjects ref silhouettes
silhouette <- readMat(paste0(yeo7_ref_dir,"1000subjects_clusters007_ref.mat"), drop = )
silhouette_lh <-  silhouette$lh.s #this is in fsaverage5 space, so 10k vertices
silhouette_rh <-  silhouette$rh.s

#clip values below 0 to be equal to 0
lower_bound <- ecdf(silhouette_rh)(0)#find the value
clipped_silhouette_rh <- clip.data(silhouette_rh, lower_bound,0.6)
lower_bound <- ecdf(silhouette_lh)(0)#find the value
clipped_silhouette_lh <- clip.data(silhouette_lh, lower_bound,0.6)

#Visuzalize them
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette.png"))
vis.data.on.subject(subjects_dir, 'fsaverage5',clipped_silhouette_lh, clipped_silhouette_rh, "inflated", colormap = colorRampPalette(c("mediumpurple","gray", "white")),  views="t4", rgloptions = rgloptions, rglactions = rglactions)

#Rewrite these to upsample to fsaverage6 so that they're at the same resolution as yeo dev
write.fs.curv(paste0(subjects_dir, "fsaverage5/label/lh.silhouette.fsaverage5.curv"), silhouette_lh[,])
write.fs.curv(paste0(subjects_dir, "fsaverage5/label/rh.silhouette.fsaverage5.curv"), silhouette_rh[,])

#####################
##### Use the bash script to upsample using mri_surf2surf here #####
###############

#Then reload them below
silhouette_lh_fs6 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/lh.silhouette.fsaverage6.curv"))
silhouette_rh_fs6 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/rh.silhouette.fsaverage6.curv"))

#clip them
lower_bound <- ecdf(silhouette_rh_fs6)(0)#find the value
clipped_silhouette_rh <- clip.data(silhouette_rh_fs6, lower_bound,0.6)
lower_bound <- ecdf(silhouette_lh_fs6)(0)#find the value
clipped_silhouette_lh <- clip.data(silhouette_lh_fs6, lower_bound,0.6)

#Visualize them on fsaverage6
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo7/"
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette_fsaverage6.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',clipped_silhouette_lh, clipped_silhouette_rh, "inflated", colormap = colorRampPalette(c("burlywood4","burlywood3","white")),  views="t4", rgloptions = rgloptions, rglactions = rglactions)

# Viz confidence maps for Yeo dev ------------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"

silhouette_lh <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh <-  yeo_dev_partition$rh.s

#clip values below 0 to be equal to 0
lower_bound <- ecdf(silhouette_rh)(0)#find the value
clipped_silhouette_rh <- clip.data(silhouette_rh, lower_bound,0.6)
lower_bound <- ecdf(silhouette_lh)(0)#find the value
clipped_silhouette_lh <- clip.data(silhouette_lh, lower_bound,0.6)

#Visuzalize them
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',clipped_silhouette_lh, clipped_silhouette_rh, "inflated", colormap = colorRampPalette(c("burlywood4","burlywood3","white")),  views="t4", rgloptions = rgloptions, rglactions = rglactions)

# Overlay confidence for adults and kids ----------------------------------
silhouette_lh_dev <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh_dev <-  yeo_dev_partition$rh.s

subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
silhouette_lh_yeo7 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/lh.silhouette.fsaverage6.curv"))
silhouette_rh_yeo7 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/rh.silhouette.fsaverage6.curv"))

silhouette_lh_dev_bin <- silhouette_lh_dev < quantile(silhouette_lh_dev, 0.20) #Thresholded at 20% lowest confidence
silhouette_rh_dev_bin <- silhouette_rh_dev < quantile(silhouette_rh_dev, 0.20)
silhouette_lh_yeo7_bin <- silhouette_lh_yeo7 < quantile(silhouette_lh_yeo7, 0.20)
silhouette_rh_yeo7_bin <- silhouette_rh_yeo7 < quantile(silhouette_rh_yeo7, 0.20)

overlap_lh <- as.numeric((silhouette_lh_dev_bin+silhouette_lh_yeo7_bin)==2)
overlap_rh <- as.numeric((silhouette_rh_dev_bin+silhouette_rh_yeo7_bin)==2)

rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette_overlap_yeo7_yeodev.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',overlap_lh, overlap_rh, "inflated", colormap = colorRampPalette(c("gray", "blue","lightblue", "purple")),  views="t4", rgloptions = rgloptions, rglactions = rglactions)

## RAND or NMI COMPARE
library(igraph)
compare(cbind(silhouette_lh_dev,silhouette_rh_dev), cbind(silhouette_lh_yeo7, silhouette_rh_yeo7), method = "nmi")
clustComp(c(silhouette_lh_dev,silhouette_rh_dev), c(silhouette_lh_yeo7, silhouette_rh_yeo7))

# Sum and difference of confidence maps -----------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"
#clip it at 0 first, then add together both!
silhouette_lh_yeo7[silhouette_lh_yeo7<0] <- 0
silhouette_rh_yeo7[silhouette_rh_yeo7<0] <- 0
silhouette_lh_dev[silhouette_lh_dev<0] <- 0
silhouette_rh_dev[silhouette_rh_dev<0] <- 0

sum_silhouette_devyeo_lh <- silhouette_lh_dev+silhouette_lh_yeo7
sum_silhouette_devyeo_rh <- silhouette_rh_dev+silhouette_rh_yeo7

#Plot sum of confidence maps
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette_summed_yeo7_yeodev.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',sum_silhouette_devyeo_lh, sum_silhouette_devyeo_rh, "inflated", colormap = colorRampPalette(c("burlywood4","burlywood3","white","white","white")),  views="t4", rglactions=rglactions,rgloptions = rgloptions,draw_colorbar = T)

#Calculate differences between the two confidence maps
#Where is there more uncertainty in kids than adults?
difference_sil_lh <- silhouette_lh_dev-silhouette_lh_yeo7 #if confidence in kids is low and adults high, this is negative, if kids high and adults low, this is positive
difference_sil_rh <- silhouette_rh_dev-silhouette_rh_yeo7

rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette_diff_yeo7_yeodev_redhighadults.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',difference_sil_lh, difference_sil_rh, "inflated", colormap = colorRampPalette(c("coral","white","cyan3")),  views="t4",rgloptions = rgloptions, rglactions=rglactions,draw_colorbar = T)

# Plot confidence in community assignment by Yeo7 system -------------------------------------
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("~/Documents/tools/raincloud.R")

silhouette_lh_dev <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh_dev <-  yeo_dev_partition$rh.s
silhouette_dev <- rbind(silhouette_lh_dev,silhouette_rh_dev)
yeo7 <- c(yeo7_lh, yeo7_rh)
data <- data.frame(silhouette_dev, as.character(yeo7))
data <- data[data$as.character.yeo7 != "0",]#remove the medial wall
colnames(data) <- c("silhouette", "yeo7")

## Raincloud plot
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

sumld<- ddply(data, ~yeo7, summarise, mean = mean(silhouette), median = median(silhouette), lower = lb(silhouette), upper = ub(silhouette))

head(sumld)
g <- ggplot(data = data, aes(y = silhouette, x = yeo7, fill = yeo7)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = silhouette, color = yeo7), position = position_jitter(width = .15), size = .08, alpha = 0.15) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme

g

# Plot confidence in community assignment by Yeodev system (Supplement)-------------------------------------
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("~/Documents/tools/raincloud.R")

silhouette_lh_dev <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh_dev <-  yeo_dev_partition$rh.s
silhouette_dev <- rbind(silhouette_lh_dev,silhouette_rh_dev)
yeodev <- c(yeo_dev_lh, yeo_dev_rh)
data <- data.frame(silhouette_dev, as.character(yeodev))
data <- data[data$as.character.yeodev. != "0",]#remove the medial wall
colnames(data) <- c("silhouette", "yeodev")

## Raincloud plot
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)

sumld<- ddply(data, ~yeodev, summarise, mean = mean(silhouette), median = median(silhouette), lower = lb(silhouette), upper = ub(silhouette))

head(sumld)
g <- ggplot(data = data, aes(y = silhouette, x = yeodev, fill = yeodev)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = silhouette, color = yeodev), position = position_jitter(width = .15), size = .08, alpha = 0.15) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme

g
####################################
############# WSBM #################
####################################
# Plot WSBM consensus partition on brain ----------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"), "snapshot_png"=paste0(output_image_directory,"communities.png"))

#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
wsbm_consensus_lh=as.list(setNames(c(0, consensus_iterative_labels[1:200]), schaefer_atlas_region_names_lh))
wsbm_consensus_rh=as.list(setNames(c(0, consensus_iterative_labels[201:400]), schaefer_atlas_region_names_rh))
#make colormap of Yeo colors
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53")) #these are corr
barplot(1:8, col=yeo_colors(8))

rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, colormap=yeo_colors, "inflated", views="t4",rgloptions = rgloptions, rglactions = rglactions)

rgloptions=list("windowRect"=c(50,50,200,200));
rglactions=list("movie"=paste0(output_image_directory,"communities.gif"))
#need only one hemisphere here
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, colormap=yeo_colors,"inflated", views="sr", rgloptions = rgloptions, rglactions = rglactions)

# Compare assignments in Yeo7 to WSBM -------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
yeo_nodes_schaefer400=as.numeric(read.delim('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt', header = F))
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";

#Where do they differ?
comparison_wsbm_yeo <- yeo_nodes_schaefer400==consensus_iterative_labels

#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
comparison_wsbm_yeo_lh=as.list(setNames(c(1, comparison_wsbm_yeo[1:200]), schaefer_atlas_region_names_lh))
comparison_wsbm_yeo_rh=as.list(setNames(c(1, comparison_wsbm_yeo[201:400]), schaefer_atlas_region_names_rh))

#Make a gray and purple color palette
gray_colors=colorRampPalette(c("#9900FF", "#D3D3D3")) #purple for 0, gray for 1

#Visualize them
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"diff_in_assignment.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  comparison_wsbm_yeo_lh, 
                             comparison_wsbm_yeo_rh, colormap=gray_colors, "inflated", views="t4", rgloptions = rgloptions, rglactions = rglactions)

perc_parcels_unmatched=(sum(comparison_wsbm_yeo==0))/(400) #33 % of parcels don't match

## ZRAND or NMI
library(igraph)
compare(yeo_nodes_schaefer400[,1], consensus_iterative_labels, method = "adjusted.rand")
compare(yeo_nodes_schaefer400[,1], consensus_iterative_labels, method = "nmi") #NMI is lower than Yeo7-Yeo dev, but adjusted Rand is higher than Yeo7-Yeo dev

# Partition overlap statistics Adjusted-Rand and NMI ------------------------------
library(igraph)
library(aricode)

## ZRAND or NMI
compare(consensus_iterative_labels, yeo_nodes_schaefer400$V1, method = "adjusted.rand") #these give the same values as zrand in matlab, but there's no zrand.
compare(consensus_iterative_labels, yeo_nodes_schaefer400$V1, method = "nmi")
#use aricode 
clustComp(consensus_iterative_labels[,1], yeo_nodes_schaefer400$V1)
NMI_value <- clustComp(consensus_iterative_labels[,1], yeo_nodes_schaefer400$V1)$NMI

#permutation testing
distr <- vector()
for (i in 1:1000){
  shuffled <- sample(consensus_iterative_labels[,1], replace = FALSE)
  distr[i] <- clustComp(shuffled, yeo_nodes_schaefer400$V1)$NMI
}
hist(distr)
abline(v =NMI_value)
sum(abs(distr) > NMI_value)/100  # no values this high

# Examine the switches in community assignment from Yeo7 to WSBM ----------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
#Where do they differ?
brain <- paste0(as.character(yeo_nodes_schaefer400$V1),as.character(consensus_iterative_labels))
table(brain)
remove=c("00","11","22","33","44","55","66","77") #remove the same assignments and look at what is most common
table(brain[!brain %in% remove]) #132 parcels changed assignment
#42 31 65 57 are the most common reassignments from Yeo to Yeo dev on lh and rh, so assign these their own categories
perc <- table(brain[!brain %in% remove])/132 #perc of total switches, total vertices that switched is 29588
switches <- rep(5,length(brain)) #create a new vector to compare to
switches[brain=="42"] <- 1 #4 to 2
switches[brain=="31"] <- 2 #3 to 1
switches[brain=="65"] <- 3 #6 to 5
switches[brain=="57"] <- 4 #5 to 7
switches[brain %in% remove ] <- 0
#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
switches_lh=as.list(setNames(c(0, switches[1:200]), schaefer_atlas_region_names_lh))
switches_rh =as.list(setNames(c(0, switches[201:400]), schaefer_atlas_region_names_rh))

switches_lh <- switches[1:200]
switches_rh <- switches[201:400]

#Make a palette for this
switch_colors=colorRampPalette(c("#D3D3D3", "#8659FF", "#69FFA7","#FDD5A3", "#FFB3C0", "#cec2b7")) #gray for 0(same),  ,tan for other switch
barplot(1:6,col=switch_colors(6))

#Plot this vector of switches on fsaverage6 surface
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"switches_in_assignment.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  switches_lh, 
                             switches_rh, colormap=switch_colors, "inflated", views="t4", rgloptions = rgloptions, rglactions = rglactions)

# Surface area in WSBM vs. Yeo7 by community ----------------------------
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

#WSBM
#LH
wsbm_surfarea_lh <- read.table("/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/lh.surfarea.stats")
wsbm_surfarea_lh <- wsbm_surfarea_lh %>% select(V4:V5)
#RH
wsbm_surfarea_rh <- read.table("/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/lh.surfarea.stats")
wsbm_surfarea_rh <- wsbm_surfarea_rh %>% select(V4:V5)
#combine them
wsbm_surf_area <- data.frame((wsbm_surfarea_lh$V4+wsbm_surfarea_rh$V4), wsbm_surfarea_lh$V5)
colnames(wsbm_surf_area) <- c("surf_area_mm","community")

#Make them into a long datafram
surf_area <- data.frame(yeo7_surf_area, wsbm_surf_area )
colnames(surf_area) <- c("yeo7", "community", "wsbm","comm" )
longdata <- surf_area %>% select(-comm) %>% melt()

ggplot(data=longdata, aes(x=community, y= value))+geom_bar(stat="identity", aes(fill = variable), position = "dodge")+labs(y="Surface area (mm^2)")+theme(axis.text.x=element_text(angle=90,hjust=1))

#Double-checking the weird large values for Default!

#copy WSBM annotation into local CBIG subjects dir
get.atlas.region.names("wsbm.consensus.fsaverage6", template_subjects_dir = subjects_dir,template_subject='fsaverage6', hemi="lh");
subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.consensus.fsaverage6')

rglactions=list("snapshot_png"=paste0(output_image_directory,"communities_annotation.png"))
vis.subject.annot(subjects_dir, 'fsaverage6', 'wsbm.consensus.fsaverage6', 'both',  'inflated', views=c('t4'), rgloptions = rgloptions, rglactions = rglactions);
#THESE do look right!

# Areas of inconsistent assignment in WSBM --------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
freq=abs(freq-670) #change into the number of non-modal assignments out of 670

#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
freq_lh=as.list(setNames(c(0, freq[1:200]), schaefer_atlas_region_names_lh))
freq_rh=as.list(setNames(c(0, freq[201:400]), schaefer_atlas_region_names_rh))

rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"inconsistent_assignment.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  freq_lh, 
                             freq_rh, colormap=colorRampPalette(c("white","burlywood3","burlywood4")), "inflated", views="t4", rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

# Plot variance in community assignment by Yeo7 system ------------------------
source("~/Documents/tools/raincloud.R")
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)

yeo_nodes=read.table('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt', col.names = "yeo_nodes")
yeo_nodes <- as.character(yeo_nodes$yeo_nodes)
#Read in the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
freq<- t(z$freq)
freq_continuous <- abs(freq-670) #change into the number of non-modal assignments out of 670
data <- data.frame(yeo_nodes,freq_continuous)
colnames(data) <- c('yeo7', 'freq')

## Raincloud plot
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)
sumld<- ddply(data, ~yeo7, summarise, mean = mean(freq), median = median(freq), lower = lb(freq), upper = ub(freq))

head(sumld)
g <- ggplot(data = data, aes(y = freq, x = yeo7, fill = yeo7)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = freq, color = yeo7), position = position_jitter(width = .08), size = .5, alpha = 0.5) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values=c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme

g
# Plot variance in community assignment by WSBM system (supplement)------------------------
#Read in the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
freq<- t(z$freq)
freq_continuous <- abs(freq-670) #change into the number of non-modal assignments out of 670
partitions <- readMat(paste0(wsbm_datadir,"n670_training_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled
data <- data.frame(as.character(consensus_iterative_labels),freq_continuous)
colnames(data) <- c('wsbm', 'freq')
## Raincloud plot
lb <- function(x) mean(x) - sd(x)
ub <- function(x) mean(x) + sd(x)
sumld<- ddply(data, ~wsbm, summarise, mean = mean(freq), median = median(freq), lower = lb(freq), upper = ub(freq))

head(sumld)
g <- ggplot(data = data, aes(y = freq, x = wsbm, fill = wsbm)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = freq, color = wsbm), position = position_jitter(width = .08), size = .5, alpha = 0.5) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme

g

# Overlap of variability in WSBM and Yeo Dev ------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"

#First relabel Schaefer400 annotation with WSBM frequency of inconsistent assignment in matlab
#Then import that annotation, compare both vectors in fsaverage6 space as we did before
#Inconsistent assignment WSBM
lh_freq <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.freq.fsaverage6' )
lh_freq <- as.numeric(as.character(lh_freq$label_names))
rh_freq <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.freq.fsaverage6' )
rh_freq <- as.numeric(as.character(rh_freq$label_names))

#Confidence for Yeo Dev
silhouette_lh <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh <-  yeo_dev_partition$rh.s

silhouette_lh_dev_bin <- silhouette_lh_dev < quantile(silhouette_lh, 0.20) #Threshold yeo dev silhouette at 20% lowest confidence
silhouette_rh_dev_bin <- silhouette_rh_dev < quantile(silhouette_rh, 0.20)

#threshold WSBM community assignment above 0
rh_freq_bin <- rh_freq>0
lh_freq_bin <- lh_freq>0

overlap_lh <- as.numeric((silhouette_lh_dev_bin+lh_freq_bin)==2)
overlap_rh <- as.numeric((silhouette_rh_dev_bin+rh_freq_bin)==2)

rglactions=list("snapshot_png"=paste0(output_image_directory,"overlap_yeodev_confidence_inconsistence_wsbm.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',overlap_lh, overlap_rh, "inflated", colormap = colorRampPalette(c("burlywood4","burlywood3","white","white")),  views="t4", rgloptions = rgloptions, rglactions = rglactions)

# Sum of variability in WSBM and Yeo Dev ------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"

#normalize the frequency of inconsistent assignments to 0-1 and invert so that lower is less consistent assignment
rh_freq_perc <- 1-(rh_freq/334)
lh_freq_perc <- 1-(lh_freq/325)
#sum them together
sum_lh <- silhouette_lh_dev+lh_freq_perc
sum_rh <- silhouette_rh_dev+rh_freq_perc

rglactions=list("snapshot_png"=paste0(output_image_directory,"sum_yeodev_confidence_inconsistence_wsbm.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',sum_lh, sum_rh, "inflated", colormap = colorRampPalette(c("burlywood4","burlywood3","white","white")),  views="t4", rgloptions = rgloptions, rglactions = rglactions)
