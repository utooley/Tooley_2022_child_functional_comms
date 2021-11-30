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
#yeo7_ref_dir="/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/"

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

##########################
####### PLOTTING #########
##########################

# Plot Yeo7 atlas on brain ---------------------------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo7/"

rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.subject.annot(yeo7_ref_dir, 'fsaverage','Yeo2011_7Networks_N1000','both', "inflated", views="t4",
                  rgloptions = rgloptions,rglactions = rglactions)

# Plot the Yeo developmental partition on brain ---------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"

yeo_dev_lh <- yeo_dev_partition$lh.labels
yeo_dev_rh <- yeo_dev_partition$rh.labels

rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
makecmap_options=list('colFn'=yeo_colors)
vis.data.on.subject(subjects_dir, 'fsaverage6',yeo_dev_lh, yeo_dev_rh, "inflated", makecmap_options = makecmap_options, views="t4", rgloptions = rgloptions)

# Plot WSBM consensus partition on brain ----------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"

#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
wsbm_consensus_lh=as.list(setNames(c(0, consensus_iterative_labels[1:200]), schaefer_atlas_region_names_lh))
wsbm_consensus_rh=as.list(setNames(c(0, consensus_iterative_labels[201:400]), schaefer_atlas_region_names_rh))

#make colormap of Yeo colors
yeo_colors=colorRampPalette(c("#000000", "#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#dcf8a4", "#EF9C23", "#E34A53")) #these are corr
makecmap_options=list('colFn'=yeo_colors)
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"communities.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  wsbm_consensus_lh, 
                             wsbm_consensus_rh, makecmap_options = makecmap_options, "inflated", views="t4",rgloptions = rgloptions, rglactions = rglactions)

#make schaefer parcellation figure
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"no_flat_pial.png"))
vis.subject.annot(subjects_dir, 'fsaverage', 'Schaefer2018_400Parcels_7Networks_order', hemi = "lh", surface = "pial",
                  rgloptions = rgloptions, rglactions = rglactions, outline=list('background'= "gray",'outline_color'='gray', 'expand_inwards' = 0L))

# Search over k: Yeo dev -------------------------------------
k_dir="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/search_over_k/"
data <- read.csv(paste0(k_dir,"k2_to_9_tries100_rand100.znorm1.dim1.csv")) %>% select(-contains(c("_1", "_9")))
data2 <-  read.csv(paste0(k_dir,"k9_to_17_tries100_rand100.znorm1.dim1.csv")) %>% select(-contains(paste0("_.",1:8)))
data <- cbind(data,data2)
longdata <- melt(data, varnames = c("x","k"))
longdata$k <-as.numeric(gsub("[^0-9-]", "", longdata$variable))
longdata$type <- str_extract(longdata$variable, "[a-z]+_[a-z]+" )
longdata <- arrange(longdata, as.numeric(k))
longdata$value[longdata$value==0] <- NA
#plot stability
dat <- longdata %>% filter(type=="stability_true") %>% filter(k!=1)
p<-ggplot(dat, aes(x=k, y=value)) + geom_jitter(position=position_jitter(0.2), cex=1, alpha=0.05)+stat_summary(aes(group=type),
                                                                                                               fun.y = mean, geom="line") + theme_classic()+scale_x_continuous(breaks=seq(2,17))+ylab(label = "Instability")
p

# Search over k: WSBM -------------------------------------
k_dir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/search_over_k/"
logevidence=matrix(NA,670,18)
for (i in 1:17){
  k=i+1
  try(mine <- read.csv(paste0(k_dir, "wsbm_search_over_k",k,"_n670_site16_30trials.csv")))
  logevidence[,k] <- with(mine, if("subjlist_.4" %in% colnames(mine)) subjlist_.4 else Var4)
}
#remove the outlier subject sub-NDARINVK0U3PRRR
logevidence <- logevidence[mine$Var1!="sub-NDARINVK0U3PRRR",]
logevidence[logevidence==0] <- NA
longdata_training <- melt(logevidence, varnames = c("x","k")) %>% filter(k != 1 & k != 18)
p<-ggplot(longdata_training, aes(x=k, y=value)) + geom_jitter(position=position_jitter(0.2), cex=1, alpha=0.5)+stat_summary(aes(group=1),
                                                                                                                            fun.y = mean, geom="line", col=rgb(101, 58, 150, maxColorValue = 255)) + theme_classic()
p

# Replication dataset, add plot overlaid
k_dir="/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample/brains/"
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
#add different data file for k=7, "n544_test_sample_",k,"consensus_log_evidence_replicate_7try.csv"
partition_logevidence=matrix(NA,544,18)
num_comms=numeric()
for (i in 2:17){
  k=i
  print(k)
  try(mine <- read.csv(paste0(k_dir, "n544_test_sample_",k,"consensus_log_evidence_replicate.csv")))
  partition_logevidence[,i] <- mine$LP
  partition <- readMat(paste0(wsbm_datadir,"n670_k",k,"training_sample_consensus_partition.mat"), drop = )
  print(range(unique(partition$consensus.iter.mode)))
  num_comms[i] <- max(unique(partition$consensus.iter.mode))
}
partition_logevidence[partition_logevidence==0] <- NA
longdata_testing <- melt(partition_logevidence, varnames = c("x","k")) %>% filter(k != 1 & k != 18) %>% arrange(k)

p<-ggplot(data=longdata_training, aes(x=k, y=value))+geom_point(position=position_jitter(0.2), cex=1,alpha=0.05,col="black") + 
  stat_summary(data=longdata_training, aes(group=1),fun.y = mean, geom="line", col="black") + 
  theme_classic() +
  labs(x="k (number of communities)", y="Log-likelihood (main sample)")
p
ggplot_build(p)$layout$panel_scales_y[[1]]$range$range

#put it on top of the training sample with 2 y-axes!
longdata_testing$value_trans <- longdata_testing$value*(4.5)+(290000)
l <- p+geom_point(data=longdata_testing, aes(x=k, y=value_trans),position=position_jitter(0.2), cex=1, alpha=0.05, col="#5e195e") +
  scale_y_continuous(limits=c(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range),
                     sec.axis = sec_axis(~ . -(290000/4.5), name = "Log-likelihood (replication sample)"))+
  stat_summary(data=longdata_testing, aes(y=value_trans,group=1), fun.y = mean, geom="line", col="#5e195e") +
  scale_x_continuous(breaks=seq(2,17))
l

b<-ggplot(data=longdata_testing, aes(x=k, y=value))+geom_point(position=position_jitter(0.2), cex=1,alpha=0.1) + 
  stat_summary(data=longdata_testing, aes(group=1),fun.y = mean, geom="line", col=rgb(101, 58, 150, maxColorValue = 255)) + 
  theme_classic() 
b
ggplot_build(b)$layout$panel_scales_y[[1]]$range$range

#plot the number of communities detected at each k
num_comms <- na.omit(num_comms)
longdata <- data.frame(c(2:17),num_comms) %>% filter(num_comms != 1);colnames(longdata) <- c("k","num_comms")
p<-ggplot(data=longdata, aes(x=k, y=num_comms))+geom_point(alpha=0.5)+stat_summary(aes(group=1), fun.y = mean, geom="line", col=rgb(101, 58, 150, maxColorValue = 255)) + theme_classic()+scale_x_continuous(breaks=seq(2,17))
p

####################################
######## YEO DEV ###################
####################################

# Compare assignments in Yeo7 to Yeo-dev ----------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"

#read in the Freesurfer version of Yeo7 in  fsaverage6 space
#this is more theoretically motivated, because this is how we generate the developmental clustering partition
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

perc_vertices_unmatched=(sum(compare_yeo7_dev_1_is_same_lh==0)+sum(compare_yeo7_dev_1_is_same_rh==0))/(40962*2)
#actually, ignore medial wall in calculating percentage of vertices, which is 3544 + 3534 = 7078 vertices
perc_vertices_unmatched=(sum(compare_yeo7_dev_1_is_same_lh==0)+sum(compare_yeo7_dev_1_is_same_rh==0))/((40962*2)-7078) #39.5 of vertices unmatched

# Partition overlap statistics Adjusted-Rand and NMI ------------------------------
## ZRAND or NMI
yeo_dev_full <- c(yeo_dev_lh,yeo_dev_rh)
yeo_full <- c(yeo7_lh, yeo7_rh)
compare(yeo_dev_full, yeo_full, method = "adjusted.rand") #these give the same values as zrand in matlab, but there's no zrand.
compare(yeo_dev_full, yeo_full, method = "nmi")
#use aricode 
NMI_value <- clustComp(yeo_dev_full, yeo_full)$NMI

#permutation testing
distr <- vector()
for (i in 1:10000){
  shuffled <- sample(yeo_dev_full, replace = FALSE)
  distr[i] <- clustComp(shuffled, yeo_full)$NMI
}
hist(distr)
abline(v =NMI_value)
sum(abs(distr) > NMI_value)/10000  # no values this high

# Examine the switches in community assignment from Yeo7 to Yeo dev ----------------------------
compare_vector <- paste0(as.character(yeo7_lh),as.character(yeo_dev_lh))
table(compare_vector)
remove=c("00","11","22","33","44","55","66","77") #remove the same assignments and look at what is most common
table(compare_vector[!compare_vector %in% remove])
#RH
compare_vector2 <- paste0(as.character(yeo7_rh),as.character(yeo_dev_rh))
table(compare_vector2)
remove=c("00","11","22","33","44","55","66","77") #remove the same assignments and look at what is most common
table(compare_vector2[!compare_vector %in% remove])
brain <- c(compare_vector,compare_vector2) 
perc <- table(brain[!brain %in% remove])/29588 #perc of total switches, total vertices that switched is 29588
29588/74846 #number of switches over number of total vertices ignoring medial wall
switches <- rep(5,length(brain)) #create a new vector to compare to
#75 24 67 46 are the most common reassignments from Yeo to Yeo dev on lh and rh, so assign these their own categories
switches[brain=="75"] <- 1 #7 to 5
switches[brain=="24"] <- 2 #2 to 4
switches[brain=="67"] <- 3 #6 to 7
switches[brain=="46"] <- 4 #4 to 6
switches[brain %in% remove ] <- 0
switches_lh <- switches[1:40962]
switches_rh <- switches[40963:81924]

switches_yeodev <- ifelse(switches==0,0,1)

#Make a palette for this
switch_colors=colorRampPalette(c("#f2f2f2", "#FF999A", "#C381FF", "#df613a","#d5668f", "#cec2b7")) #gray for 0(same), tan for other switch
barplot(1:6,col=switch_colors(6))

#Plot this vector of switches on fsaverage6 surface
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"switches_in_assignment.png"), 'shift_hemis_apart'=list('min_dist'=20))
makecmap_options=list('colFn'=switch_colors)
vis.data.on.subject(subjects_dir, 'fsaverage6',switches_lh, switches_rh, 
                     "inflated", views="t9", makecmap_options=makecmap_options, rgloptions = rgloptions, rglactions = rglactions)

# Surface area in Yeo7 vs Yeo-dev by community ----------------------------
library(reshape2)
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

#Make them into a long dataframe
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

##### Use the bash script 2_upsample_to_fsaverage6.sh to upsample these here #####

#Then reload them below
silhouette_lh_fs6 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/lh.silhouette.fsaverage6.curv"))
silhouette_rh_fs6 <- read.fs.curv(paste0(subjects_dir, "fsaverage6/label/rh.silhouette.fsaverage6.curv"))
clipped_silhouette_lh <- silhouette_lh_fs6
clipped_silhouette_rh <- silhouette_rh_fs6
clipped_silhouette_lh[clipped_silhouette_lh==1] <- NA
clipped_silhouette_rh[clipped_silhouette_rh==1] <- NA

#No longer needed, now use the 'range' argument to cmap_options below
#clip them for under 0 is 0, above 0.6 to be equal to max
# lower_bound <- ecdf(silhouette_rh_fs6)(0)#find the value
# clipped_silhouette_rh <- clip.data(silhouette_rh_fs6, lower_bound,ecdf(silhouette_rh_fs6)(0.6))
# lower_bound <- ecdf(silhouette_lh_fs6)(0)#find the value
# clipped_silhouette_lh <- clip.data(silhouette_lh_fs6, lower_bound,ecdf(silhouette_lh_fs6)(0.6))

#Visualize them on fsaverage6
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo7/"
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette_fsaverage6.png"))
colFn_diverging = colorRampPalette(c("burlywood4","burlywood3","white"));
makecmap_options=list('range'=c(0,0.5),'colFn'=colFn_diverging)
vis.data.on.subject(subjects_dir, 'fsaverage6',clipped_silhouette_lh, clipped_silhouette_rh, "inflated",  views="t4", rgloptions = rgloptions, rglactions = rglactions,
                    makecmap_options = makecmap_options,draw_colorbar = TRUE)

# Viz confidence maps for Yeo dev ------------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo_dev/"

silhouette_lh <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh <-  yeo_dev_partition$rh.s
clipped_silhouette_lh <- silhouette_lh;clipped_silhouette_rh <- silhouette_rh
clipped_silhouette_lh[clipped_silhouette_lh==1] <- NA;clipped_silhouette_rh[clipped_silhouette_rh==1] <- NA

#No longer needed, now use the 'range' argument to cmap_options below
#clip values below 0 to be equal to 0, above 0.6 to be equal to max
# lower_bound <- ecdf(silhouette_rh)(0)#find the value
# clipped_silhouette_rh <- clip.data(silhouette_rh, lower_bound,ecdf(silhouette_rh)(0.8))
# lower_bound <- ecdf(silhouette_lh)(0)#find the value
# clipped_silhouette_lh <- clip.data(silhouette_lh, lower_bound,ecdf(silhouette_lh)(0.8))

#Visuzalize them
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette.png"))
colFn_diverging = colorRampPalette(c("burlywood4","burlywood3","white"));
makecmap_options=list('range'=c(0,0.5),'colFn'=colFn_diverging)
vis.data.on.subject(subjects_dir, 'fsaverage6',clipped_silhouette_lh, clipped_silhouette_rh, "inflated", makecmap_options = makecmap_options,  views="t4", rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

# Supplemental fig 2--negative values of confidence -----------------------
#visualize only the 0 confidence values
silhouette_lh_fs6_min <- ifelse(silhouette_lh>0, 0,silhouette_lh)
silhouette_rh_fs6_min <- ifelse(silhouette_rh>0, 0,silhouette_rh)
vis.data.on.subject(subjects_dir, 'fsaverage6',silhouette_lh_fs6_min, silhouette_rh_fs6_min, "inflated", colormap = colorRampPalette(c("blue","blue","white")),  views="t4", draw_colorbar = T, rgloptions = rgloptions)
#overlay on the communities
yeodev_lh_mask=yeo_dev_lh, yeodev_rh_mask=yeo_dev_rh
yeodev_lh_mask[silhouette_lh_fs6_min<0] <- NA
yeodev_rh_mask[silhouette_rh_fs6_min<0] <- NA
vis.data.on.subject(subjects_dir, 'fsaverage6',yeodev_lh_mask, yeodev_rh_mask, "inflated", colormap = yeo_colors,  views="t4", draw_colorbar = T, rgloptions = rgloptions)

#Yeo7 values of confidence
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/yeo7/"
#Read in the 1000subjects ref silhouettes
silhouette <- readMat(paste0(yeo7_ref_dir,"1000subjects_clusters007_ref.mat"), drop = )
silhouette_lh <-  silhouette$lh.s #this is in fsaverage5 space, so 10k vertices
silhouette_rh <-  silhouette$rh.s
#mask and visualize
silhouette_lh_fs5_min <- ifelse(silhouette_lh>0, 0,silhouette_lh)
silhouette_rh_fs5_min <- ifelse(silhouette_rh>0, 0,silhouette_rh)
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette_neg_values_fsaverage5.png"))
makecmap_options=list('colFn'=colorRampPalette(c("blue","blue","white")))
vis.data.on.subject(subjects_dir, 'fsaverage5',silhouette_lh_fs5_min, silhouette_rh_fs5_min, "inflated", makecmap_options = makecmap_options,  views="t4", draw_colorbar = T,rglactions = rglactions, rgloptions = rgloptions)

#overlay on the communities
#read in the Freesurfer version of Yeo7 in  fsaverage5 (!) space
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";
lh <- subject.annot(subjects_dir, 'fsaverage5', 'lh','Yeo2011_7Networks_N1000' )
yeo7_lh <- as.numeric(str_split(lh$label_names, "_", simplify = TRUE)[,2])
yeo7_lh[is.na(yeo7_lh)] <- 0
rh <- subject.annot(subjects_dir, 'fsaverage5', 'rh','Yeo2011_7Networks_N1000' )
yeo7_rh <- as.numeric(str_split(rh$label_names, "_", simplify = TRUE)[,2])
yeo7_rh[is.na(yeo7_rh)] <- 0
yeo7_lh_mask=yeo7_lh; yeo7_rh_mask=yeo7_rh
yeo7_lh_mask[silhouette_lh_fs5_min<0] <- NA
yeo7_rh_mask[silhouette_rh_fs5_min<0] <- NA
makecmap_options=list('colFn'=yeo_colors)
rglactions=list("snapshot_png"=paste0(output_image_directory,"silhouette_neg_and_comms_fsaverage5.png"))
vis.data.on.subject(subjects_dir, 'fsaverage5',yeo7_lh_mask, yeo7_rh_mask, "inflated", makecmap_options = makecmap_options,  views="t4", draw_colorbar = T,rglactions = rglactions, rgloptions = rgloptions)

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

# Plot confidence in community assignment by Yeo7 system-supplementary fig 1a -------------------------------------
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("~/Documents/tools/raincloud.R")

silhouette_lh_dev <-  yeo_dev_partition$lh.s #this is in fsaverage6 space, so 40k vertices
silhouette_rh_dev <-  yeo_dev_partition$rh.s
silhouette_dev <- c(silhouette_lh_dev,silhouette_rh_dev) 
yeo7 <- c(yeo7_lh, yeo7_rh)
data <- data.frame(silhouette_dev, as.character(yeo7))
data <- data[data$as.character.yeo7 != "0",]#remove the medial wall
colnames(data) <- c("silhouette", "yeo7")
data <- data %>% filter(silhouette!=1) #remove the 1's for silhouette around the medial wall

## Raincloud plot
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

#Stats
kruskal.test(silhouette~yeo7, data = data)
pairwise.wilcox.test(data$silhouette, data$yeo7,
                     p.adjust.method = "bonferroni")


# Plot confidence in community assignment by Yeodev system-Fig 3d-------------------------------------
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
data <- data %>% filter(silhouette!=1) #remove the 1's for silhouette around the medial wall

## Raincloud plot
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

#Stats
kruskal.test(silhouette~yeodev, data = data)
pairwise.wilcox.test(data$silhouette, data$yeodev,
                     p.adjust.method = "bonferroni")

# Plot adult confidence in community assignment by Yeo7 system- Fig 3c ------------------------------------
silhouette_yeo <- c(silhouette_lh_fs6,silhouette_rh_fs6)
yeo7 <- c(yeo7_lh, yeo7_rh)
data <- data.frame(silhouette_yeo, as.character(yeo7))
data <- data[data$as.character.yeo7. != "0",]#remove the medial wall
colnames(data) <- c("silhouette", "yeo7")
data <- data %>% filter(silhouette!=1) #remove the 1's for silhouette around the medial wall

#get y-axis limits from the plot above
ylim <- ggplot_build(g)$layout$panel_scales_y[[1]]$range$range
## Raincloud plot
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
  coord_cartesian(ylim=c(-0.2932246,0.8143016))+
  raincloud_theme

g

#Stats
kruskal.test(silhouette~yeo7, data = data)
pairwise.wilcox.test(data$silhouette, data$yeodev,
                     p.adjust.method = "bonferroni")

####################################
############# WSBM #################
####################################
# Compare assignments in Yeo7 to WSBM -------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
yeo_nodes_schaefer400=read.delim('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt', header = F)
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";

#Where do they differ?
comparison_wsbm_yeo <- yeo_nodes_schaefer400==consensus_iterative_labels

#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
comparison_wsbm_yeo_lh=as.list(setNames(c(1, comparison_wsbm_yeo[1:200]), schaefer_atlas_region_names_lh))
comparison_wsbm_yeo_rh=as.list(setNames(c(1, comparison_wsbm_yeo[201:400]), schaefer_atlas_region_names_rh))

perc_parcels_unmatched=(sum(comparison_wsbm_yeo==0))/(400) #33 % of parcels don't match

## ZRAND or NMI
library(igraph)
compare(yeo_nodes_schaefer400[,1], consensus_iterative_labels, method = "adjusted.rand")
compare(yeo_nodes_schaefer400[,1], consensus_iterative_labels, method = "nmi") #NMI is lower than Yeo7-Yeo dev, but adjusted Rand is higher than Yeo7-Yeo dev

# Partition overlap statistics Adjusted-Rand and NMI ------------------------------
#first, copy WSBM annotation created with wsbm/viz_write_annot.m into local CBIG subjects dir
wsbm_rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.consensus.fsaverage6')
wsbm_rh<- as.numeric(as.character(wsbm_rh$label_names))
wsbm_lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.consensus.fsaverage6')
wsbm_lh<- as.numeric(as.character(wsbm_lh$label_names))
yeo7 <- c(yeo7_lh,yeo7_rh)

## ZRAND or NMI ON PARCELS
compare(consensus_iterative_labels, yeo_nodes_schaefer400$V1, method = "adjusted.rand") #these give the same values as zrand in matlab, but there's no zrand.
compare(consensus_iterative_labels, yeo_nodes_schaefer400$V1, method = "nmi")
#use aricode 
clustComp(consensus_iterative_labels[,1], yeo_nodes_schaefer400$V1)
NMI_value <- clustComp(consensus_iterative_labels[,1], yeo_nodes_schaefer400$V1)$NMI

## ZRAND or NMI ON VERTICES
compare(c(wsbm_lh,wsbm_rh), yeo7, method = "adjusted.rand") #these give the same values as zrand in matlab, but there's no zrand.
compare(c(wsbm_lh,wsbm_rh), yeo7, method = "nmi")
#use aricode 
clustComp(c(wsbm_lh,wsbm_rh), yeo7)
NMI_value <- clustComp(c(wsbm_lh,wsbm_rh), yeo7)$NMI

#permutation testing
distr <- vector()
for (i in 1:1000){
  shuffled <- sample(c(wsbm_lh,wsbm_rh), replace = FALSE)
  distr[i] <- clustComp(shuffled, yeo7)$NMI
}
hist(distr)
abline(v =NMI_value)
sum(abs(distr) > NMI_value)/100  # no values this high

# Look at switches in community assignment from Yeo7 to WSBM in VERTEX-Fig 4e --------
#copy WSBM annotation into local CBIG subjects dir
wsbm_rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.consensus.fsaverage6')
wsbm_rh<- as.numeric(as.character(wsbm_rh$label_names))
wsbm_lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.consensus.fsaverage6')
wsbm_lh<- as.numeric(as.character(wsbm_lh$label_names))
yeo7 <- c(yeo7_lh,yeo7_rh)
#examine switches vertex-wise
brain <- paste0(yeo7,as.character(c(wsbm_lh,wsbm_rh)))
table(brain)
remove=c("00","11","22","33","44","55","66","77") #remove the same assignments and look at what is most common
table(brain[!brain %in% remove]) #132 parcels changed assignment
#42 31 65 57 are the most common reassignments from Yeo to Yeo dev on lh and rh, so assign these their own categories
perc <- table(brain[!brain %in% remove])/27866 #perc of total switches, total vertices that switched is 27866
sort(perc)
switches <- rep(6,length(brain)) #create a new vector to compare to
switches[brain=="42"] <- 1 #4 to 2
switches[brain=="57"] <- 2 #5 to 7
switches[brain=="31"] <- 3 #3 to 1
switches[brain=="56"] <- 4 #5 to 6
switches[brain=="67"] <- 5 #6 to 7
switches[brain %in% remove ] <- 0
#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
switches_lh=as.list(setNames(c(0, switches[1:40962]), schaefer_atlas_region_names_lh))
switches_rh =as.list(setNames(c(0, switches[40963:81924]), schaefer_atlas_region_names_rh))

switches_lh <- switches[1:40962]
switches_rh <- switches[40963:81924]

#Make a palette for this
switch_colors=colorRampPalette(c("#f2f2f2", "#8659FF","#FFB3C0", "#69FFA7","#FDD5A3", "#e1552f","#cec2b7")) #gray for 0(same),  ,tan for other switch
barplot(1:7,col=switch_colors(7))

#Plot this vector of switches on fsaverage6 surface
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"switches_in_assignment.png"), 'shift_hemis_apart'=list('min_dist'=20))
makecmap_options=list('colFn'=switch_colors)
vis.data.on.subject(subjects_dir, 'fsaverage6',switches_lh, switches_rh, 
                    "inflated", views="t9", makecmap_options=makecmap_options, rgloptions = rgloptions, rglactions = rglactions)

rglactions=list("snapshot_png"=paste0(output_image_directory,"communities_annotation.png"))
vis.subject.annot(subjects_dir, 'fsaverage6', 'wsbm.consensus.fsaverage6', 'both',  'inflated', views=c('t4'), rgloptions = rgloptions);
#THESE do look right! The sum of vertices here is the same as in surfarea.stats.

# Gradient of switches-Supp Fig S7 ----------------------------------------------------
switches_wsbm <- ifelse(switches==0,0,1)
sum <- switches_wsbm+switches_yeodev #switches from adults

output_image_directory="/cbica/projects/spatial_topography/output/images/brains/"
switch_colors=colorRampPalette(c("white","#BFD3E6","#8C96C6")) #gray for 0(same), tan for other switch
barplot(1:3,col = switch_colors(3))
makecmap_options=list('colFn'=switch_colors)
rglactions=list("snapshot_png"=paste0(output_image_directory,"gradient_of_switches_from_adult.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',sum[1:40962], sum[40963:81924], 
                    "inflated", views="t4", makecmap_options=makecmap_options, rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

#What switches between Yeo-dev and WSBM?
brain <- paste0(yeo_dev,as.character(c(wsbm_lh,wsbm_rh)))
table(brain)
remove=c("00","11","22","33","44","55","66","77") #remove the same assignments and look at what is most common
sum(table(brain[!brain %in% remove]))
#42 31 65 57 are the most common reassignments from Yeo to Yeo dev on lh and rh, so assign these their own categories
perc <- table(brain[!brain %in% remove])/43013 #perc of total switches, total vertices that switched is 27866
sort(perc)
switches <- rep(6,length(brain)) #create a new vector to compare to
switches[brain %in% remove ] <- 0
switches_yeodev_wsbm <- ifelse(switches==6,1,0)

sum_all <- switches_wsbm+switches_yeodev+switches_yeodev_wsbm #switches between all three
switch_colors=colorRampPalette(c("white","#BFD3E6","#8C96C6","#88419D")) #gray for 0(same),  ,tan for other switch
rglactions=list("snapshot_png"=paste0(output_image_directory,"gradient_of_switches_all_3.png"))
makecmap_options=list('colFn'=switch_colors)
vis.data.on.subject(subjects_dir, 'fsaverage6',sum_all[1:40962], sum_all[40963:81924], 
                    "inflated", views="t4", makecmap_options=makecmap_options, rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

# Surface area in WSBM vs. Yeo7 by community-Fig 4d ----------------------------
#YEO7
#LH
yeo7_surfarea_lh <- read.table("/cbica/projects/spatial_topography/tools/parcellations/yeo7_from_freesurfer/fsaverage6/lh.7nets.surfarea.stats")
#RH
yeo7_surfarea_rh <- read.table("/cbica/projects/spatial_topography/tools/parcellations/yeo7_from_freesurfer/fsaverage6/rh.7nets.surfarea.stats")
#combine them
yeo7_surf_area <- data.frame((yeo7_surfarea_lh$V4+yeo7_surfarea_rh$V4), yeo7_surfarea_lh$V5)
colnames(yeo7_surf_area) <- c("surf_area_mm","community")

#WSBM
#LH
wsbm_surfarea_lh <- read.table("/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/lh.surfarea.stats")
#RH
wsbm_surfarea_rh <- read.table("/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/rh.surfarea.stats")
#combine them
wsbm_surf_area <- data.frame((wsbm_surfarea_lh$V4+wsbm_surfarea_rh$V4), wsbm_surfarea_lh$V5)
colnames(wsbm_surf_area) <- c("surf_area_mm","community")

#Make them into a long dataframe
surf_area <- data.frame(yeo7_surf_area, wsbm_surf_area)
colnames(surf_area) <- c("yeo7", "community", "wsbm","comm" )
longdata <- surf_area %>% select(-comm) %>% melt()

ggplot(data=longdata, aes(x=community, y= value))+geom_bar(stat="identity", aes(fill = variable), position = "dodge")+labs(y="Surface area (mm^2)")+theme(axis.text.x=element_text(angle=90,hjust=1))

# Areas of inconsistent assignment in WSBM-Fig 5a --------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
freq=abs(freq-670) #change into the number of non-modal assignments out of 670

#Create a list from the vector of Schaefer structure names and the labels from WSBM consensus partition
freq_lh=as.list(setNames(c(0, freq[1:200]), schaefer_atlas_region_names_lh))
freq_rh=as.list(setNames(c(0, freq[201:400]), schaefer_atlas_region_names_rh))

colormap=list('colFn'=colorRampPalette(c("white","burlywood3","burlywood4")))
rgloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"inconsistent_assignment_thresh05.png"), 'shift_hemis_apart'=list('min_dist'=20))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  freq_lh, 
                             freq_rh, "inflated", makecmap_options = colormap, views="t9", rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

# Areas of low confidence in the WSBM-Supp Fig S6a -------------------------------------
#load silhouette data
z <- readMat(paste0(wsbm_datadir,"silhouettes_wsbm.mat"), drop = )
silhouette_wsbm <- z$silhouette.wsbm.1.diag

silhouette_wsbm_lh=as.list(setNames(c(1, silhouette_wsbm[1:200]), schaefer_atlas_region_names_lh))
silhouette_wsbm_rh=as.list(setNames(c(1, silhouette_wsbm[201:400]), schaefer_atlas_region_names_rh))
silhouette_wsbm_lh

gloptions=list("windowRect"=c(50,50,1000,1000));
rglactions=list("snapshot_png"=paste0(output_image_directory,"confidence_wsbm.png"))
vis.region.values.on.subject(subjects_dir, 'fsaverage6', 'Schaefer2018_400Parcels_7Networks_order',  silhouette_wsbm_lh,
                             silhouette_wsbm_rh, colormap=colorRampPalette(c("burlywood4","burlywood3","white","white", "white")), "inflated", views="t4", rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

# Plot variance in community assignment by Yeo7 system-Supp Fig S8b ------------------------
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

#Do it by VERTEX instead
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
get.atlas.region.names("wsbm.consensus.fsaverage6", template_subjects_dir = subjects_dir,template_subject='fsaverage6', hemi="rh");
wsbm_rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.freq.fsaverage6')
wsbm_rh<- as.numeric(as.character(wsbm_rh$label_names))
wsbm_lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.freq.fsaverage6')
wsbm_lh<- as.numeric(as.character(wsbm_lh$label_names))

yeo7 <- c(yeo7_lh,yeo7_rh)
freq <- c(wsbm_lh, wsbm_rh)
data <- data.frame(as.character(yeo7),freq)
colnames(data) <- c('yeo7', 'freq')
data <- data[data$yeo7 != "0",]#remove the medial wall

## Raincloud plot
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

#Stats
kruskal.test(freq~yeo7, data = data)
pairwise.wilcox.test(data$freq, data$yeo7,
                     p.adjust.method = "bonferroni")
# Plot percentage of varying comm assignment (WSBM) by Yeo7 system-Supp Fig S1b ----------------------------------------

sumld<- ddply(data, ~yeo7, summarise, mean = mean(freq), median = median(freq), count = sum(freq!=0), total = length(freq), perc_varying_assign=sum(freq!=0)/length(freq))
g <- ggplot(data = sumld, aes(y = perc_varying_assign, x = yeo7, fill = yeo7)) +
  geom_bar(stat="identity", size = .5, alpha = 0.80) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme
g

# Plot variance in community assignment by WSBM systems-Supp Fig S8a ------------------------
#Read in the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
freq<- t(z$freq)
freq_continuous <- abs(freq-670) #change into the number of non-modal assignments out of 670
partitions <- readMat(paste0(wsbm_datadir,"n670_training_sample_consensus_partitions_yeorelabeled.mat"), drop = )
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled
data <- data.frame(as.character(consensus_iterative_labels),freq_continuous)
colnames(data) <- c('wsbm', 'freq')
#data <- data[data$freq != "0",]#remove the 0's in inconsistent assignment
## Raincloud plot
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

# Plot % of varying comm assignment (WSBM) by WSBM system-Fig 5b --------
g <- ggplot(data = sumld, aes(y = perc_varying_assign, x = wsbm, fill = wsbm)) +
  geom_bar(stat="identity", position = position_jitter(width = .08), size = .5, alpha = 0.80) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme
g

#Stats
kruskal.test(freq~wsbm, data = data)
pairwise.wilcox.test(data$freq, data$wsbm,
                     p.adjust.method = "bonferroni")

# Overlap of uncertainty in WSBM and Yeo Dev-Fig 5c ------------------------------
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"

#First relabel Schaefer400 annotation with WSBM frequency of inconsistent assignment in matlab
#Then import that annotation, compare both vectors in fsaverage6 space as we did before
#Inconsistent assignment WSBM
lh_freq <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.freq.fsaverage6' )
lh_freq <- as.numeric(as.character(lh_freq$label_names))
rh_freq <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.freq.fsaverage6' )
rh_freq <- as.numeric(as.character(rh_freq$label_names))

#Confidence for Yeo Dev, cut at 0 and mapping the 1's to the max value of confidence
silhouette_lh <-  ifelse(yeo_dev_partition$lh.s==1, max(yeo_dev_partition$lh.s[yeo_dev_partition$lh.s!=1]),yeo_dev_partition$lh.s) #this is in fsaverage6 space, so 40k vertices
silhouette_rh <-  ifelse(yeo_dev_partition$rh.s==1, max(yeo_dev_partition$rh.s[yeo_dev_partition$rh.s!=1]),yeo_dev_partition$rh.s) #this is in fsaverage6 space, so 40k vertices
silhouette_yeo7 <- c(silhouette_lh, silhouette_rh)

#want to make one go up and one go down in different colors, and where they overlap is both
assign_wsbm <- c(lh_freq, rh_freq)#this ends up on 0-0.5, but I want 0-1 for comparison with the silhouette measure
assign_wsbm_scaled <- (assign_wsbm-max(assign_wsbm))+max(assign_wsbm) #1 is highest variability in assignment, 0 is lowest
assign_wsbm_scaled <- assign_wsbm_scaled/max(assign_wsbm_scaled)

#scale silhouette so that the high confidence of each hemisphere is 0, lowest confidence (most negative) is 1
silhouette_scaled <- abs(silhouette_yeo7-max(silhouette_yeo7))
silhouette_scaled <- silhouette_scaled/max(silhouette_scaled) #lowest confidence is 1, not 1.3

#plot on brain to check
colormap = colorRampPalette(c("white","pink","blue"))
makecmap_options=list('colFn'=colormap)
vis.data.on.subject(subjects_dir, 'fsaverage6',assign_wsbm_scaled[1:40962], assign_wsbm_scaled[40963:81924],"inflated", makecmap_options = makecmap_options,  views="t4", rgloptions = rgloptions, draw_colorbar = T)

#continuous version
overlap <- silhouette_scaled+assign_wsbm_scaled
hist(assign_wsbm[yeo7!=0])
hist(silhouette_scaled[yeo7!=0])
hist(overlap[yeo7!=0])
rglactions=list("snapshot_png"=paste0(output_image_directory,"summed_symmetric_yeodev_confidence_inconsistency_wsbm.png"))
colormap = colorRampPalette(c("white", "#c48a69"))
makecmap_options=list('colFn'=colormap)
vis.data.on.subject(subjects_dir, 'fsaverage6',overlap[1:40962], overlap[40963:81924], "inflated", makecmap_options = makecmap_options,  views="t4", 
                    rgloptions = rgloptions, rglactions = rglactions, draw_colorbar = T)

# Overlap of silhouette measure in both partitions-Supp Fig S6b ------------------------
#need silhouette for WSBM vertex-wise
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/wsbm_consensus/"
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
subject_id = 'fsaverage6';       # for function which use one subject only
atlas='Schaefer2018_400Parcels_7Networks_order'

#Get the Schaefer atlas region names
schaefer_atlas_region_names_lh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="lh");
schaefer_atlas_region_names_rh = get.atlas.region.names(atlas, template_subjects_dir = subjects_dir,template_subject=subject_id, hemi="rh");

z <- readMat(paste0(wsbm_datadir,"silhouettes_wsbm.mat"), drop = )
silhouette_wsbm <- z$silhouette.wsbm.0.diag
silhouette_wsbm_lh=as.list(setNames(c(NA, silhouette_wsbm[1:200]), schaefer_atlas_region_names_lh))
silhouette_wsbm_rh=as.list(setNames(c(NA, silhouette_wsbm[201:400]), schaefer_atlas_region_names_rh))
#write it out to an atlas file with the regions
# write.region.values(subjects_dir,'fsaverage6', 'lh', 'Schaefer2018_400Parcels_7Networks_order', silhouette_wsbm_lh,
#                     'silhouette_wsbm_fsaverage6', output_path=output_image_directory)
# write.region.values(subjects_dir,'fsaverage6', 'rh', 'Schaefer2018_400Parcels_7Networks_order', silhouette_wsbm_rh,
#                     'silhouette_wsbm_fsaverage6', output_path=output_image_directory)
#Confidence for Yeo Dev, mapping the 1's to the max value of confidence
silhouette_lh <-  ifelse(yeo_dev_partition$lh.s==1, max(yeo_dev_partition$lh.s[yeo_dev_partition$lh.s!=1]),yeo_dev_partition$lh.s) #this is in fsaverage6 space, so 40k vertices
silhouette_rh <-  ifelse(yeo_dev_partition$rh.s==1, max(yeo_dev_partition$rh.s[yeo_dev_partition$rh.s!=1]),yeo_dev_partition$rh.s) #this is in fsaverage6 space, so 40k vertices
silhouette_yeo7 <- c(silhouette_lh, silhouette_rh)

#read it in
silhouette_wsbm_lh <- read.fs.mgh(paste0(output_image_directory,'lh.silhouette_wsbm_fsaverage6.mgz'), flatten = T)
silhouette_wsbm_rh <- read.fs.mgh(paste0(output_image_directory,'rh.silhouette_wsbm_fsaverage6.mgz'), flatten = T)
silhouette_wsbm_lh[is.na(silhouette_wsbm_lh)] <- max(silhouette_wsbm_lh[!is.na(silhouette_wsbm_lh)])
silhouette_wsbm_rh[is.na(silhouette_wsbm_rh)] <- max(silhouette_wsbm_rh[!is.na(silhouette_wsbm_rh)])

#want to make one go up and one go down in different colors, and where they overlap is both
silhouette_wsbm <- c(silhouette_wsbm_lh, silhouette_wsbm_rh)#this ends up on -0.64-0.79, but I want 0-1 for comparison with the silhouette measure
silhouette_wsbm_scaled <- abs(silhouette_wsbm-max(silhouette_wsbm))
silhouette_wsbm_scaled <- silhouette_wsbm_scaled/max(silhouette_wsbm_scaled)

#scale silhouette so that the high confidence of each hemisphere is 0, lowest confidence (most negative) is 1
silhouette_scaled <- abs(silhouette_yeo7-max(silhouette_yeo7))
silhouette_scaled <- silhouette_scaled/max(silhouette_scaled) #lowest confidence is 1, not 1.3

#continuous version
overlap <- silhouette_scaled+silhouette_wsbm_scaled
hist(silhouette_wsbm_scaled[yeo7!=0])
hist(silhouette_scaled[yeo7!=0])
hist(overlap[yeo7!=0])
rglactions=list("snapshot_png"=paste0(output_image_directory,"summed_symmetric_yeodev_confidence_wsbm_confidence.png"))
colormap = colorRampPalette(c("white", "#ab7150"), space = 'rgb') #one shade darker
makecmap_options=list('colFn'=colormap)
vis.data.on.subject(subjects_dir, 'fsaverage6',overlap[1:40962], overlap[40963:81924], "inflated", makecmap_options = makecmap_options,  views="t4", 
                    rgloptions = rgloptions, draw_colorbar = T)

silhouette_wsbm_lh_norm<- (silhouette_wsbm_lh-min(silhouette_wsbm_lh))/(max(silhouette_wsbm_lh)-min(silhouette_wsbm_lh))
silhouette_wsbm_rh_norm<- (silhouette_wsbm_rh-min(silhouette_wsbm))/(max(silhouette_wsbm)-min(silhouette_wsbm))
silhouette_yeodev_lh_norm<- (silhouette_lh-min(silhouette_lh))/(max(silhouette_lh)-min(silhouette_lh))
silhouette_yeodev_rh_norm<- (silhouette_rh-min(silhouette_rh))/(max(silhouette_rh)-min(silhouette_rh))

overlap_lh <- silhouette_yeodev_lh_norm+silhouette_wsbm_lh_norm
overlap_rh <- silhouette_yeodev_rh_norm+silhouette_wsbm_rh_norm

#anything above 0.82 is uncertain in both partitions!
rgloptions=list("windowRect"=c(50,50,1000,1000));

# SNR: Supplemental Figure S3b ------------------------------------------------
output_image_directory="/cbica/projects/spatial_topography/output/images/brains/"
#Read in adult SNR
snr_dir="~/Documents/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/"
l <- read_nisurfacefile(paste0(snr_dir,"lh.1000subjects.snr.nii.gz"))
l <- readNifti(paste0(snr_dir,"lh.1000subjects.snr.nii.gz")) #there is some issue with reading this in...
#this is the only function that works and doesn't give an error
adult_snr_lh <- l[1:length(l)]
l <- readNifti(paste0(snr_dir,"rh.1000subjects.snr.nii.gz")) #there is some issue with reading this in...
adult_snr_rh <- l[1:length(l)]
#mask out medial wall by reading in yeo7 in fsaverage5 and removing those vertices
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";
lh <- subject.annot(subjects_dir, 'fsaverage5', 'lh','Yeo2011_7Networks_N1000' );yeo7_lh <- as.numeric(str_split(lh$label_names, "_", simplify = TRUE)[,2])
adult_snr_lh[is.na(yeo7_lh)] <- 0 #remove the medial wall
rh <- subject.annot(subjects_dir, 'fsaverage5', 'rh','Yeo2011_7Networks_N1000' );yeo7_rh <- as.numeric(str_split(rh$label_names, "_", simplify = TRUE)[,2])
adult_snr_rh[is.na(yeo7_rh)] <- 0 #remove the medial wall

#Plot adult SNR
colormap= grDevices::colorRampPalette(c("white",RColorBrewer::brewer.pal(11, name="YlOrRd")));
#colormap = colorRampPalette(c("white",rev(rainbow(10)))) #check this against yeo image, freeview projection
makecmap_options=list('colFn'=colormap)
rglactions=list("snapshot_png"=paste0(output_image_directory,"GSP_1000subjects_SNR.png"))
vis.data.on.subject(subjects_dir, 'fsaverage5',adult_snr_lh, adult_snr_rh, "inflated", makecmap_options = makecmap_options,draw_colorbar = T, views="t4", rglactions = rglactions, rgloptions = rgloptions)

#Read in child SNR
snr_dir="/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/surfaces/"
l <- read_nisurfacefile(paste0(snr_dir,"lh.fs5_average_tsnr_avg_2runs_trunc_volsdropped.nii.gz"))
l <- readNifti(paste0(snr_dir,"lh.fs5_average_tsnr_avg_2runs_trunc_volsdropped.nii.gz")) #there is some issue with reading this in...
#this is the only function that works and doesn't give an error
child_snr_lh <- l[1:length(l)]
l <- readNifti(paste0(snr_dir,"rh.fs5_average_tsnr_avg_2runs_trunc_volsdropped.nii.gz")) #there is some issue with reading this in...
child_snr_rh <- l[1:length(l)]
#mask out medial wall
child_snr_lh[is.na(yeo7_lh)] <- 0 #remove the medial wall
child_snr_rh[is.na(yeo7_rh)] <- 0 #remove the medial wall

#Plot child SNR
colormap= grDevices::colorRampPalette(c("white",RColorBrewer::brewer.pal(5, name="YlOrRd")));
#colormap = colorRampPalette(c("white",rev(rainbow(10)))) #check this against yeo image, freeview projection
makecmap_options=list('colFn'=colormap)
vis.data.on.subject(subjects_dir, 'fsaverage5',child_snr_lh, child_snr_rh, "inflated", makecmap_options = makecmap_options,draw_colorbar = T, views="t4", rgloptions = rgloptions)

#make colormap range same as adult, which is 0-96.34
child_snr_lh[1333] <- 96.34887 #make a random, non-visible vertex the max value
colormap= grDevices::colorRampPalette(c("white",RColorBrewer::brewer.pal(11, name="YlOrRd")));
makecmap_options=list('colFn'=colormap)
rglactions=list("snapshot_png"=paste0(output_image_directory,"ABCD_670subjects_SNR.png"))
vis.data.on.subject(subjects_dir, 'fsaverage5',child_snr_lh, child_snr_rh, "inflated", makecmap_options = makecmap_options,draw_colorbar = T, views="t4",rglactions = rglactions, rgloptions = rgloptions)

# Overlay SNR map on dev communities: Supp Fig S3b ------------------------------------------
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','yeodev.fsaverage6' )
rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','yeodev.fsaverage6' )
#Read in child SNR in fsaverage6
snr_dir="/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/surfaces/"
child_snr_lh <- read.fs.curv("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/freesurfer/fsaverage6/surf/lh.fs6_average_tsnr_avg_2runs_trunc_volsdropped.curv") #there is some issue with reading this in...
child_snr_rh  <- read.fs.curv("/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/freesurfer/fsaverage6/surf/rh.fs6_average_tsnr_avg_2runs_trunc_volsdropped.curv") #there is some issue with reading this in...

child_snr_lh[lh$label_names=="NONAME0"] <- 0 #remove the medial wall
child_snr_rh[rh$label_names=="NONAME0"] <- 0 #remove the medial wall
#read in fs6 surface
fs6_lh <- read.fs.surface(paste0(subjects_dir, "fsaverage6/surf/lh.inflated"))
fs6_rh <- read.fs.surface(paste0(subjects_dir, "fsaverage6/surf/rh.inflated"))
#outline the yeodev communities
yeodev_comm_outlines_lh <- annot.outline(lh,fs6_lh, outline_color = "black")
yeodev_comm_outlines_rh <- annot.outline(rh,fs6_rh, outline_color = "black") # limit_to_regions = "7Networks_6")
#overlay with snr
yeodev_snr_outlines_lh <- ifelse(yeodev_comm_outlines_lh=="black",0, child_snr_lh)
yeodev_snr_outlines_rh <- ifelse(yeodev_comm_outlines_rh=="black",0, child_snr_rh)
child_snr_lh[1333] <- 96.34887 #make a random, non-visible vertex the max value
colormap= grDevices::colorRampPalette(c("white",RColorBrewer::brewer.pal(11, name="YlOrRd")));
makecmap_options=list('colFn'=colormap)
rglactions=list("snapshot_png"=paste0(output_image_directory,"ABCD_SNR_comm_outlines.png"))
vis.data.on.subject(subjects_dir, 'fsaverage6',yeodev_snr_outlines_lh,yeodev_snr_outlines_rh, "inflated", makecmap_options = makecmap_options,
                    views="t4", rgloptions = rgloptions, rglactions = rglactions)

# Violin plot of adult SNR by adult community: Supp Fig S3a -----------------------------
#read in fsaverage6 space
adult_snr_rh <- read.fs.curv("~/Documents/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/rh.fs6.1000subjects.snr.nii.gz")
adult_snr_lh <- read.fs.curv("~/Documents/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/lh.fs6.1000subjects.snr.nii.gz")

adult_snr <- c(adult_snr_lh,adult_snr_lh) 
yeo7 <- c(yeo7_lh, yeo7_rh)
data <- data.frame(adult_snr, as.character(yeo7)) #in fsaverage5 space
#data <- data[is.na(data$as.character.yeo7.),]#remove the medial wall
colnames(data) <- c("snr", "yeo7")
data <- data %>% filter(snr!=0) %>% filter(!is.na(yeo7)) #remove the 0's for snr on medial wall
g <- ggplot(data = data, aes(y = snr, x = yeo7, fill = yeo7)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = snr, color = yeo7), position = position_jitter(width = .15), size = .08, alpha = 0.15) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme +
  coord_cartesian(ylim=c(5,95))
g

# Violin plot of dev SNR by dev community: Supp Fig S3b ---------------------------------
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("~/Documents/tools/raincloud.R")

#mask out medial wall by reading in yeo7 in fsaverage6 and removing those vertices
subjects_dir = "/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/";
lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','Yeo2011_7Networks_N1000' );yeo7_lh <- as.numeric(str_split(lh$label_names, "_", simplify = TRUE)[,2])
rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','Yeo2011_7Networks_N1000' );yeo7_rh <- as.numeric(str_split(rh$label_names, "_", simplify = TRUE)[,2])

child_snr <- c(child_snr_lh,child_snr_rh) 
data <- data.frame(child_snr, as.character(yeo_dev_full))
data <- data[data$as.character.yeo_dev_full. != "0",]#remove the medial wall
colnames(data) <- c("snr", "yeodev")
#data <- data %>% filter(silhouette!=1) #remove the 1's for silhouette around the medial wall

## Raincloud plot
g <- ggplot(data = data, aes(y = snr, x = yeodev, fill = yeodev)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
  geom_point(aes(y = snr, color = yeodev), position = position_jitter(width = .15), size = .08, alpha = 0.15) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  scale_fill_manual(values= c("#7B287E", "#5CA1C8", "#0A9045","#C33AF8", "#B6C988", "#EF9C23", "#E34A53")) +
  # coord_flip() +
  theme_bw() +
  raincloud_theme +
  coord_cartesian(ylim=c(5,95))
g

#Stats
kruskal.test(silhouette~yeo7, data = data)
pairwise.wilcox.test(data$silhouette, data$yeo7,
                     p.adjust.method = "bonferroni")

# Stats for overleaf output -----------------------------------------------
#Load ABCD demographic data
data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
output_data_dir='~/Documents/projects/in_progress/spatial_topography_CUBIC/data/subjData/Release2_fixed/'
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/'
net_data_dir='/cbica/projects/spatial_topography/data/imageData/net_stats/'
raw_data_dir='/data/jux/mackey_group/public_data/ABCD/'
sites <- read.delim(paste0(data_dir,"abcd_lt01.txt"), stringsAsFactors = F)# change strings as factors = FALSE next time
sites = sites[-1,]
sites = droplevels(sites)
sites <- select(sites, -c(collection_id:dataset_id))
sites$ID <- sites$subjectkey
socio <- read.delim(paste0(data_dir,"acspsw03.txt"), stringsAsFactors = F)
socio = socio[-1,]
socio = droplevels(socio)
socio <- select(socio, -c(collection_id:dataset_id)) 
socio$ID <- socio$subjectkey
income <- read.delim(paste0(data_dir,"pdem02.txt")) 
income = income[-1,]
income = droplevels(income)
income <- select(income, -c(collection_id:dataset_id))
income$ID <- income$subjectkey

main<- left_join(sites,socio, by=c("ID", "eventname"))
main <- left_join(main, income, by=c("ID", "eventname"))

#Only baseline visits from cognition
main <- filter(main, eventname=="baseline_year_1_arm_1") %>% filter(.,site_id_l=="site16") #filter out only the baseline visits
#for site 16 only
subjlist1 <- read.table(paste0(subjlist_dir, "/n670_subjects_only_filtered_runs_site16_postprocess.txt"), col.names = c("ID"))
#Keep only n670 subject list
main$ID <- str_remove(main$ID, "_") #subjectid has a _ in it after NDAR in the NDAR data, doesn't match subject list
main$ID <- paste0("sub-",main$ID)
#make age numeric
main$age <- as.numeric(main$interview_age.x)/12
#make sex a factor
main$gender <- as.factor(main$gender.x)

main_yeodev <- left_join(subjlist1, main, by="ID")

#Education
main_yeodev$demo_prtnr_ed_v2
#parent edu variable
summary(main_yeodev$demo_prnt_ed_v2)
main_yeodev$demo_prnt_ed_num_v2 <- as.numeric(main_yeodev$demo_prnt_ed_v2)
main_yeodev$demo_prnt_ed_num_v2_recoded <- car::recode(as.character(main_yeodev$demo_prnt_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_yeodev$demo_prnt_ed_num_v2_recoded <- as.numeric(main_yeodev$demo_prnt_ed_num_v2_recoded)

#partner edu variable
summary(main_yeodev$demo_prtnr_ed_v2)
main_yeodev$demo_prtnr_ed_num_v2 <- as.numeric(main_yeodev$demo_prtnr_ed_v2)
main_yeodev$demo_prtnr_ed_num_v2_recoded <- car::recode(as.character(main_yeodev$demo_prtnr_ed_v2), "0=0; 1=1; 2=2; 3=3; 4=4; 5=5; 6=6; 7=7; 8=8; 9=9; 10=10; 11=11; 12=12; 13=12; 14=12; 15=13; 16=14; 17=14; 18=16; 19=18; 20=20; 21=20; 777= NA; 999 = NA")
main_yeodev$demo_prtnr_ed_num_v2_recoded <- as.numeric(main_yeodev$demo_prtnr_ed_num_v2_recoded)

#creating a combined parent education variable
main_yeodev$demo_comb_parent_edu <- rowMeans(main_yeodev[c('demo_prnt_ed_num_v2_recoded', 'demo_prtnr_ed_num_v2_recoded')], na.rm=TRUE)
View(main_yeodev$demo_comb_parent_edu)
hist(as.numeric(main_yeodev$demo_comb_parent_edu), xlab = "Combined averaged years of parent education", breaks = 15, label = TRUE, col = "red")
median(main_yeodev$demo_comb_parent_edu)

#Is motion associated with age?
cor.test(main_yeodev$age, main_yeodev$fd_mean_avg)
cor.test(main_yeodev$age, main_yeodev$pctVolsCensored)

save(main_yeodev, file="~/Downloads/Demo_data.Rdata")
load("Demo_data.Rdata")

wsbm_full <- c(wsbm_lh,wsbm_rh)
save(yeo_dev_full, yeo_full, wsbm_full,NMI_yeodev_yeo7, NMI_wsbm_yeo7, NID_wsbm_yeo7, NID_yeodev_yeo7, file="~/Downloads/Partition_comparison_stats.RData")

#NMI
NMI_yeodev_yeo7<- compare(yeo_dev_full, yeo_full, method = "nmi")
NID_yeodev_yeo7<- clustComp(yeo_dev_full, yeo_full)$NID

NMI_wsbm_yeo7<- compare(wsbm_full, yeo_full, method = "nmi")
NID_wsbm_yeo7<- clustComp(wsbm_full, yeo_full)$NID

#Confidence & variability
silhouette_yeo <- c(silhouette_lh_fs6,silhouette_rh_fs6)
silhouette_dev <- c(silhouette_lh_dev,silhouette_rh_dev)
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
wsbm_rh <- subject.annot(subjects_dir, 'fsaverage6', 'rh','wsbm.freq.fsaverage6');wsbm_rh<- as.numeric(as.character(wsbm_rh$label_names))
wsbm_lh <- subject.annot(subjects_dir, 'fsaverage6', 'lh','wsbm.freq.fsaverage6');wsbm_lh<- as.numeric(as.character(wsbm_lh$label_names))
freq <- c(wsbm_lh, wsbm_rh)

#WSBM variability by parcel
yeo_nodes <- as.character(yeo_nodes$yeo_nodes)
#Read in the frequency of ties in the group WSBM partition
z <- readMat(paste0(wsbm_datadir,"consensus_iter_mode.mat"), drop = )
freq<- t(z$freq)
consensus_iterative_labels <- partitions$consensus.iter.mode.yeorelabeled

save(silhouette_dev, yeo_dev_full, silhouette_yeo, freq, yeo_nodes,consensus_iterative_labels, file="~/Downloads/confidence_variability.RData")
