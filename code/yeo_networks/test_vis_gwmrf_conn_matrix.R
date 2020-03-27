library(R.matlab)

# SETUP -------------------------------------------------------------------
#CUBIC cluster filepaths, mounted locally
wsbm_datadir="/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/"
subjects_dir = "/Users/utooley/Documents/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/";
yeo_dev_dir="/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_mot_outliers/"
yeo7_ref_dir="/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/"
gwmrf_dir="/cbica/projects/spatial_topography/data/imageData/gwMRF_schaefer400_site16_n670/mult_mat/"
# Read in files -----------------------------------------------------------
#Vector of WSBM community assignments
avg_matrix <- read.csv(paste0(gwmrf_dir,"/averaged_cov_mat_from_gwmrf_n670.csv"), header = FALSE)


avgmatrix <- as.matrix(avg_matrix)
# palf <- colorRampPalette(c("red","white", "blue")) 
# heatmap(avgmatrix, Rowv = NA, Colv = NA, col = palf(20))
plot(vizu.mat(avgmatrix, fill.scale.limits = c(-1.3, 1.3), x.lab=pipeline))