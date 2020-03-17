#!/bin/csh

#run this with qsub -l h_vmem=50.0G,s_vmem=49.7G -j y -o ${SPATIAL_TOPOGRAPHY}/output/job_output/ -cwd /cbica/projects/spatial_topography/code/yeo_networks/yeo_clustering_n670_1000x_mot_outliers.csh

set output_dir = "/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries_mot_outliers_2" #change this!
set cluster_out = "${output_dir}/yeo7_n670_2runsonly_1000tries_mot_outliers_2" #CHANGE THIS!
set orig_data_dir = "${ABCD_DATA}/bids_release2_site16/derivatives/surfaces"
set sub_list = "/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n670_site16_subjlist_tworunsonly.txt"
set surf_list = "/cbica/projects/spatial_topography//data/subjLists/release2/site16/yeo_networks/n670_site16_fullpaths_tworunsonly.csv"
set sub_dir = "${ABCD_DATA}/bids_release2_site16/derivatives/surfaces"
set subjects = `cat $sub_list`
set code_dir = "${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering"
set surf_stem = "fs6_sm4.8_task-rest_uncensor_trunc_reshape"
set outlier_stem = "FDRMS0.2_DVARS50_motion_outliers"
set outlier_list = "/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n670_site16_tworunsonly_motion_filepaths.csv"
set target = fsaverage6
set roi = fsaverage3
set num_clusters = 7
set num_tries = 1000
set scrub_flag = 1

## Create folder structure within output_dir, and make soft links of input files to orig_data_dir
foreach s ($subjects)
	set s_id = `echo $s | cut -d '_' -f 1`
	set sess_id = `echo $s | cut -d '_' -f 2`
	mkdir -p $output_dir/subjects/$s
#	ln -s $orig_data_dir/$s_id/$s/qc $output_dir/subjects/$s/
	ln -s $orig_data_dir/$s_id/surf $output_dir/subjects/$s/
#	ln -s $orig_data_dir/$s_id/logs $output_dir/subjects/$s/
end


## Make a directory
mkdir -p $output_dir/clustering
#skip the creation of the subject list since I already have them

#call the compute function for FC profiles for each subject
#add outlier list with 	-outlier_ls ${outlier_list}
${code_dir}/CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list} -surf_ls ${surf_list} -outlier_ls ${outlier_list} -target $target -roi $roi

#Call the clustering function to cluster FC profiles for each subject
#add -scrub_flag 1 to include the motion outliers
${code_dir}/CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list} -n ${num_clusters} -out ${cluster_out} -tries ${num_tries} -scrub_flag ${scrub_flag} -mesh $target -roi ${roi}
