#!/bin/csh

set orig_data_dir = "${ABCD_DATA}/bids_release2_site16/derivatives/surfaces"
set output_dir = $1
set sub_list = "/cbica/projects/spatial_topography/dropbox/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/yeo_networks/n2_subjlist_test_CUBIC.txt"
set surf_list = "/cbica/projects/spatial_topography/dropbox/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/yeo_networks/n2_fullpaths_test_CUBIC.csv"
set sub_dir = "${ABCD_DATA}/bids_release2_site16/derivatives/surfaces"
set subjects = `cat $sub_list`
set code_dir = "${CBIG_CODE_DIR}/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering"
set out_dir = "/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/"
set surf_stem = "fs6_sm4.8_task-rest_uncensor_trunc_reshape"
set outlier_stem = "_FDRMS0.2_DVARS50_motion_outliers"
set target = fsaverage6
set roi = fsaverage3

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
${code_dir}/CBIG_Yeo2011_compute_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list} -surf_ls ${surf_list} -target $target -roi $roi

#Call the clustering function to cluster FC profiles for each subject
#${code_dir}/CBIG_Yeo2011_cluster_fcMRI_surf2surf_profiles_subjectlist.csh -sd ${sub_dir} -sub_ls ${sub_list} -n ${num_clusters} -out ${cluster_out} -tries ${num_tries} -mesh $target
