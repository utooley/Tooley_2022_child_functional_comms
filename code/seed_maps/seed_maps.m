%% Seed maps
%Need to edit CBIG_ComputeSeedBasedCorrelation such that it takes the whole
%surface as a mask, besides the medial wall

%Need to transform seed first, and make a list of volumetric seeds for each
%subject

%%%% LH %%%%
seed_file_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n670_vols_site16_fullpaths_tworunsonly.txt" % is list of files that are time series to be used as seed regions
target_file_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n670_site16_fullpaths_tworunsonly_lh.txt" %is list of files that are time series to be used as target regions
seed_mask_file='/cbica/projects/spatial_topography/tools/parcellations/Yeo_ThalamusParcellation/thalamus_5_24mm_bin0.5.nii.gz'  %is a single file that is a binary mask defining the seed
target_mask_file='' %is a single file that is a binary mask defining the region of target time series to perform corrrelation
output_file='/cbica/projects/spatial_topography/data/imageData/seed_analyses/lh.xxx.nii.gz'
type='average'

%this does what I want, but only does one set of subject hemispheres at a
%time! That's fine, I can just read in all the lh and then all the rh and
%write them out separately.
CBIG_ComputeSeedBasedCorrelation(seed_file_txt, target_file_txt, seed_mask_file, target_mask_file, type, output_file)

%%%% RH %%%%%
seed_file_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_vol_fullpaths_tworunsonly.txt" % is list of files that are time series to be used as seed regions
target_file_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_fullpaths_tworunsonly_rh.txt" %is list of files that are time series to be used as target regions
seed_mask_file='/cbica/projects/spatial_topography/tools/parcellations/Yeo_ThalamusParcellation/thalamus_5_24mm_bin0.5.nii.gz'  %is a single file that is a binary mask defining the seed
target_mask_file='' %is a single file that is a binary mask defining the region of target time series to perform corrrelation
output_file='/cbica/projects/spatial_topography/data/imageData/seed_analyses/rh.xxx.nii.gz'

CBIG_ComputeSeedBasedCorrelation(seed_file_txt, target_file_txt, seed_mask_file, target_mask_file, type, output_file)



%% Testing

%compute seed-based correlation?
roi_file='/cbica/projects/spatial_topography/tools/parcellations/Yeo_ThalamusParcellation/thalamus_5_24mm_bin0.5.nii.gz'
output_file='/cbica/projects/spatial_topography/data/imageData/test.nii'
surf_input_files_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_fullpaths_tworunsonly.txt"
vol_input_files_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_vol_fullpaths_tworunsonly.txt"

%save subject lists as Tab-delimited text.
%Need to upsample the thalamus parcellation to 2.4 mm isotropic, which is
%ABCD voxel size, then binarize.
%Also need to take either the subject volume data in MNI space, or
%transform this parcel to subject space for each subject?

CBIG_ComputeCorrelationProfileVol2Surf(roi_file, output_file, 1, surf_input_files_txt, vol_input_files_txt)

output_dir="/cbica/projects/spatial_topography/data/imageData/"   %output file name
seed_mask_file='/cbica/projects/spatial_topography/tools/parcellations/Yeo_ThalamusParcellation/thalamus_5_24mm_bin0.5.nii.gz'      %surface mask
lh_varargin_text="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_fullpaths_tworunsonly.txt"
rh_varargin_text="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_fullpaths_tworunsonly_rh.txt" %a text file containing input surface data list of left and right hemisphere.
vol_varargin_text="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_vol_fullpaths_tworunsonly.txt"  %a text file containing volumetric data list
regress_bool=0 %   negative or 0 for no regression; 1 for Cerebellum and Striatum regression; 2 for Putamen regression
start_index='1'
stop_index='10000' %only consider indices from start_index to stop_index that have not been already produced
pval='0.01' %If pval > 0 then this function will compute statistics

%This one does something weird with taking the seed-based corr of each
%index separately.
CBIG_ComputeFullVol2SurfCorrelationIndex(output_dir, lh_varargin_text, rh_varargin_text, vol_varargin_text, seed_mask_file, regress_bool, start_index, stop_index, pval)

mask_file='/cbica/projects/spatial_topography/tools/parcellations/Yeo_ThalamusParcellation/thalamus_5_24mm_bin0.5.nii.gz'      %surface mask
output_file="/cbica/projects/spatial_topography/data/imageData/test"   %output file name, a matlab corr matrix.

%this one does what I want eventually, but outputs a correlation
%matrix that is num of seed voxels x number of surface vertices
% so then you have to average across the
%correlation with each of the voxels within the seed mask after it's
%calculated, which is potentially different than the average correlation
%with the whole seed timeseries
CBIG_ComputeFullVol2SurfCorrelation(output_file, lh_varargin_text, rh_varargin_text, vol_varargin_text, mask_file, regress_bool, pval) 

seed_file_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_vol_fullpaths_tworunsonly.txt" % is list of files that are time series to be used as seed regions
target_file_txt="/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks/n2_test_site16_fullpaths_tworunsonly.txt" %is list of files that are time series to be used as target regions
seed_mask_file='/cbica/projects/spatial_topography/tools/parcellations/Yeo_ThalamusParcellation/thalamus_5_24mm_bin0.5.nii.gz'  %is a single file that is a binary mask defining the seed
target_mask_file='' %is a single file that is a binary mask defining the region of target time series to perform corrrelation
output_file='/cbica/projects/spatial_topography/data/imageData/test.nii.gz'
output_file='test.nii.gz' 
type='average'

%this does what I want, but only does one set of subject hemispheres at a
%time! That's fine, I can just read in all the lh and then all the rh and
%write them out separately.
CBIG_ComputeSeedBasedCorrelation(seed_file_txt, target_file_txt, seed_mask_file, target_mask_file, type, output_file)
