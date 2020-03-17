%% PATHS
nets_dir='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries'
data_dir='/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek'
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks'

file='yeo7_n670_2runsonly_1000tries.mat'
load(fullfile(nets_dir, file))
clustered=load(fullfile(nets_dir, file))

my_colors = [190   190  190
160    32  240
0     0  255
127   255  212
255   246  143
255     0  255
255   165    0
255     0    0]
%% write Yeo-dev community partitition to an annotation or label file

CBIG_generate_ROIlabel_from_parcellation(lh_labels',rh_labels', 'test_ROIlabel_fsaverage6', 'fsaverage6', nets_dir)
%this is not what I want because these are each separate files for a label!
CBIG_SaveParcellationToFreesurferAnnotation(fullfile(nets_dir, file), 'lh.nets.fsaverage.annot', 'rh.nets.fsaverage.annot')
%this uses the old colors

CBIG_WriteParcellationToAnnotation(lh_labels,fullfile(nets_dir,'lh.yeonets.fsaverage6.annot'), my_colors)
CBIG_WriteParcellationToAnnotation(rh_labels,fullfile(nets_dir,'rh.yeonets.fsaverage6.annot'), my_colors)
%this is what I want! With my colors and in fsaverage6 space.

%% get subject correlation matrices based on this partition
yeo_dev = fullfile(nets_dir,'lh.yeonets.fsaverage6.annot')
subjlist1 = fullfile(subjlist_dir, 'n3_site16_tworunsonly_volume_filepaths.txt')
motionlist = fullfile(subjlist_dir, 'n3_site16_tworunsonly_motion_filepaths.txt')

% CBIG_ComputeROIs2ROIsCorrelationMatrix(output_file, subj_text_list1, subj_text_list2, discard_frames_list, ROIs1, ROIs2, regression_mask1, regression_mask2, all_comb_bool, avg_sub_bool)
CBIG_ComputeROIs2ROIsCorrelationMatrix('testoutputfile.mat', subjlist1, subjlist1 , 'NONE', yeo_dev,yeo_dev,'NONE','NONE', 1, 0)
%corr mat is file of 7 x 7 correlation matrices

surflist1 = dlmrefullfile(subjlist_dir,'n1_surf2surf_lh.txt')
surflist2 = fullfile(subjlist_dir,'n1_surf2surf_rh.txt')
%CBIG_ComputeFullSurfaceCorrelation('testoutputfullsurf.mat' , varargin_text1, varargin_text2, pval)
CBIG_ComputeFullSurfaceCorrelation('testoutputfullsurf.mat' , surflist1, surflist2, 0)

%% zscore and average together