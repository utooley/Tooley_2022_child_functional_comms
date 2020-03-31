%% PATHS
nets_dir='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries'
data_dir='/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek'
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks'

file='yeo7_n670_2runsonly_1000tries.mat'
load(fullfile(nets_dir, file))
clustered=load(fullfile(nets_dir, file))
labels=[clustered.lh_labels
    clustered.rh_labels]
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
CBIG_ComputeROIs2ROIsCorrelationMatrix('testoutputfile.mat', subjlist1, subjlist1 , motionlist, yeo_dev,yeo_dev,'NONE','NONE', 1, 0)
%corr mat is file of 7 x 7 correlation matrices

surflist1 = readtable(fullfile(subjlist_dir,'n2_surf2surf_lh.csv'), 'ReadVariableNames', false, 'Delimiter', ' ')
surflist2 = readtable(fullfile(subjlist_dir,'n2_surf2surf_rh.csv'), 'ReadVariableNames', false, 'Delimiter', ' ')
surlist1 = fullfile(subjlist_dir,'n2_surf2surf_lh2.txt')
surflist2=fullfile(subjlist_dir,'n2_surf2surf_rh2.txt')
%CBIG_ComputeFullSurfaceCorrelation('testoutputfullsurf.mat' , varargin_text1, varargin_text2, pval)
CBIG_ComputeFullSurfaceCorrelation('testoutputfullsurf.mat' , surflist1, surflist2, 0)

Ci=labels;
nCi = unique(Ci);
M=corr_mat;

for i = 1:length(nCi) % loop through communities
    for j = 1:length(nCi)
       Wi = Ci == nCi(i); % find index for within communitiy edges
       Bi = Ci == nCi(j); % find index for between communitiy edges to specific community

       Wv_temp = M(Wi,Wi); % extract within communitiy edges
       Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community

       %calculate standard deviation of edges within blocks here
       %deviation_edge_weights_yeo(n,1)=deviation_edge_weights_yeo(n,1)+std(Wv_temp(Wv_temp~=0))^2; %Gu paper calculation, finished below
       %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
       Bv = [Bv_temp(:)'];
       system_connectivity(i,j)=mean(Bv(Bv~=0));
       %system_between(i,1)=mean(Bv);
       %if i==j
       %else
    end
end

%% zscore and average together