%x='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/yeo7_n670_2runsonly_1000tries_lh_profile.txt'
addpath('/cbica/projects/spatial_topography/code/yeo_networks')
%input profiles data.
profile1='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/lh.yeo7_n670_2runsonly_1000tries.avg_profiles007.nii.gz'
profile2='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/rh.yeo7_n670_2runsonly_1000tries.avg_profiles007.nii.gz'
mesh_name='fsaverage6'
mask='cortex'
num_smooth=0;
normalize=0;

% %to get data into the right matrix form may need below! This is pulled from
% %CBIG_VonMisesSeriesClustering_fix_bessel_randnum_bsxfun.m
% 
% % read data (voxels x N subjects)
% series = read_fmri(profile1);
% if(~isempty(strfind(mesh_name, 'fsaverage')))
%     series2 = read_fmri(profile2);
%     series = [series; series2];
% end
% % read mask
% if(~isempty(strfind(mesh_name, 'fsaverage')))
%     lh_avg_mesh = CBIG_ReadNCAvgMesh('lh', mesh_name, 'inflated', mask);
%     l1 = find(lh_avg_mesh.MARS_label == 2); lh_num_verts = size(lh_avg_mesh.vertices, 2);
%     rh_avg_mesh = CBIG_ReadNCAvgMesh('rh', mesh_name, 'inflated', mask);
%     l2 = find(rh_avg_mesh.MARS_label == 2); rh_num_verts = size(rh_avg_mesh.vertices, 2);
%     l = [l1 l2+length(lh_avg_mesh.MARS_label)];
% else
%     % mesh.l: 0 - medial wall, 1 - cortex
%     [mesh.v, mesh.l, mesh.ct] = read_annotation(fullfile(getenv('CBIG_CODE_DIR'), 'data', 'templates', 'surface', 'fs_LR_32k', 'label', 'medialwall.annot'));
%     lh_num_verts = length(mesh.l) / 2;
%     rh_num_verts = lh_num_verts;
%     cort_label = mesh.ct.table(2, 5);
%     l1 = 1:lh_num_verts;          l1 = l1(mesh.l(l1)==cort_label);
%     l2 = 1:rh_num_verts;          l2 = l2(mesh.l(l2+lh_num_verts)==cort_label);
%     l = [l1 l2+lh_num_verts];
% end
% 
% % We (Hesheng, Mert, Jorge and I) decided that random initialization can be same across groups of subjects.
% % rand('twister',5489) was moved to after MRIread is because MRIread calls/apps/arch/Linux_x86_64/freesurfer/4.5.0/matlab/load_nifti.m, which calls rand('state', sum(100*clock));
% rand('twister',5489)
% 
% % smooth, only applied for data in fsaverage* space
% if(~isempty(strfind(mesh_name, 'fsaverage')))
%     series(1:end/2, :)     = transpose(MARS_AverageData(lh_avg_mesh, transpose(series(1:end/2, :)), 0, num_smooth));
%     series(end/2+1:end, :) = transpose(MARS_AverageData(rh_avg_mesh, transpose(series(end/2+1:end, :)), 0, num_smooth));
% end
% 
% % extract mask voxels series
% series = series(l, :);
% 
% % remove voxels that are completely not correlated with any rois. 
% non_zero_corr_index = (sum(series, 2) ~= 0);
% series = series(non_zero_corr_index, :);
% % znormalize (series assumed to be voxels x subjects or voxels x profile)
% % can take out as don't normalize
% % if(normalize)
% %     mean_series = nanmean(series, 1);
% %     std_series = nanstd(series, 1, 1);
% %     series = bsxfun(@minus, series, mean_series);
% %     series = bsxfun(@times, series, 1./(std_series+eps) );
% % end
% 
% % Normalize to zero mean across subjects
% series = bsxfun(@minus, series, mean(series, 2) );
% %tic, clustered = direcClus_fix_bessel_bsxfun(series, num_clusters, size(series, 2) - 1, num_tries, lambda, 0, 0, 1e-4, 1, max_iter, 1); toc
% 
% %% Actual function
% Kmax = 22
% Kmin = 1
% B = 10 %number of different random permutations for each k
% d = 0 % the dimensionality. If enter zero, the actual dimensionality of the data
% nruns = 10 %number of random initializations for each run B of a given k
% epsilon = 1e-4

%out = CBIG_runresamplingK_single(series,Kmax,Kmin,B,d,nruns,epsilon)

% save('resamplingk.mat')
% save('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/resampling_out.mat')

%% Estimate the stability
load('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/resampling_out.mat')

CBIG_determineK_single(series,out)

save('stability.mat')
save('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/stability.mat')

