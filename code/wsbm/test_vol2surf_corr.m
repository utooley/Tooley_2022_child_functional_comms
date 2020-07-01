

%% Use silhouette measure on WSBM partition
%paths
zavg_dir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'

%what is series?
load('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/resamplingk.mat')
size(series)
%it's 1175 x 74k (num vertices - medial wall)
%Load training sample site 16 WSBM average corr matrix
load(fullfile(zavg_dir,'averaged_FC_mat_n670.mat'))
%load assignments
load('/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/n670_training_sample_consensus_partitions_yeorelabeled.mat')
consensus_nodes=consensus_iter_mode_yeorelabeled;

%maybe I should replace the 0's on the diagonal with 1's?
tic; s = silhouette(meanMatrix,consensus_nodes, 'correlation'); toc

%try replacing 0's with 1's, see if different
    for x=1:400
        meanMatrix(x,x)=NaN;
    end
tic; s_2 = silhouette(meanMatrix,consensus_nodes, 'correlation'); toc


if(num_clusters > 1 && ~isempty(strfind(mesh_name, 'fsaverage')) && no_silhouette==0)
    tic; s = silhouette(series, cidx(non_zero_corr_index), 'correlation'); toc
    if(sum(non_zero_corr_index) < length(non_zero_corr_index)) %there are vertices with no correlation
        new_s = ones(length(non_zero_corr_index), 1);
        new_s(non_zero_corr_index) = s;
        s = new_s; 
    end
else
    s = zeros(length(l), 1);
end
lh_s = ones(lh_num_verts, 1);
lh_s(l1) = s(1:length(l1));

rh_s = ones(rh_num_verts, 1);
rh_s(l2) = s(length(l1)+1:end);


%compute seed-based correlation?
roi_file='~/Documents/tools/parcellations/Yeo_ThalamusParcellation/test.nii'
output_file='~/test.nii'


CBIG_ComputeCorrelationProfileVol2Surf(roi_file, output_file, 1, surf_input_files_txt, vol_input_files_txt, regress_bool)
