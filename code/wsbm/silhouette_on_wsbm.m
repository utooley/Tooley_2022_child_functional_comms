%% Use silhouette measure on WSBM partition
%paths
zavg_dir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'

load(fullfile(zavg_dir,'averaged_FC_mat_n670.mat'))
%load assignments
load('/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/n670_training_sample_consensus_partitions_yeorelabeled.mat')
consensus_nodes=consensus_iter_mode_yeorelabeled;

%maybe I should replace the 0's on the diagonal with 1's?
tic; s = silhouette(meanMatrix,consensus_nodes, 'correlation'); toc

outfile=('/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/silhouettes_wsbm.mat')
save(outfile, 's')

