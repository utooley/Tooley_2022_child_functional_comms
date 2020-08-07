
cd /cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks
load averaged_FC_mat_n670.mat %my 400-node network
load /cbica/projects/spatial_topography/data/imageData/wsbm/% wsbm labels

% generate minimum spanning tree
mst = graphminspantree(sparse(max(meanMatrix(:)) - meanMatrix),'method','kruskal');

% node properties
rngSz = [5,50];                     % range of node sizes
nodeSz = fcn_sz(sum(meanMatrix,2),rngSz);    % node sizes

% visualization parameters
mapSz = 1000;    % size of map
sigma = 10;      % radius of smoothed colors

% generate thresholded network and union with mst
thr = 0.02;
A = double(threshold_proportional(meanMatrix,thr) | mst | mst');

% make figure
[f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(A,lab,cmap,nodeSz,mapSz,sigma);

[X,Y,Z] = fcn_adjacency_plot(aij,xy)