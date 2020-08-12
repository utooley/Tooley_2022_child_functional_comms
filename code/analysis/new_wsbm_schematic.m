
cd /cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks
load averaged_FC_mat_n670.mat %my 400-node network
load /cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/n670_training_sample_consensus_partitions_yeorelabeled.mat% wsbm labels
load ~/Documents/tools/MATLAB/Rick_betzel_layoutCode/layoutData.mat

% generate minimum spanning tree
mst = graphminspantree(sparse(max(meanMatrix(:)) - meanMatrix),'method','kruskal');

% node properties
rngSz = [5,50];                     % range of node sizes
nodeSz = fcn_sz(sum(meanMatrix,2),rngSz);    % node sizes

% visualization parameters
mapSz = 1000;    % size of map
sigma = 20;      % radius of smoothed colors

% generate thresholded network and union with mst
thr = 0.025;
A = double(threshold_proportional(meanMatrix,thr) | mst | mst');

% make figure
[f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(A,consensus_iter_mode_yeorelabeled,cmap,nodeSz,mapSz,sigma);

%customize colormap
colors=["000000", "7B287E", "5CA1C8", "0A9045","C33AF8", "dcf8a4", "EF9C23", "E34A53"]
cmap(1,:)=hex2rgb('000000')
cmap(2,:)=hex2rgb('7B287E')
cmap(3,:)=hex2rgb('5CA1C8');cmap(4,:)=hex2rgb('0A9045');cmap(5,:)=hex2rgb('C33AF8');cmap(6,:)=hex2rgb('dcf8a4');cmap(7,:)=hex2rgb('EF9C23');
cmap(8,:)=hex2rgb('E34A53');
cmap=cmap(2:8,:)%take out the gray at first

%% Cycle through different assortative and non-assortative structures
outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
load(fullfile(outdir,'cmap.mat'))
cd /cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks
load averaged_FC_mat_n670.mat %my 400-node network

colormap(cm)
imagesc(meanMatrix) %raw

labels=[repmat(1,[1,100]) repmat(2,[1,100]) repmat(3,[1,100]) repmat(4,[1,100])]
%assortative-100 in each of 4 communities
for i=1:4;
    idx = labels == i;  
    meanMatrix(idx, idx)= abs(meanMatrix(idx,idx)* 2)
end
imagesc(meanMatrix)
mst = graphminspantree(sparse(max(meanMatrix(:)) - meanMatrix),'method','kruskal');
thr = 0.01;
A = double(threshold_proportional(meanMatrix,thr) | mst | mst');
[f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(A,labels,cmap,nodeSz,mapSz,sigma);

%disassortative
load averaged_FC_mat_n670.mat %my 400-node network
for i=1:4;
    idx = labels == i;  
    idx2= labels==i+1;
    meanMatrix(idx, idx2)= abs(meanMatrix(idx,idx2)* 2)
    meanMatrix(idx2, idx)= abs(meanMatrix(idx2,idx)* 2)
    %meanMatrix(idx, idx)= -(meanMatrix(idx,idx)/2)
    %meanMatrix(idx, ~idx2)= -(meanMatrix(idx,~idx2)/2)
end 
imagesc(meanMatrix)

%core-periphery
load averaged_FC_mat_n670.mat %my 400-node network
idx = labels == 2;  
meanMatrix(idx, idx)= abs(meanMatrix(idx,idx)* 2)
meanMatrix(idx, ~idx)= abs(meanMatrix(idx,~idx)* 2)
meanMatrix(~idx, idx)= abs(meanMatrix(~idx,idx)* 2)
imagesc(meanMatrix)
% mst = graphminspantree(sparse(max(meanMatrix(:)) - meanMatrix),'method','kruskal');
% thr = 0.015;
% A = double(threshold_proportional(meanMatrix,thr) | mst | mst');
% [f,ax,nodeHandles,edgeHandles] = fcn_community_boundaries(A,labels,cmap,nodeSz,mapSz,sigma);

