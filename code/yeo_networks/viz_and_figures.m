%% Make an average connectivity matrix for Schaefer 400
zavg_net_dir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/parcellation'
outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess.csv')) %the first two runs here are those that were input into gwMRF

% Only include the averaged weighted connectivity matrix for each subject
unique_subjlist= subjlist.id
stackedMatrix = zeros(400, 400, size(unique_subjlist,1));
for n=1:size(unique_subjlist,1)
    sub=char(unique_subjlist(n))%look at this
    file=fullfile(zavg_net_dir,strcat(num2str(sub),'_avg_Schaefer400x7_znetwork.txt'))
    %file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400MNI_znetwork.txt'))
    subfcmat = load(file);
    for x=1:400
        subfcmat(x,x)=0;
    end
    stackedMatrix(:,:, n)=subfcmat;
end
meanMatrix = mean(stackedMatrix,3); %doc mean for more info.
%make the diagonal 0's
for x=1:400
    meanMatrix(x,x)=0;
end
%export mean matrix
csvwrite(fullfile(outdir,strcat('averaged_FC_mat_n', num2str(size(unique_subjlist,1)),'.csv')), meanMatrix)
save(fullfile(outdir, strcat('averaged_FC_mat_n', num2str(size(unique_subjlist,1)),'.mat')), 'meanMatrix')
cm = colormap(gca)
save(fullfile(outdir,'cmap.mat'),'cm')

%% Make the ROI x vertex-wise corr profile matrix for clustering

load('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/resamplingk.mat')
colormap(cm)
imagesc(series)

%Checking Yeo code outputs
cluster_dir='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/'
file_name='yeo7_n670_2runsonly_1000tries.mat'

clustered = load(fullfile(cluster_dir,file_name));

ref_file = ['/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/1000subjects_clusters007_ref.mat'];
ref = load(ref_file);

%% Save developmental communities to Freesurfer annotation
%CBIG_SaveParcellationToFreesurferAnnotation(fullfile(cluster_dir,file_name), 'lh.yeodev.fsaverage.annot', 'rh.yeodev.fsaverage.annot')
%this uses the colors in the wrong order.

CBIG_WriteParcellationToAnnotation(clustered.lh_labels,fullfile(cluster_dir,'lh.yeodev.fsaverage6.annot'), ref.colors)
CBIG_WriteParcellationToAnnotation(clustered.rh_labels,fullfile(cluster_dir,'rh.yeodev.fsaverage6.annot'), ref.colors)
%this is what I want! With my colors and in fsaverage6 space.

%and save them as an annotation
CBIG_WriteParcellationToAnnotation(ref.lh_s,fullfile(yeo_ref_dir,'lh.silhouette.fsaverage5.annot'), mycmap)
CBIG_WriteParcellationToAnnotation(ref.rh_s,fullfile(yeo_ref_dir,'rh.silhouette.fsaverage5.annot'), ref.colors)
CBIG_WriteParcellationToAnnotation(clustered.rh_labels,fullfile(cluster_dir,'rh.yeodev.fsaverage6.annot'), ref.colors)

