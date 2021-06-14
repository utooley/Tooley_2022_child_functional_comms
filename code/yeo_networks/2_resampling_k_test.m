%% Use suggested code by Ruby Kong, use ConsistencySurf.m instead of resamplingk.m and determinek.m
addpath('/cbica/projects/spatial_topography/code/yeo_networks')
%input profiles data.
profile1='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/lh.yeo7_n670_2runsonly_1000tries.avg_profiles007.nii.gz'
profile2='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/rh.yeo7_n670_2runsonly_1000tries.avg_profiles007.nii.gz'
mesh_name='fsaverage6'
mask='cortex'
num_smooth=0;
num_tries=100 %1000
rand_num=100
dim=1
normalize=1
output_dir='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/search_over_k/'

for k = 2:17
CBIG_VonmisesSeriesConsistencySurf(mesh_name, mask, k, output_dir, profile1, profile2, num_smooth, num_tries, rand_num, dim, normalize)
end

%Aggregate results into a csv
clear consistency_true consistency_rand stability_rand stability_true
for k = 2:17
    if k >= 10
        load(fullfile(output_dir,  strcat('Cluster0',num2str(k),'.s00.tries100.rand100.znorm1.dim1..mat')))
    else
         load(fullfile(output_dir,  strcat('Cluster00',num2str(k),'.s00.tries100.rand100.znorm1.dim1..mat')))
    end
    consistency_true(:,k)=con_struct.orig_overlap'
    consistency_rand(:,k)=con_struct.rand_overlap'
    try
    stability_true(:,k)=con_struct.stab
    stability_rand(:,k)=con_struct.rand_stab
    catch
    end
end
outfile=dataset(consistency_true, consistency_rand,stability_true,stability_rand)
export(outfile,'File',strcat(output_dir,'/k2_to_17_tries100_rand100.znorm1.dim1.csv'),'Delimiter',',')