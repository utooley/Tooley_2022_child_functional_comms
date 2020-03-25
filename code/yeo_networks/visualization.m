%Checking Yeo code outputs
cluster_dir='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/'
file_name='yeo7_n670_2runsonly_1000tries.mat'

clustered = load(fullfile(cluster_dir,file_name));

ref_file = ['/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/1000subjects_clusters007_ref.mat'];
ref = load(ref_file); 

my_colors = [190   190  190
160    32  240
0     0  255
127   255  212
255   246  143
255     0  255
255   165    0
255     0    0]
    
    

%PLOT 7 NETWORKS-can change the colors to match WSBM colors?
CBIG_DrawSurfaceMaps(clustered.lh_labels, clustered.rh_labels, 'fsaverage6','inflated',0, 7, ref.colors)

CBIG_DrawSurfaceMaps(clustered.lh_labels, clustered.rh_labels, 'fsaverage5','inflated',0, 17, ref.colors)

%read in WSBM data and relabel Schaefer400 annot file with WSBM assignments
%can use CBIG_read_annotation.m
%then use CBIG_VisualizeSurfaceAnnotationInFreeview to view.

%% Save developmental communities to Freesurfer annotation
%CBIG_SaveParcellationToFreesurferAnnotation(fullfile(cluster_dir,file_name), 'lh.yeodev.fsaverage.annot', 'rh.yeodev.fsaverage.annot')
%this uses the colors in the wrong order.

CBIG_WriteParcellationToAnnotation(clustered.lh_labels,fullfile(cluster_dir,'lh.yeodev.fsaverage6.annot'), ref.colors)
CBIG_WriteParcellationToAnnotation(clustered.rh_labels,fullfile(cluster_dir,'rh.yeodev.fsaverage6.annot'), ref.colors)
%this is what I want! With my colors and in fsaverage6 space.

%% Confidence maps
%Plot Yeo7 confidence maps
yeo_ref_dir='/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering/1000subjects_reference/'

CBIG_DrawSurfaceMaps(ref.lh_s, ref.rh_s, 'fsaverage5','inflated',0, 0.4)
%make a colormap by hand in the GUI to save with the annotation

%and save them as an annotation
CBIG_WriteParcellationToAnnotation(ref.lh_s,fullfile(yeo_ref_dir,'lh.silhouette.fsaverage5.annot'), mycmap)
CBIG_WriteParcellationToAnnotation(ref.rh_s,fullfile(yeo_ref_dir,'rh.silhouette.fsaverage5.annot'), ref.colors)
CBIG_WriteParcellationToAnnotation(clustered.rh_labels,fullfile(cluster_dir,'rh.yeodev.fsaverage6.annot'), ref.colors)

abs_path_to_lh_annot_file = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage5', 'label', 'lh.aparc.annot');
abs_path_to_rh_annot_file = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage5', 'label', 'rh.aparc.annot');

CBIG_DrawSurfaceDataAsAnnotation(ref.lh_s, ref.rh_s, abs_path_to_lh_annot_file, abs_path_to_rh_annot_file, 0, 'fsaverage5', yeo_ref_dir, 'test')

[test colortable]=CBIG_ConvertSingleHemiSurfaceDataToDiscretizedAnnotation(ref.lh_s, 15, 'hsv',0,1)  
read_annotation('/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/fsaverage5/lh.Schaefer2018_400Parcels_7Networks_order.annot')


%Silhouettes/confidence maps-can change the colors?
load('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/yeo7_n670_2runsonly_1000tries.mat')
load('/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n131_5tries/silhou_cmap.mat')

CBIG_DrawSurfaceMaps(clustered.lh_s, clustered.rh_s, 'fsaverage6','pial',0, 0.4, color)


