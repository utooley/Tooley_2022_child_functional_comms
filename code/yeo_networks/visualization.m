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

%% ZRand of partitions
%Yeo7 vs. developmental Yeo
%read in annotation of Schaefer in fsaverage6
[vertices yeo7_lh colortable]=read_annotation('/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/fsaverage6/label/lh.Yeo2011_7Networks_N1000.annot')
[vertices yeo7_rh colortable]=read_annotation('/Users/utooley/Documents/tools/parcellations/Yeo_from_freesurfer/fsaverage6/label/rh.Yeo2011_7Networks_N1000.annot')
yeo7=vertcat(yeo7_lh, yeo7_rh)
yeo_dev=vertcat(clustered.lh_labels,clustered.rh_labels)

[zr sr sar vi]=zrand(yeo7, yeo_dev)
%take out the 0's
yeo_7=yeo7(yeo7 ~= 65793);
yeo_dev=yeo_dev(yeo_dev~=0);
%Compare Yeo7 to Yeo developmental in fsaverage6

[zr sr sar vi]=zrand(yeo_7, yeo_dev)

%load WSBM consensus iterative vs. yeo assignments
outdir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/'
wsbm_consensus=(fullfile(outdir, 'n670_training_sample_consensus_partitions_yeorelabeled.mat'))
load(wsbm_consensus)
%Yeo in Schaefer nodes
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')

[zr sr sar vi]=zrand(consensus_iter_mode_yeorelabeled, yeo_nodes)

%% Relabel Schaefer annotation file with WSBM labels
%load WSBM consensus iterative labels
outdir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/'
wsbm_consensus=(fullfile(outdir, 'n670_training_sample_consensus_partitions_yeorelabeled.mat'))
load(wsbm_consensus)
values=consensus_iter_mode_yeorelabeled;
%split into lh and rh, and add a 0 at the beginning for the medial wall
%label
values_lh=vertcat(0,values(1:200));
values_rh=vertcat(0,values(201:400));

%% LH ANNOTATION FILE
schaefer_parcellation_dir='/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/'
[num_vertices label colortable]=read_annotation(fullfile(schaefer_parcellation_dir,'fsaverage6/label/lh.Schaefer2018_400Parcels_7Networks_order.annot'));

%c%create a numbered list with region names
mapping = table((1:1:length(colortable.struct_names))', colortable.struct_names)
%make empty vector for new labels
newlabel=label;
%for each region in the lh annotation
for n=1:numel(unique(colortable.struct_names));
    name=char(colortable.struct_names(n));
    %the row of mapping that matches the structure name
    ind=strmatch(name, mapping{:,2}, 'exact'); %match exactly, since otherwise there is overlap in names!
    region_num=mapping{ind,1};
    comm_assign=values_lh(region_num); %find the value for that region in the mapping
    %go into the colortable struct-names and find the colortable.table(:,5)
    %that matches the vertex ID.
    vertind=find(label == colortable.table(n,5));
    newlabel(vertind)=comm_assign;
end

%create a new empty struct named copy for the colortable and fill it
copy = struct();
copy.numEntries=numel(unique(newlabel)); %num entries=num new labels, plus one for the unknown label needed at the beginning
copy.orig_tab='wsbm_consensus_community' %name the new colortable
copy.struct_names=cell(numel(unique(newlabel)),1); %%make enough structure names for num of labels
ranked=sort(unique(newlabel))%rank the values (community assignments) from low to high
copy.struct_names{1}='Unknown'  %make the first struct_name the empty unknown
copy.table(1,:)=[1 1 1 0 0];  %make the first colortable entry the empty unknown
for i=1:length(unique(newlabel))
    %make the struct_names the clust co, in order of magnitude
    %make sure this works
    copy.struct_names{i}=[char(num2str(ranked(i)))];
end
%make an empty colortable of 0's
copy.table=zeros(numel(unique(newlabel)),5)

% COLORMAP
clear cmap
n = numel(unique(newlabel));  %// number of colors
cmap(1,:) = [51 0 102];   %// color first row - dark purple
cmap(2,:) = [153 51 255];   %// color mid row - mid purple
cmap(3,:) = [255 255 255];   %// color last row - light purple

[X,Y] = meshgrid([1:3],[1:numel(unique(newlabel))]);  %// mesh of indices,

%interpolate from 1 to 4 to whatever the number of unique labels is!
cmap = interp2(X([1,4,numel(unique(newlabel))],:),Y([1,4,numel(unique(newlabel))],:),cmap,X,Y); %// interpolate colormap
cmap(:,4)=zeros(numel(unique(newlabel)),1)
cmap(:,5)=zeros(numel(unique(newlabel)),1)
copy.table(1:numel(unique(newlabel)), :)=round(cmap)
copy.table(:,5)=copy.table(:,1) + (copy.table(:,2).*(2.^8)) + (copy.table(:,3).*(2.^16))

% RELABEL
%now transfer the RGB values in the fifth col back to the values
%want each value of the colortable (low to high) to be mapped to the value
%in values
for n=1:(numel(unique(newlabel)));
    vertind=find(newlabel==ranked(n));
    %vertind=find(newlabel == uniqlabels(n));
    newlabel(vertind)=copy.table(n,5);
    %the row of mapping that matches the structure name
end

%write out the annotation file
filename=fullfile(outdir, 'lh.wsbm.consensus.fsaverage6.annot')
write_annotation(filename, num_vertices, newlabel, copy)

%% RH ANNOTATION FILE
schaefer_parcellation_dir='/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/FreeSurfer5.3/'
[num_vertices label colortable]=read_annotation(fullfile(schaefer_parcellation_dir,'fsaverage6/label/rh.Schaefer2018_400Parcels_7Networks_order.annot'));

%create a numbered list with region names
mapping = table((1:1:length(colortable.struct_names))', colortable.struct_names)
%make empty vector for new labels
newlabel=label;
%for each region in the lh annotation
for n=1:numel(unique(colortable.struct_names));
    name=char(colortable.struct_names(n));
    %the row of mapping that matches the structure name
    ind=strmatch(name, mapping{:,2}, 'exact'); %match exactly, since otherwise there is overlap in names!
    region_num=mapping{ind,1};
    comm_assign=values_rh(region_num); %find the value for that region in the mapping
    %go into the colortable struct-names and find the colortable.table(:,5)
    %that matches the vertex ID.
    vertind=find(label == colortable.table(n,5));
    newlabel(vertind)=comm_assign;
end

%create a new empty struct named copy for the colortable and fill it
copy = struct();
copy.numEntries=numel(unique(newlabel)); %num entries=num new labels, plus one for the unknown label needed at the beginning
copy.orig_tab='wsbm_consensus_community' %name the new colortable
copy.struct_names=cell(numel(unique(newlabel)),1); %%make enough structure names for num of labels
ranked=sort(unique(newlabel))%rank the values (community assignments) from low to high
copy.struct_names{1}='Unknown'  %make the first struct_name the empty unknown
copy.table(1,:)=[1 1 1 0 0];  %make the first colortable entry the empty unknown
for i=1:length(unique(newlabel))
    %make the struct_names the clust co, in order of magnitude
    %make sure this works
    copy.struct_names{i}=[char(num2str(ranked(i)))];
end
%make an empty colortable of 0's
copy.table=zeros(numel(unique(newlabel)),5)

% COLORMAP
clear cmap
n = numel(unique(newlabel));  %// number of colors
cmap(1,:) = [51 0 102];   %// color first row - dark purple
cmap(2,:) = [153 51 255];   %// color mid row - mid purple
cmap(3,:) = [255 255 255];   %// color last row - light purple

[X,Y] = meshgrid([1:3],[1:numel(unique(newlabel))]);  %// mesh of indices,

%interpolate from 1 to 4 to whatever the number of unique labels is!
cmap = interp2(X([1,4,numel(unique(newlabel))],:),Y([1,4,numel(unique(newlabel))],:),cmap,X,Y); %// interpolate colormap
cmap(:,4)=zeros(numel(unique(newlabel)),1)
cmap(:,5)=zeros(numel(unique(newlabel)),1)
copy.table(1:numel(unique(newlabel)), :)=round(cmap)
copy.table(:,5)=copy.table(:,1) + (copy.table(:,2).*(2.^8)) + (copy.table(:,3).*(2.^16))

% RELABEL
%now transfer the RGB values in the fifth col back to the values
%want each value of the colortable (low to high) to be mapped to the value
%in values
for n=1:(numel(unique(newlabel)));
    vertind=find(newlabel==ranked(n));
    %vertind=find(newlabel == uniqlabels(n));
    newlabel(vertind)=copy.table(n,5);
    %the row of mapping that matches the structure name
end

%write out the annotation file
filename=fullfile(outdir, 'rh.wsbm.consensus.fsaverage6.annot')
write_annotation(filename, num_vertices, newlabel, copy)