%% Writing WSBM partition to a Freesurfer annotation file (same process for uncertainty)
%load WSBM consensus iterative labels
outdir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/'
wsbm_consensus=(fullfile(outdir, 'n670_training_sample_consensus_partitions_yeorelabeled.mat'))
load(wsbm_consensus)
%for WSBM partition labels
values=consensus_iter_mode_yeorelabeled;

%for variability in community assignment labels
% freq=(fullfile(outdir, "consensus_iter_mode.mat"))
% load(freq)
% values=abs(freq-670)';

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
copy.orig_tab='wsbm_consensus' %name the new colortable
copy.struct_names=cell(numel(unique(newlabel)),1); %%make enough structure names for num of labels
ranked=sort(unique(newlabel))%rank the values (community assignments) from low to high
copy.struct_names{1}='Unknown'  %make the first struct_name the empty unknown
copy.table(1,:)=[1 1 1 0 0];  %make the first colortable entry the empty unknown
for i=1:length(unique(newlabel))
    %make the struct_names the label values, in order of magnitude
    %make sure this works
    copy.struct_names{i}=[char(num2str(ranked(i)))];
end
%make an empty colortable of 0's
copy.table=zeros(numel(unique(newlabel)),5)

% COLORMAP--this doesn't matter if visualizing in R afterwards
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
  1     7Networks_1 120  18 134   0
  2     7Networks_2  70 130 180   0
  3     7Networks_3   0 118  14   0
  4     7Networks_4 196  58 250   0
  5     7Networks_5 220 248 164   0
  6     7Networks_6 230 148  34   0
  7     7Networks_7 205  62  78   0
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
copy.orig_tab='wsbm_consensus' %name the new colortable
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