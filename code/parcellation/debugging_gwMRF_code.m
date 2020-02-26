%left hemisphere
l=load('/Volumes/home/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/gwMRF_schaefer400_site16_n670/mult_mat/lh_mult_matrix.mat')
time_data/lh_time_matrix.mat
scans=l.scans
files_used_lh=l.files_used_lh

%who has NaNs in their code?
none=isnan(l.lh_time_mat);
whichsubjs=sum(none,1);
index=whichsubjs==0
new_lh_time_mat=l.lh_time_mat(:,index);
lh_time_mat=new_lh_time_mat;
save('~/Downloads/time_data/lh_time_matrix2.mat','lh_time_mat','scans','files_used_lh','-v7.3');

%right hemisphere
load('~/Downloads/time_data/rh_time_matrix.mat')
%who has NaNs in their code?
none=isnan(rh_time_mat);
whichsubjs=sum(none,1);
index=whichsubjs==0
new_rh_time_mat=rh_time_mat(:,index);
rh_time_mat=new_rh_time_mat;
save('~/Downloads/time_data/rh_time_matrix2.mat','rh_time_mat','scans','files_used_rh','-v7.3');
%figure out which subjects it is
bad=find(whichsubjs>=1);
subjnumbers=bad/740;
%figure out which vertices it is.
whichvertices=sum(none,2);

%searching for the right k0, put a breakpoint in CBIG_Cldn at k0 and search for it on the cluster 
for k0=3000:51200
k(k0)=k0
f(ko)=log(besseli(d/2-1,k0))
end