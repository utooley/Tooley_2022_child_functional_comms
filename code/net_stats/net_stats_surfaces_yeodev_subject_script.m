function net_stats_surfaces_yeodev_subject_script(sub)
%% Paths
outdir='/cbica/projects/spatial_topography/data/imageData/net_stats/site16_subjectwise'
nets_dir='/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/'
%data_dir='/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek'
subjlist_dir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks'
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/WSBM_v1.2'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/system_matrix_tools/'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/BCT/'))

%% Load input data
file='yeo7_n670_2runsonly_1000tries.mat'
clustered=load(fullfile(nets_dir, file))
%read in yeo communities, should be 81k total
labels=[clustered.lh_labels
    clustered.rh_labels]

listdir=strcat('/cbica/projects/spatial_topography/data/subjLists/release2/site16')
fc_matrices_dir=strcat('/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/fsaverage6_full_corrmats/')
%load subject list
%% Load in surfaces and compute Yeo-dev connectivity matrices
%clear variables that might have overlapped with previous subject list
%load in the FC matrix for each subject

%The more principled way to do this is to create the full
%surf-to-surf correlation matrix, and then average w/in Yeodev ROIs
%with the code below

%for each subject, write their left and right hemispheres into two text
%files using bash script 
surflist1=fullfile(listdir, 'yeo_networks','subjects', strcat(sub, '_fullpaths_surf2surf_lh.txt'))
surflist2=fullfile(listdir, 'yeo_networks','subjects', strcat(sub, '_fullpaths_surf2surf_rh.txt'))

%then use this to compute the full surface 80k x 80k correlation for that subject

%CBIG_ComputeFullSurfaceCorrelation('testoutputfullsurf.mat' , varargin_text1, varargin_text2, pval)
outfile=fullfile(fc_matrices_dir, strcat(sub, '_fsaverage_corr_mat.mat'))
%load the file if it exists already
if exist(outfile)
       load(outfile)
   else
    fprintf('Calculating corr matrix');
    CBIG_ComputeFullSurfaceCorrelation(outfile , surflist1, surflist2, 0)
    load(outfile)
end
    %take out the medial wall labels!
    [medwall]=ismember(labels, 0);
    a=corr_mat((medwall==0),:);
    corr_mat=a(:,(medwall==0));
    labels=labels(labels~=0);
    
    avgweight=mean(corr_mat(corr_mat~=0));
    %then average connectivity within the Yeodev communities
    [S, W, B] = segregation(corr_mat,labels);
    system_segreg_yeodev=S;
    mean_within_sys_yeodev=W;
    mean_between_sys_yeodev=B;
    %Connectivity between each of the 7 blocks
    Ci=labels;
    nCi = unique(Ci);
    M=corr_mat;

    for i = 1:length(nCi) % loop through communities
        for j = 1:length(nCi)
           Wi = Ci == nCi(i); % find index for within communitiy edges
           Bi = Ci == nCi(j); % find index for between communitiy edges to specific community

           Wv_temp = M(Wi,Wi); % extract within communitiy edges
           Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community

           %calculate standard deviation of edges within blocks here
           %deviation_edge_weights_yeo(n,1)=deviation_edge_weights_yeo(n,1)+std(Wv_temp(Wv_temp~=0))^2; %Gu paper calculation, finished below
           %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
           Bv = [Bv_temp(:)'];
           system_connectivity(i,j)=mean(Bv(Bv~=0));
           %system_between(i,1)=mean(Bv);
           %if i==j
           %else
        end
    end

    %transpose these (or something) so they can be saved out on a subject basis.
    system_connectivity;
    system_connectivity_all=system_connectivity;  %put this all into a matrix for everyone so that we can see the average
    yeodev_conn_vector = reshape(system_connectivity',[],1)';
    yeodev_connectivity=yeodev_conn_vector;
    %other metrics work? Take too much memory, so skipping them for now.
    %[Ppos Pneg]=participation_coef_sign(corr_mat, labels);
    %each row is a subject, with part coef for 359 nodes
    %sub_partcoef_pos_yeodev=mean(Ppos);
    %sub_partcoef_neg_yeodev=mean(Pneg);
    %estimate modularity with this a-priori partition?
     %% CHECK THAT THESE TWO Q's ARE THE SAME %%
    %[M Q]=modul_only(corr_mat, [], double(labels), 'negative_asym');
    %QFModul(double(yeo_nodes), subfcmat) %this gives the same answer as modul_only above when using on a pos-only matrix, but 
    %we want the negative asymmteric weighting, so ignore for now.
    % [M Q]=community_louvain(subfcmat, [], double(yeo_nodes), 'negative_asym'); %this optimizes modularity, which is not what we ant.
    %modul_yeodev(n,1)=Q;


%save outfile
%outfile=dataset(sub, avgweight, system_segreg_yeodev, mean_within_sys_yeodev, mean_between_sys_yeodev, yeodev_connectivity)
%export(outfile,'File',fullfile(outdir,strcat(sub,'_site16_fsaverage6_yeodev_network_stats.csv')),'Delimiter',',')

%try to save outfile with correct header
header={'ID', 'avgweight', 'system_segreg_yeodev', 'mean_within_sys_yeodev', 'mean_between_sys_yeodev','sys1to1','sys1to2','sys1to3','sys1to4','sys1to5','sys1to6','sys1to7','sys2to1','sys2to2','sys2to3','sys2to4','sys2to5','sys2to6','sys2to7','sys3to1','sys3to2','sys3to3','sys3to4','sys3to5','sys3to6','sys3to7','sys4to1','sys4to2','sys4to3','sys4to4','sys4to5','sys4to6','sys4to7','sys5to1','sys5to2','sys5to3','sys5to4','sys5to5','sys5to6','sys5to7','sys6to1','sys6to2','sys6to3','sys6to4','sys6to5','sys6to6','sys6to7','sys7to1','sys7to2','sys7to3','sys7to4','sys7to5','sys7to6','sys7to7'}

outfile=table(sub, avgweight, system_segreg_yeodev, mean_within_sys_yeodev, mean_between_sys_yeodev, yeodev_connectivity)
save(fullfile(outdir, strcat(sub,'_site16_fsaverage6_yeodev_network_stats.csv')), 'outfile')
writetable(outfile,fullfile(outdir,strcat(sub,'_site16_fsaverage6_yeodev_network_stats.csv')))

outfile2=splitvars(outfile, 'yeodev_connectivity')
outfile2.Properties.VariableNames=header

save(fullfile(outdir, strcat(sub,'_site16_fsaverage6_yeodev_network_stats.csv')), 'outfile2')
writetable(outfile2,fullfile(outdir,strcat(sub,'_site16_fsaverage6_yeodev_network_stats.csv')))

%outfile=table(char(unique_subjlist.id0), avgweight, modul, avgclustco_both, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg, mean_within_sys, mean_between_sys, system_conn)

end

%% UNUSED

% yeo_dev = fullfile(nets_dir,'lh.yeodev.fsaverage6.annot')
% subjlist1 = fullfile(subjlist_dir, 'n3_site16_tworunsonly_volume_filepaths.txt')
% motionlist = fullfile(subjlist_dir, 'n3_site16_tworunsonly_motion_filepaths.txt')
% 
% % CBIG_ComputeROIs2ROIsCorrelationMatrix(output_file, subj_text_list1, subj_text_list2, discard_frames_list, ROIs1, ROIs2, regression_mask1, regression_mask2, all_comb_bool, avg_sub_bool)
% CBIG_ComputeROIs2ROIsCorrelationMatrix(fullfile(outdir,'testoutputfile.mat'), subjlist1, subjlist1 , motionlist, yeo_dev,yeo_dev,'NONE','NONE', 1, 0)
% %this corr mat is file of 7 x 7 x n subjects correlation matrices, just
% %averaging within each of the 7 Yeodev ROIs for the left hemisphere--this
% %is problematic because I can only do left and rh separately.
% 
% %I m
% 
% %load the file back in
% 
% %turn it into an n x 7 matrix that we can write to a csv
% for n=1:size(subjlist,1)
%     sub=subjlist{n,:}
%     yeo_conn_vector = reshape(system_connectivity(:,:,n)',[],1)';
%     yeo_conn_all(n,:)=yeo_conn_vector
% end
% 
% write(subjlist, yeo_conn_all)

