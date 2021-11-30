%PATHS
datadir=fullfile('/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/')
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/parcellation'
z_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zNetworks'
noz_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400Networks'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample'
addpath(genpath('/data/picsl/mackey_group/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/WSBM_v1.2'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/BCT'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/NCT'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/'))

%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess_long.csv')) %the first two runs here are those we use for WSBM and clustering

% Z-score FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.id(n))
    run1=char(subjlist.var1(n)) %%% FOR RUN 1 %%
    file=fullfile(datadir,sub, run1, strcat('fcon/schaefer400x7/',sub,'_',run1,'_schaefer400x7_network.txt'));
    try
    subfcmat=load(file);
    %make into adjacency matrix and save out
    size_vec=tril(ones(400,400),-1);
    adj_mat=size_vec;
    adj_mat(adj_mat==1)=subfcmat;
    subfcmat=adj_mat+adj_mat';
    outfile=fullfile(noz_outdir,strcat(sub,'_',run1,'_schaefer400x7_network.txt'));
    csvwrite(outfile, subfcmat);
    for x=1:400
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc1=[];
    for i=1:400
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc1(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(z_outdir, strcat(sub,'_', run1,'_Schaefer400x7_znetwork.txt'));
    csvwrite(outfile, zfc1);
    run2=char(subjlist.var2(n))  %%% FOR RUN 2 %%
    file=fullfile(datadir,sub, run2, strcat('fcon/schaefer400x7/',sub,'_',run2,'_schaefer400x7_network.txt'));
    subfcmat=load(file);
    size_vec=tril(ones(400,400),-1); %make into adjacency matrix and save out
    adj_mat=size_vec;
    adj_mat(adj_mat==1)=subfcmat;
    subfcmat=adj_mat+adj_mat';
    outfile=fullfile(noz_outdir,strcat(sub,'_',run2,'_schaefer400x7_network.txt'));
    csvwrite(outfile, subfcmat);
    for x=1:400
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc2=[];
    for i=1:400
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc2(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(z_outdir, strcat(sub,'_', run2,'_Schaefer400x7_znetwork.txt'));
    csvwrite(outfile, zfc2);
    %% AVERAGE THE TWO RUNS %%
    aggregmat=(zfc2+zfc1)/2;
    outfile=fullfile(z_avg_outdir, strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
    csvwrite(outfile, aggregmat);
    catch
    fprintf('Cant read sub %s run %s, skipped. \n', sub);
  end
end

%% Run WSBM with k=7

k=7
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'

for n=1:height(subjlist)
    sub=char(subjlist.id(n)) %look at this
    outfile=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'))
    if (exist(outfile)==2) %if it's already written don't do it again
        fprintf('Sub %s already exists. \n', sub);
    else
    file=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
    try
    subfcmat=load(file);
    [Labels Model]=wsbm(subfcmat, k,'E_Distr','None', 'W_Distr', 'Normal', 'numTrials', 50, 'verbosity', 1, 'alpha', 0, 'parallel', 0);
    %logHw - Additive Log-likelihood constant for W_Distr
    sub_log_lik_for_weights(n,1)=Model.Data.logHw;
    %logHe - Additive Log-likelihood constant for E_Distr
    sub_log_lik_for_edges(n,1)=Model.Data.logHe;
    %LogEvidence - Marginal Log-likelihood (aka Log-Evidence), a model selection criterion 
    sub_log_evidence(n,1)=Model.Para.LogEvidence;
    %a row for each subject's labels for the k=7 partition
    sub_labels(n,:)=Labels;
    sub_unique_labels(n,1)=numel(unique(Labels));
    models(n, :)=Model;
    %save the model
    outfile= fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'))
    save(outfile, 'Model', 'Labels')
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
    end
end

outfile=dataset(subjlist, sub_log_lik_for_weights, sub_log_lik_for_edges, sub_log_evidence, sub_unique_labels, sub_labels)
export(outfile,'File',fullfile(wsbm_dir,'wsbm_k7_n670_site16_50trials.csv'),'Delimiter',',')

%% Relabel community partitions, create a consensus partition
yeo_nodes=dlmread('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')

%read in each subject's wsbm partition
for n=1:height(subjlist)
    n
    sub=char(subjlist.id(n));
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'));
    %file=fullfile(wsbm_dir,strcat('thresh_',density,'/',sub,'_wsbm_thresh.mat'))
    try 
        load(file);
        sub_log_evidence(n,1)=Model.Para.LogEvidence;
        sub_unique_labels(n,1)=numel(unique(Labels));
        models(n, :)=Model;
        partition=Labels;%load in the partition from the WSBM for this subject
        part_matrix(:,n)=partition;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end

%save subject-level variables if not already done.
outfile=dataset(subjlist.id, sub_log_evidence, sub_unique_labels, sub_labels)
export(outfile,'File',fullfile(wsbm_dir,'wsbm_k7_n670_site16_thresh_50trials.csv'),'Delimiter',',')

%input=400 by p (nodes by partitions) matrix, each subj is a column, n
%subjects
noyeo_opt_part_matrix=multislice_pair_labeling(part_matrix) %this is not actually necessary, can be skipped

%% Create consensus partition
%use consensus iterative procedure (no relabeling is necessary, can skip above 2 lines)
gamma=2 %when keeping the thresholding of the consensus matrix, this is the right gamma
[consensus_mat_iter_noyeo Q2 X_new3 qpc cooccurence_matrix]=consensus_iterative(part_matrix', gamma);
[consensus_iter_mode freq ties]=mode(consensus_mat_iter_noyeo) %extract assignment variance across optimizations

%save these outfiles
outfile=('/cbica/projects/spatial_topography//data/imageData/wsbm/site16_training_sample/brains/freq.mat')
save(outfile, 'freq')

%% Look at versatility of the co-occurence matrix
% load part_matrix
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'

gamma=2 %when keeping the thresholding of the consensus matrix, this is the right gamma
[consensus_mat_iter_noyeo Q2 X_new3 qpc cooccurence_matrix]=consensus_iterative(part_matrix', gamma);
[consensus_iter_mode freq ties]=mode(consensus_mat_iter_noyeo) %extract assignment variance across optimizations

% for x=1:400
%     cooccurence_matrix(x,x)= NaN;
% end
cooccurence_matrix_prop = cooccurence_matrix/670;
%This is the association matrix, so now we can look at versatility.

%adapted from versatility code at https://github.com/mwshinn/versatility/tree/master/matlab
Cs = sin(pi*cooccurence_matrix_prop);
V = sum(Cs, 1,'omitnan')./sum(cooccurence_matrix_prop, 1,'omitnan'); % CM/Cs are symmetric so axis doesn't matter
V(V<1e-10) = 0; % Stupid floats

outfile=('/cbica/projects/spatial_topography//data/imageData/wsbm/site16_training_sample/brains/versatility_n670.mat')
save(outfile, 'V')

%% Relabel and save the consensus partitions so they are visually comparable to Yeo7
outdir='/cbica/projects/spatial_topography//data/imageData/wsbm/site16_training_sample/brains/'
%relabel the consensus partitions so that they are visually comparable to Yeo
consensus_iter_mode_wsbmrelabeled=multislice_pair_labeling([consensus_iter_mode_yeorelabeled consensus_iter_mode']);
consensus_iter_mode_wsbmrelabeled=consensus_iter_mode_wsbmrelabeled(:,2);

outfile=(fullfile(outdir, 'n670_training_sample_consensus_partitions_yeorelabeled.mat'))
save(outfile, 'consensus_iter_mode_yeorelabeled')

%% SAVE OUTFILES
outfile=(fullfile(outdir, 'n670_training_sample_consensus_mat_and_nodal_variance.mat'))
save(outfile, 'noyeo_opt_part_matrix', 'part_matrix', 'freq', 'ties', 'consensus_iter_mode')