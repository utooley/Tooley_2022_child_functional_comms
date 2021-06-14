%% Calculate log-likelihood of different k's of WSBM calculated on replication sample
%Higher k's could just be overfitting the training data, test the consensus partition on the replication sample
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/WSBM_v1.2'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/NCT'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/GenLouvain-2.1'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/'))
addpath(genpath('/cbica/projects/spatial_topography/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/BCT'))
addpath(genpath('/cbica/projects/spatial_topography/code/wsbm/'))

%% First must run WBSM model with different values of k on main sample
datadir=fullfile('/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek/')
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/parcellation'
z_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zNetworks'
noz_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400Networks'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/search_over_k'

%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess.csv')) %the first two runs here are those that were input into gwMRF

k=[2:17] %change this to be the right range
names=strcat('k',string(k))
for k=1:16 %change this to be the right range
    num_comms=k+1
    sub_log_lik_for_weights.(names{k})=zeros(height(subjlist),1);
    sub_log_lik_for_edges.(names{k})=zeros(height(subjlist),1);
    sub_log_evidence.(names{k})=zeros(height(subjlist),1);
    sub_labels.(names{k})=zeros(height(subjlist),400);
    temp_weights=zeros(height(subjlist),1);
    temp_edges=zeros(height(subjlist),1);
    temp_evidence=zeros(height(subjlist),1);
    temp_labels=zeros(height(subjlist),400);
    %parfor n=1:height(subjlist)
    parfor n=1:height(subjlist)
        sub=char(subjlist.id(n)) %look at this
        [temp_weights(n,1),temp_edges(n,1),temp_evidence(n,1),temp_labels(n,:)]=wsbm_function(sub, num_comms,wsbm_dir,z_avg_outdir)
    end
    sub_log_lik_for_weights.(names{k})=temp_weights
    sub_log_lik_for_edges.(names{k})=temp_edges
    sub_log_evidence.(names{k})=temp_evidence
    sub_labels.(names{k})=temp_labels
    
       outfile1=dataset(subjlist.id, sub_log_lik_for_weights.(names{k}), sub_log_lik_for_edges.(names{k}), sub_log_evidence.(names{k}), sub_labels.(names{k}))
       export(outfile1,'File',fullfile(wsbm_dir,strcat('wsbm_search_over_k', num2str(num_comms),'_n670_site16_30trials.csv')),'Delimiter',',')
    end

%% Then create a consensus partition at each k from main sample
yeo_nodes=dlmread('/cbica/projects/cbpd_main_data/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
wsbm_dir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/search_over_k/'
%get the replication subject list
listdir='/cbica/projects//spatial_topography/data/subjLists/release2/site14site20/'
test_subjlist=readtable(fullfile(listdir,'n544_filtered_runs_site14site20_postprocess.csv')) %the first two runs here are those that were input into gwMRF
test_subjlist=test_subjlist.id;
%get main subject list
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/parcellation'
subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess.csv')) %the first two runs here are those that were input into gwMRF

%read in each subject's wsbm partition
for k = [2:17]
for n=1:height(subjlist)
    sub=char(subjlist.id(n));
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub,'_k',num2str(k),'_wsbm.mat'));
    %file=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'));
    try 
        load(file);
        partition=Labels;%load in the partition from the WSBM for this subject
        part_matrix(:,n)=partition;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end

%either use multislice_pair_label or just pair_label
%input=400 by p (nodes by partitions) matrix, each subj is a column, n
%subjects
opt_part_matrix=part_matrix
noyeo_opt_part_matrix=multislice_pair_labeling(opt_part_matrix)

%% Create consensus partition
%consensus iterative partition
gamma=2 %when keeping the thresholding of the consensus matrix, this is the right gamma
[consensus_mat_iter_noyeo Q2 thresholded_association_matrix qpc]=consensus_iterative(noyeo_opt_part_matrix', gamma);
[consensus_iter_mode freq ties]=mode(consensus_mat_iter_noyeo) %can also look at nodes that still aren't able to be assigned , which nodes are have the most variance across nodes in the coocurrence matrix.
%save out
outdir='/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/'
outfile=(fullfile(outdir, strcat('n670_k',num2str(k),'training_sample_consensus_partition.mat')))
save(outfile, 'consensus_iter_mode', 'qpc')

%% Calculate log-evidence of this consensus partition on the replication sample
%load subjlist
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'

%load in the FC matrix for each subject, estimate the log-evidence for the
%consensus partition at that k
LP=zeros(size(test_subjlist,1),1);
LE=zeros(size(test_subjlist,1),1);
for n=1:size(test_subjlist,1)
    sub=test_subjlist{n,:}
    try 
        %% Load FC matrix
        fcfile=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
        subfcmat = load(fcfile);

        %perform the following line just to get a basic model structure, run only
        %for a single trial
        [S(n,:),M] = wsbm(subfcmat,k,'W_Distr','Normal','numTrials',1,'E_Distr','None');

        %want to get the log-likelihood for some partition fitting the network A
        LP(n) = calc_nullEvidence(M,consensus_iter_mode,subfcmat);

    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end

%save outfile
outdir='/cbica/projects/spatial_topography/data/imageData/wsbm/site14site20_test_sample/brains/'
outfile=dataset(test_subjlist, LE,LP)
export(outfile,'File',fullfile(outdir,strcat('n544_test_sample_',num2str(k),'consensus_log_evidence_replicate.csv')),'Delimiter',',')

end
