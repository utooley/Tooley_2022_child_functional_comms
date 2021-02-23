%PATHS
datadir='~/Downloads/Tooley_2021_child_functional_comms/data/example/'
listdir='~/Downloads/Tooley_2021_child_functional_comms/data/example/'
subjlist='~/Downloads/Tooley_2021_child_functional_comms/data/example/n2_example_subjlist' %2 subjects example data
wsbm_dir='~/Downloads/Tooley_2021_child_functional_comms/data/example/wsbm_output/'
addpath(genpath('~/Documents/tools/MATLAB/WSBM_v1.2')) %add toolboxes
addpath(genpath('~/Documents/tools/MATLAB/NCT')) %add toolboxes
addpath(genpath('~/Documents/tools/MATLAB/GenLouvain-2.1'))

%get the subject list
subjlist=readtable(fullfile(listdir,'n2_example_subjlist.csv'))

%% Run WSBM with k=7

k=7
for n=1:height(subjlist)
    sub=char(subjlist.id(n)) %look at this
    outfile=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'))
    if (exist(outfile)==2) %if it's already written don't do it again
        fprintf('Sub %s already exists. \n', sub);
    else
    file=fullfile(datadir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
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
export(outfile,'File',fullfile(wsbm_dir,'wsbm_50trials.csv'),'Delimiter',',')

%% Relabel community partitions, create a consensus partition
yeo_nodes_schaefer_label=dlmread('~/Downloads/Tooley_2021_child_functional_comms/data/example/schaefer400x7CommunityAffiliation.1D.txt')
clear part_matrix

%read in each subject's wsbm partition
for n=1:height(subjlist)
    sub=char(subjlist.id(n));
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'));
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

%% Create consensus partition
%for the relabeling step, you need a reasonably large number of rows as input to
%the iterative consensus algorithm, so we create a "200-subject" partition matrix
part_matrix=repmat(part_matrix,1,100);
%use consensus iterative procedure (no relabeling is necessary before this
%step)
gamma=2 %when keeping the thresholding of the consensus matrix as default, this is the right gamma
[consensus_mat_iter_noyeo Q2 X_new3 qpc cooccurence_matrix]=consensus_iterative(part_matrix', gamma);
[consensus_iter_mode freq ties]=mode(consensus_mat_iter_noyeo) %extract assignment variance across optimizations

%save these outfiles
outfile=('~/Downloads/Tooley_2021_child_functional_comms/data/example/wsbm_output/consensus_iter_freq.mat')
save(outfile, 'freq')

%% Relabel and save the consensus partitions so they are visually comparable to Yeo7
outdir='~/Downloads/Tooley_2021_child_functional_comms/data/example/wsbm_output'
%relabel the consensus partition so that it is visually comparable to Yeo
%community assignments with Hungarian algorithm
consensus_iter_mode_yeorelabeled=multislice_pair_labeling([yeo_nodes_schaefer_label consensus_iter_mode']);
consensus_iter_mode_yeorelabeled=consensus_iter_mode_yeorelabeled(:,2);

outfile=(fullfile(outdir, 'n2_example_consensus_partition_yeorelabeled2.mat'))
save(outfile, 'consensus_iter_mode_yeorelabeled')
