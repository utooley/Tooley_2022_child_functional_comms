%Running on the cluster
datadir=fullfile('/data/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_despike_onerun')
listdir='/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/'
outdir='/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/Schaefer400zNetworks'

%running with the cluster mounted locally
%CHECK WHERE ON THE CLUSTER IS MOUNTED EXACTLY-THAT MAY MESS THIS UP
datadir=fullfile('~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_36p_gsr_multruns')
listdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/parcellation'
z_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/Schaefer400zNetworks'
noz_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/Schaefer400Networks'
wsbm_dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/full_sample'

%get the subject list,excluding those who have NAs
%subjlist=readtable(fullfile(listdir,'n27_cohort_file_one_run_only_21019.csv'),'Delimiter',',')
subjlist=readtable(fullfile(listdir,'n63_train_sample.txt'))

subjlist=subjlist(:,1);

%% Make a consensus model from only the training sample
%load each training subject's labels for Schaefer 400 parcels from WSBM k=7
for n=1:height(subjlist)
    sub=subjlist{n,:}
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub{1},'_wsbm.mat'));
    try 
        load(file);
        sub_log_evidence(n,1)=Model.Para.LogEvidence;
        %a row for each subject's labels for the k=7 partition
        part_matrix(:,n)=Labels;
    catch
        fprintf('Cant read sub %s, skipped. \n', sub{1});
    end
end
%% Relabel community partitions to be most parsimonious across subjects, create a consensus partition
%either use multislice_pair_label or just pair_label
%input=n by p (nodes by partitions) matrix, each subj is a column
opt_part_matrix=multislice_pair_labeling(part_matrix)

% Create a consensus partition
consensus_mat=consensus_similarity(part_matrix'); %maybe transposed?
%see z-score of the Rand coefficient for the consensus community vs. Yeo
%communities
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
consensus_wsbm_yeo=multislice_pair_labeling([yeo_nodes consensus_mat'])
consensus_wsbm_yeo=consensus_wsbm_yeo(:,2)
zrandconsensus=zrand(yeo_nodes,consensus_mat) 

node_var=var(opt_part_matrix, 0,2) %calculate node-wise SD (how much the labeling varies across subjects)
node_entropy=Entropy(opt_part_matrix') %first entropy from Entropy.m
[node_mode freq ties]=mode(opt_part_matrix, 2)%calculate nodewise mode (most common label), and export ties as well.

%% Look at node variance by Yeo network
dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/training_sample/'
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')

node_var=load(fullfile(dir,'node_var.mat'));
comms=unique(yeo_nodes)
for i = 1:length(comms) % for each of the yeo communities
    Wi = yeo_nodes == comms(i); % find index for parcels with that community
    Wv_temp = node_var.node_var(Wi,1); % extract node var for those parcels
    %average it
    var_by_yeocomm(i)=mean(Wv_temp); 
end
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/training_sample/'
export(dataset(var_by_yeocomm),'File',fullfile(outdir,'variance_in_wsbm_assign_by_yeonet.csv'),'Delimiter',',')
%% SAVE OUTFILES
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/training_sample/'
outfile=(fullfile(outdir, 'optim_partition_matrix.mat'))
save(outfile, 'opt_part_matrix')
outfile=(fullfile(outdir, 'node_var.mat'))
save(outfile, 'node_var')
outfile=(fullfile(outdir, 'consensus_mat.mat'))
save(outfile, 'consensus_wsbm_yeo')

%% Save a .nii file so can see WSBM on the brain 
%read into nifti?
templateVolume = '~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400MNI.nii.gz'
nii = load_nii(templateVolume);
image = double(nii.img);
spacing = nii.hdr.dime.pixdim(2:4);

%% read in the mapping of template nifti voxel label numbers to actual brain regions
mapping = readtable('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400NodeIndex.1D.txt','ReadVariableNames',false);
%mapping = readtable('~/Documents/tooleyEnviNetworks/parcels/Glasser/glasser_lookup.csv');
mapping_nums=table2array(mapping(:,1));

parc = zeros(size(image)); %create a new image matrix
labels=node_var %use the most commonly assigned community or the variance as the label
labels=node_entropy'
labels=node_mode %or consensus mat
labels=consensus_wsbm_yeo' %or consensus mat
%% ASSIGN LABELS TO NODE REGIONS

for r=1:length(labels) %go through all possible voxel values and assign labels
    parc(find(image==mapping_nums(r))) = labels(r);
end
%% WRITE OUTFILE

parc = double(parc);
orig = nii.hdr.hist.originator; %get the origin of the original image
orig = orig(1:3);
niiNew = make_nii(parc,spacing,orig); %write out the new nifti
niiNew.hrd.dime.bitpix=16; %set the datatype
save_nii(niiNew,'~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/training_sample/n63_training_sample_node_var.nii.gz');
save_nii(niiNew,'~/Downloads/test.nii.gz')
%% Get WSBM stats out using training partition on test sample
datadir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/wsbm/full_sample'
wsbm_dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/training_sample/'
consensus_wsbm_yeo=load(fullfile(wsbm_dir, 'consensus_mat.mat'))
consensus_wsbm_yeo=consensus_wsbm_yeo.consensus_wsbm_yeo;
%read in the yeo partition
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
subjlist=readtable(fullfile(listdir,'n64_test_sample.txt'))

subjlist=subjlist(:,1);
%preallocate variables
deviation_edge_weights_wsbm=zeros(height(subjlist),1);
deviation_edge_weights_yeo=zeros(height(subjlist),1);

%load in the models for each subject
for n=1:height(subjlist)
    sub=subjlist{n,:}
    %load in their functional data
    try 
        fcfile=fullfile(z_outdir,strcat(sub{1},'_run-01_Schaefer400subjectspace_znetwork.txt'));
        subfcmat = load(fcfile); 
        %calculate some output measures
        %average edge weight for the subject
        avgweight(n,1)=mean(subfcmat(subfcmat~=0));
        %examine the log-evidence of the consensus partition fitting this
        %subject, from Andrew sharing loglikelihood_example.m
        Data = subfcmat(:);
        w = 0.000001;
        [WDistr] = setup_distr('Normal',[mean(Data),var(Data)+mean(Data)^2,w],1); %do I need to do this at other times?

        %perform the following line just to get a basic model structure, run only
        %for a single trial
        [S,M] = wsbm(subfcmat,7,'W_Distr','Normal','numTrials',1,'E_Distr','None');
        log_evidence(n) = calc_logEvidence(M);
        %calc_nullEvidence(M,consensus_wsbm_yeo,subfcmat); %WAITING ON
        %ANDREW TO ANSWER
        %apply the consensus partition, should already be maximally relabeled
        %for consensus with Yeo
        %Within and between-block connectivity--use segregation code below to get this.
        [S, W, B] = segregation(subfcmat,consensus_wsbm_yeo');
        system_segreg_wsbm(n,1)=S;
        mean_within_sys_wsbm(n,1)=W;
        mean_between_sys_wsbm(n,1)=B;
        %Connectivity between each of the 7 blocks
        Ci=consensus_wsbm_yeo;
        nCi = unique(Ci);
        M=subfcmat;

        for i = 1:length(nCi) % loop through communities
            for j = 1:length(nCi)
               Wi = Ci == nCi(i); % find index for within communitiy edges
               Bi = Ci == nCi(j); % find index for between communitiy edges to specific community

               Wv_temp = M(Wi,Wi); % extract within communitiy edges
               Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community
                
               %calculate standard deviation of edges within blocks here
               deviation_edge_weights_wsbm(n,1)=deviation_edge_weights_wsbm(n,1)+std(Wv_temp(Wv_temp~=0))^2; %Gu paper calculation, finished below
               %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
               Bv = [Bv_temp(:)'];
               system_connectivity(i,j)=mean(Bv(Bv~=0));
               %system_between(i,1)=mean(Bv);
               %if i==j
               %else
            end
        end

        %transpose these (or something so they can be saved out on a subject basis.
        system_connectivity;
        block_conn_vector = reshape(system_connectivity',[],1)';

        block_connectivity(n,:)=block_conn_vector;
        %finish the std deviation of edge weights calculation
        deviation_edge_weights_wsbm(n,1)=sqrt(deviation_edge_weights_wsbm(n,1)/7);
        %some net stats on the WSBM partition
        %Participation coefficient average with WSBM partition
        [Ppos Pneg]=participation_coef_sign(subfcmat, consensus_wsbm_yeo);
        %each row is a subject, with part coef for 359 nodes
        sub_partcoef_pos_wsbm(n,:)=mean(Ppos);
        sub_partcoef_neg_wsbm(n,:)=mean(Pneg);
        %% Apply Yeo partition and calculation statistics to compare
        [S, W, B] = segregation(subfcmat,yeo_nodes);
        system_segreg_yeo(n,1)=S;
        mean_within_sys_yeo(n,1)=W;
        mean_between_sys_yeo(n,1)=B;
        %Connectivity between each of the 7 blocks
        Ci=yeo_nodes;
        nCi = unique(Ci);
        M=subfcmat;

        for i = 1:length(nCi) % loop through communities
            for j = 1:length(nCi)
               Wi = Ci == nCi(i); % find index for within communitiy edges
               Bi = Ci == nCi(j); % find index for between communitiy edges to specific community

               Wv_temp = M(Wi,Wi); % extract within communitiy edges
               Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community
                
               %calculate standard deviation of edges within blocks here
               deviation_edge_weights_yeo(n,1)=deviation_edge_weights_yeo(n,1)+std(Wv_temp(Wv_temp~=0))^2; %Gu paper calculation, finished below
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
        yeo_conn_vector = reshape(system_connectivity',[],1)';
        yeo_connectivity(n,:)=yeo_conn_vector;
        %finish the std deviation of edge weights calculation
        deviation_edge_weights_yeo(n,1)=sqrt(deviation_edge_weights_yeo(n,1)/7);
        %some net stats on the WSBM partition
        %Participation coefficient average with WSBM partition
        [Ppos Pneg]=participation_coef_sign(subfcmat, yeo_nodes);
        %each row is a subject, with part coef for 359 nodes
        sub_partcoef_pos_yeo(n,:)=mean(Ppos);
        sub_partcoef_neg_yeo(n,:)=mean(Pneg);
        %compare overlap of the WSBM partition to the Yeo partition
        [zRand SR]=zrand(yeo_nodes, consensus_wsbm_yeo) 
        zRandstat(n,1)=zRand;
        SRstat(n,1)=SR;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub{1});
    end
end
%save outfile
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/test_sample'
outfile=dataset(subjlist{:,1}, avgweight, system_segreg_wsbm, mean_within_sys_wsbm, mean_between_sys_wsbm, log_evidence', deviation_edge_weights_wsbm, sub_partcoef_pos_wsbm, sub_partcoef_neg_wsbm, system_segreg_yeo, mean_within_sys_yeo, mean_between_sys_yeo, deviation_edge_weights_yeo, sub_partcoef_pos_yeo, sub_partcoef_neg_yeo, zRandstat, SRstat, block_connectivity, yeo_connectivity)
export(outfile,'File',fullfile(outdir,'n64_test_sample_wsbm_yeo_comparisons_logevidence_k7.csv'),'Delimiter',',')

%% Run WBSM model on Data-Derived Parcellation
outdir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/full_sample'
%add a loop for different values of k (communities)
%for k=1:12, for example
%or when running on the cluster can use the function below
%wsbm_function(sub, comm) %need to set outputs in the WSBM function so that it writes to the right place.
k=7

for n=1:height(subjlist)
    %sub=char(subjlist.id0(n)) %look at this
    sub=subjlist{n,:}; %only for one run
    outfile=fullfile(wsbm_dir,strcat(sub{1},'_wsbm.mat'))
    if (exist(outfile)==2) %if it's already written don't do it again
        fprintf('Sub %s already exists. \n', sub{1});
    else
    file=fullfile(z_outdir,strcat(sub,'_run-01_Schaefer400subjectspace_znetwork.txt'));
    try
    subfcmat=load(file{1});
    [Labels Model]=wsbm(subfcmat, k,'E_Distr','None', 'verbosity', 1, 'alpha', 0, 'parallel', 0);
    %logHw - Additive Log-likelihood constant for W_Distr
    sub_log_lik_for_weights(n,1)=Model.Data.logHw;
    %logHe - Additive Log-likelihood constant for E_Distr
    sub_log_lik_for_edges(n,1)=Model.Data.logHe;
    %LogEvidence - Marginal Log-likelihood (aka Log-Evidence), a model selection criterion 
    sub_log_evidence(n,1)=Model.Para.LogEvidence;

    %a row for each subject's labels for the k=7 partition
    sub_labels(n,:)=Labels;

    %Look at how many unique labels were output, since Rick says these might be
    %different
    sub_unique_labels(n,1)=numel(unique(Labels));
    models(n, :)=Model;

    %save the model
    outfile= fullfile(outdir,strcat(sub{1},'_wsbm.mat'))
    save(outfile, 'Model', 'Labels')
    catch
    fprintf('Cant read sub %s, skipped. \n', sub{1});
    end
    end
end

%save subject-level variables.
outfile=dataset(subjlist, sub_log_lik_for_weights, sub_log_lik_for_edges, sub_log_evidence, sub_unique_labels, sub_labels)
export(outfile,'File',fullfile(outdir,'wsbm_k7_n130.csv'),'Delimiter',',')
