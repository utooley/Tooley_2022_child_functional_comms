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
wsbm_dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/'

%get the subject list,excluding those who have NAs
%subjlist=readtable(fullfile(listdir,'n27_cohort_file_one_run_only_21019.csv'),'Delimiter',',')
subjlist=readtable(fullfile(listdir,'n64_train_sample.txt'))

subjlist=subjlist(:,1);
%% Z-score FC matrices
for n=1:height(subjlist)
    %sub=char(subjlist.id0(n)) %look at this
    sub=subjlist{n,:}; %only for one run
    file=fullfile(datadir,strcat(sub,'/run-01/fcon/schaefer400/',sub,'_run-01_schaefer400_network.txt'));
    try
    %subfcmat=load(file{1});
    subfcmat=load(file);
    %make into adjacency matrix and save out
    size_vec=tril(ones(400,400),-1);
    adj_mat=size_vec;
    adj_mat(adj_mat==1)=subfcmat;
    subfcmat=adj_mat+adj_mat';
    outfile=fullfile(noz_outdir,strcat(sub,'_run-01_schaefer400_network.txt'));
    csvwrite(outfile{1}, subfcmat);
    csvwrite(outfile, subfcmat);
    %elimate parcel 52 (parcel index 103), delete row 103
    %subfcmat=removerows(subfcmat, 'ind', [103]);
    %remove column 103
   % subfcmat(:,103)=[]; %never checked parcel coverage for this.
    %replace the diagonal of 1's with 0's
    for x=1:400
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc=[];
    for i=1:400
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(z_outdir, strcat(sub,'_run-01_Schaefer400subjectspace_znetwork.txt'));
    csvwrite(outfile{1}, zfc);
    catch
    fprintf('Cant read sub %s, skipped. \n', sub{1});
  end
end

%% Run WBSM model signed matrices
outdir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/'
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

%% Run some comparisons between the Yeo partition and the WSBM partition
datadir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/wsbm/'
 %read in the yeo partition
yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
subjlist=readtable(fullfile(listdir,'n130_one_site_0.2mm_nobogus.csv'))

subjlist=subjlist(:,1);
%preallocate variables
deviation_edge_weights_wsbm=zeros(height(subjlist),1);
deviation_edge_weights_yeo=zeros(height(subjlist),1);

%load in the models for each subject
for n=1:height(subjlist)
    sub=subjlist{n,:}
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub{1},'_wsbm.mat'));
    try 
        load(file);
        sub_log_evidence(n,1)=Model.Para.LogEvidence;
        %a row for each subject's labels for the k=7 partition
        sub_labels(n,:)=Labels;
        %Look at how many unique labels were output, since Rick says these might be
        %different
        sub_unique_labels(n,1)=numel(unique(Labels));
        models(n, :)=Model;
        %% Apply the calculated partition from the WSBM and do some statistics on it
        fcfile=fullfile(z_outdir,strcat(sub{1},'_run-01_Schaefer400subjectspace_znetwork.txt'));
        %parcel 52 is already gone
        subfcmat = load(fcfile);
        partition=Labels;%load in the partition from the WSBM for this subject
        %relabel WSBM to fit numbers label to Yeo partition so can compare stats
        temp=multislice_pair_labeling([yeo_nodes partition])
        partition=temp(:,2)
        %make the matrix for relabeling for optimal persistence
        part_matrix(:,n)=partition;
        %average edge weight for the subject
        avgweight(n,1)=mean(subfcmat(subfcmat~=0));
        %Within and between-block connectivity--use segregation code below to get this.
        [S, W, B] = segregation(subfcmat,partition);
        system_segreg_wsbm(n,1)=S;
        mean_within_sys_wsbm(n,1)=W;
        mean_between_sys_wsbm(n,1)=B;
        %Connectivity between each of the 7 blocks
        Ci=partition;
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
        [Ppos Pneg]=participation_coef_sign(subfcmat, partition);
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
        [zRand SR]=zrand(yeo_nodes, partition) 
        zRandstat(n,1)=zRand;
        SRstat(n,1)=SR;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub{1});
    end
end
%save outfile
outdir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/'
outfile=dataset(subjlist{:,1}, sub_log_evidence, avgweight, system_segreg_wsbm, mean_within_sys_wsbm, mean_between_sys_wsbm, deviation_edge_weights_wsbm, sub_partcoef_pos_wsbm, sub_partcoef_neg_wsbm, system_segreg_yeo, mean_within_sys_yeo, mean_between_sys_yeo, deviation_edge_weights_yeo, sub_partcoef_pos_yeo, sub_partcoef_neg_yeo, zRandstat, SRstat, block_connectivity, yeo_connectivity)
export(outfile,'File',fullfile(outdir,'n130_wsbm_yeo_comparisons_k7.csv'),'Delimiter',',')

%% Relabel community partitions to be most parsimonious across subjects, create a consensus partition
%either use multislice_pair_label or just pair_label
%input=n by p (nodes by partitions) matrix, each subj is a column
opt_part_matrix=multislice_pair_labeling(part_matrix)

node_var=var(opt_part_matrix, 0,2) %calculate node-wise SD (how much the labeling varies across subjects)
node_entropy=Entropy(opt_part_matrix') %first entropy from Entropy.m
[node_mode freq ties]=mode(opt_part_matrix, 2)%calculate nodewise mode (most common label), and export ties as well.

%% Create a consensus partition
consensus_mat=consensus_similarity(part_matrix'); %maybe transposed?
%see z-score of the Rand coefficient for the consensus community
zrandconsensus=zrand(yeo_nodes,consensus_mat) 
%% SAVE OUTFILES
outdir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/brains/'
outfile=(fullfile(outdir, 'optim_partition_matrix.mat'))
save(outfile, 'opt_part_matrix')
outfile=(fullfile(outdir, 'node_var.mat'))
save(outfile, 'node_var')
outfile=(fullfile(outdir, 'consensus_mat.mat'))
save(outfile, 'consensus_mat')

%transpose the re-labeled partition matrix so it can be saved with each
%subj on a row
% sub_labels_opt=opt_part_matrix'
% outfile=dataset(subjlist, sub_log_evidence, sub_unique_labels, sub_labels_opt)

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
labels=consensus_mat' %or consensus mat
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
save_nii(niiNew,'~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/brains/wsbm_k7_node_variance_in_assignment.nii');
save_nii(niiNew,'~/Desktop/wsbm_modal_node_assignment.nii');
save_nii(niiNew,'~/Desktop/wsbm_k7_community_assignments.nii');
save_nii(niiNew,'~/Desktop/wsbm_k7_community_assignments.nii');