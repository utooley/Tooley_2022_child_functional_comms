
%Running on the cluster
datadir=fullfile('/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
listdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/parcellation'
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zNetworks'
z_avg_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
wsbm_dir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample'
addpath(genpath('/data/picsl/mackey_group/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal'))
addpath(genpath('/home/utooley/matlab/WSBM_v1.2'))
addpath(genpath('/home/utooley/matlab/BCT'))
addpath(genpath('/home/utooley/matlab/NCT'))
addpath(genpath('/home/utooley/matlab/'))
addpath(genpath('/home/utooley/matlab/NIfTI_20140122/'))

%running with the cluster mounted locally
%CHECK WHERE ON THE CLUSTER IS MOUNTED EXACTLY-THAT MAY MESS THIS UP
datadir=fullfile('~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
listdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/parcellation'
z_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zNetworks'
noz_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400Networks'
z_avg_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
wsbm_dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample'

%get the subject list,excluding those who have NAs
%subjlist=readtable(fullfile(listdir,'n27_cohort_file_one_run_only_21019.csv'),'Delimiter',',')
subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess.csv')) %the first two runs here are those that were input into gwMRF

wsbm_dir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample'
%on cluster
wsbm_dir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/'

k=7

n=2
sub=char(subjlist.id(n)) %look at this
file=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
subfcmat=load(file);
[Labels Model]=wsbm(subfcmat, k,'E_Distr','None', 'W_Distr', 'Normal', 'numTrials', 50, 'verbosity', 1, 'alpha', 0, 'parallel', 1);
Labels_orig=Labels; %run WSBM twice, see variation in partitions
Model_orig=Model;
l=pair_labeling(Labels_orig, Labels)
zrand(Labels, Labels_orig) %does relabeling for persistence make a difference to zRand?
zrand(
randind = randperm(size(subfcmat,1)) 
[~,previous_order]=sort(randind)
shuff_subfcmat = subfcmat(randind,randind)
[Labels Model]=wsbm(shuff_subfcmat, k,'E_Distr','None', 'W_Distr', 'Normal', 'numTrials', 50, 'verbosity', 1, 'alpha', 0, 'parallel', 1);
old_subfcmat=shuff_subfcmat(previous_order,previous_order)
shuff_Labels=Labels(previous_order)
%compare Labels_orig to shuff_Labels (they may have different numbers, but
%should be similar)

l=multislice_pair_labeling([Labels Labels_orig])

%% LOOK AT WHETHER THE ORDER OF RELABELING MATTERS

yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
yeo_nodes=dlmread('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')

%read in each subject's wsbm partition
for n=1:height(subjlist)
    sub=char(subjlist.id(n));
    %load the wsbm for a given subject
    file=fullfile(wsbm_dir,strcat(sub,'_wsbm.mat'));
    try 
        load(file);
        sub_log_evidence(n,1)=Model.Para.LogEvidence;
        %a row for each subject's labels for the k=7 partition
        sub_labels(n,:)=Labels; %skip this, before relabeling for
        %persistence with Yeo these are not right
        %Look at how many unique labels were output, since Rick says these might be
        %different
        sub_unique_labels(n,1)=numel(unique(Labels));
        models(n, :)=Model;
        partition=Labels;%load in the partition from the WSBM for this subject
        %relabel WSBM to fit numbers label to Yeo partition so can compare stats
        %relabel for similarity to Yeo partition
        %temp=multislice_pair_labeling([yeo_nodes partition]);
        %sub_labels(n,:)=temp(:,2)';
        %partition=temp(:,2);
        %generate one big partition matrix
        part_matrix(:,n)=partition;
        %save how similar a subject is to Yeo partition
        sub_similarity_to_yeo(n,1)=zrand(partition, yeo_nodes);
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end
%% SAVE OUTFILES

outfile=dataset(subjlist.id, sub_log_evidence, sub_unique_labels, sub_similarity_to_yeo,sub_labels)
export(outfile,'File',fullfile(wsbm_dir,'wsbm_k7_training_sample_n670_site16_50trials_norelabel.csv'),'Delimiter',',')

outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling'
%outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling'
outfile=(fullfile(outdir, 'n670_training_sample_no_relabeling_part_matrix_and_partitions.mat'))
save(outfile, 'part_matrix', 'yeo_opt_part_matrix','noyeo_opt_part_matrix', 'consensus_mat_noyeo', 'consensus_mat_yeo', 'yeo_mode', 'noyeo_mode')

outfile=(fullfile(outdir, 'n670_training_sample_norelabel_nodal_var.mat'))
save(outfile, 'node_var')
outfile=(fullfile(outdir, 'n670_training_sample_norelabel_nodal_entropy.mat'))
save(outfile, 'e')
outfile=(fullfile(outdir, 'consensus_mat.mat'))
save(outfile, 'consensus_mat')
%% Load in data and examine entropy and partitions
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling'
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling'
outfile=(fullfile(outdir, 'n670_training_sample_no_relabeling_part_matrix_and_partitions.mat'))
load(outfile)
outfile=(fullfile(outdir, 'consensus_iter_mode.mat'))
load(outfile)

%either use multislice_pair_label or just pair_label
%input=400 by p (nodes by partitions) matrix, each subj is a column, n
%subjects
yeo_opt_part_matrix=multislice_pair_labeling([yeo_nodes part_matrix]) %relabel starting with Yeo
yeo_opt_part_matrix=yeo_opt_part_matrix(:,2:671); %take out the first column which is Yeo nodes.
noyeo_opt_part_matrix=multislice_pair_labeling(part_matrix) %relabel just across subjects
%opt_part_matrix=part_matrix

%The modal matrix is dependent on the ordering of the subjects when
%relabeling. Trying to get around that by averaging!
num_iterations=10000
idx=[1:670]
initial=noyeo_mode
running_sum=zeros(400,1)
for i=1:num_iterations
noyeo_node_relabeled=multislice_pair_labeling(part_matrix(:,idx(randperm(length(idx))))); %permute it once
noyeo_node_mode_perm=mode(noyeo_node_relabeled,2);
relabeled=multislice_pair_labeling([initial noyeo_node_mode_perm]); %relabel it to an initial modal partition
running_sum=running_sum+relabeled(:,2);
end
[mode freq ties]=mode(noyeo_node_mode_perm,2)
avg_modal_partition2=round(running_sum./num_iterations)
%try relabeling at every go-around and then just taking the mode
for i=1:size(part_matrix,2)
noyeo_node_relabeled=multislice_pair_labeling([initial part_matrix(:,i)]); %permute it once
noeyo_node_relabeled_all(:,i)=noyeo_node_relabeled(:,2);
end

%modal community from these relabeled matrices
yeo_mode=mode(yeo_opt_part_matrix,2)
noyeo_mode=mode(noyeo_opt_part_matrix, 2)

%consensus similarity community from these relabeled matrices
consensus_mat_noyeo=consensus_similarity(noyeo_opt_part_matrix')'; %maybe transposed?
consensus_mat_yeo=consensus_similarity(yeo_opt_part_matrix')'; %maybe transposed?

%consensus iterative community (does it matter if you relabel or not? Not more than it varies with iterating over mod max)
gamma=2 %when keeping the thresholding of the consensus matrix, this is the right gamma
[consensus_mat_iter_noyeo Q2 X_new3 qpc cooccurence_matrix]=consensus_iterative(noyeo_opt_part_matrix', gamma);
[consensus_iter_mode freq ties]=mode(consensus_mat_iter_noyeo) %can also look at nodes that still aren't able to be assigned , which nodes are have the most variance across nodes in the coocurrence matrix.

%looking at the frequencies of assignments in consensus iterative across Yeo communities
freq_continuous=abs(freq-670)' %get the absolute value of how many times less than 670 it was assigned to the same community
comms=unique(yeo_nodes)
for i = 1:length(comms) % for each of the yeo communities
    Wi = yeo_nodes == comms(i); % find index for parcels with that community
    Wv_temp = freq_continuous(Wi,1); % extract node var for those parcels
    %average it
    num_ties_consensus_iter_by_yeocomm(i)=mean(Wv_temp); 
end

%need to figure out entropy on the co-occurence matrix--not sure this is
%right!
p = bsxfun(@rdivide,cooccurence_matrix,sum(cooccurence_matrix));  % probabilities
p = bsxfun(@rdivide,cooccurence_matrix,670);  % probabilities
p2=p;
p2(p==1)=NaN;
node_entropy = -nansum(p.*log2(p))';        % entropy

%entropy and variance
node_var=var(cooccurence_matrix, 0,2) %calculate node-wise SD (how much the labeling varies across subjects after each is relabeled for Yeo)
h = hist(noyeo_opt_part_matrix(:,:)',7); % node histogram
p = bsxfun(@rdivide,h,sum(h));  % probabilities
node_entropy = -nansum(p.*log2(p))';        % entropy

%see z-score of the Rand coefficient for the consensus representative
%partition
zrandconsensus=zrand(yeo_nodes,consensus_mat_yeo)
zrandconsensus_noyeo=zrand(yeo_nodes,consensus_mat_noyeo)

%see z-score of the Rand coefficient for the consensus iterative partition
zrandconsensus_iter=zrand(yeo_nodes,consensus_iter_mode)

%see z-score of the Rand coefficient for the modal community
zrandmodal_yeo=zrand(yeo_nodes,yeo_mode)
zrandmodal_noyeo=zrand(yeo_nodes,noyeo_mode)

%relabel the modal partitions so that they are visually comparable to Yeo
relabeled_yeo_mode=multislice_pair_labeling([yeo_nodes yeo_mode])
relabeled_noyeo_mode=multislice_pair_labeling([yeo_nodes noyeo_mode])

%entropy and variance
node_var=var(noyeo_opt_part_matrix, 0,2) %calculate node-wise SD (how much the labeling varies across subjects after each is relabeled for Yeo)
h = hist(noyeo_opt_part_matrix(:,:)',7); % node histogram
p = bsxfun(@rdivide,h,sum(h));  % probabilities
node_entropy2 = -nansum(p.*log2(p))';        % entropy

comms=unique(yeo_nodes)
for i = 1:length(comms) % for each of the yeo communities
    Wi = yeo_nodes == comms(i); % find index for parcels with that community
    Wv_temp = node_var(Wi,1); % extract node var for those parcels
    %average it
    var_by_yeocomm(i)=mean(Wv_temp); 
end
for i = 1:length(comms) % for each of the yeo communities
    Wi = yeo_nodes == comms(i); % find index for parcels with that community
    Wv_temp = node_entropy(Wi,1); % extract node var for those parcels
    %average it
    entropy_by_yeocomm(i)=mean(Wv_temp); 
end

%% Relabel the consensus representative and group partitions so they are visually comparable to Yeo
%relabel the consensus partitions so that they are visually comparable to Yeo
consensus_iter_mode_yeorelabeled=multislice_pair_labeling([yeo_nodes consensus_iter_mode']);
consensus_iter_mode_yeorelabeled=consensus_iter_mode_yeorelabeled(:,2);
consensus_represent_yeorelabeled=multislice_pair_labeling([yeo_nodes consensus_mat_noyeo])
consensus_represent_yeorelabeled=consensus_represent_yeorelabeled(:,2)

outfile=(fullfile(outdir, 'n670_training_sample_consensus_partitions_yeorelabeled.mat'))
save(outfile, 'consensus_iter_mode_yeorelabeled', 'consensus_represent_yeorelabeled')
%% Save a .nii file so can see WSBM on the brain 
%read into nifti?
templateVolume = '/data/picsl/mackey_group/tools/schaefer400/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_1mm.nii.gz'
nii = load_nii(templateVolume);
image = double(nii.img);
spacing = nii.hdr.dime.pixdim(2:4);

%% Assign labels to brain regions
outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling'
% read in the mapping of template nifti voxel label numbers to actual brain regions
mapping = readtable('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7NodeIndex.1D.txt','ReadVariableNames',false);
%mapping = readtable('~/Documents/tooleyEnviNetworks/parcels/Glasser/glasser_lookup.csv');
mapping_nums=table2array(mapping(:,1));
brains=struct('yeo_mode',yeo_mode, 'noyeo_mode', noyeo_mode, 'consensus_mat_noyeo', consensus_mat_noyeo, 'consensus_mat_yeo', consensus_mat_yeo);
brains=struct('avg_modal_partition',avg_modal_partition, 'avg_modal_partition2', avg_modal_partition2, 'avg_modal_partition_mode_mode', mode);

names=fieldnames(brains);
for x=1:numel(names)
    parc = zeros(size(image)); %create a new image matrix
    labels=brains.(names{x})
    for r=1:length(labels) %go through all possible voxel values and assign labels
        parc(find(image==mapping_nums(r))) = labels(r);
    end
    % WRITE OUTFILE
    parc = double(parc);
    orig = nii.hdr.hist.originator; %get the origin of the original image
    orig = orig(1:3);
    niiNew = make_nii(parc,spacing,orig); %write out the new nifti
    niiNew.hrd.dime.bitpix=16; %set the datatype
    save_nii(niiNew,strcat('/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/testing_relabeling/wsbm_k7_n670_site16_',char(names(x)),'.nii'));
end

%% COMPARE TO CONSENSUS PARTITION FOR I DID FIRST, RELABELING TO YEO AT SUBJECT LEVEL
outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/site16_training_sample/brains/'
outfile=(fullfile(outdir, 'n670_training_sample_relabel_first_consensus_mat_and_nodal_variance.mat'))
load(outfile)

%compare the optimal partition matrix to those generated without relabeling
%at the subject level
opt_part_matrix

zrand(node_mode,yeo_mode)
zrand(node_mode,noyeo_mode)
zrand(yeo_nodes,node_mode)
zrand(yeo_nodes,yeo_mode)
zrand(yeo_nodes,noyeo_mode)

zrand(consensus_mat,consensus_mat_noyeo)
zrand(consensus_mat,consensus_mat_yeo)

zrand(yeo_nodes,consensus_mat)
zrand(yeo_nodes,consensus_mat_noyeo)

multislice_pair_labeling([consensus_mat' consensus_mat_yeo])
