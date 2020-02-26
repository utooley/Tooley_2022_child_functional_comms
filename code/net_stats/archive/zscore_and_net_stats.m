%Running on the cluster
datadir=fullfile('/data/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_despike_onerun')
listdir='/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/'
outdir='/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/Schaefer400zNetworks'

%running with the cluster mounted locally
%CHECK WHERE ON THE CLUSTER IS MOUNTED EXACTLY-THAT MAY MESS THIS UP
datadir=fullfile('~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_despike_onerun')
listdir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/'
z_outdir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/Schaefer400zNetworks'
noz_outdir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/Schaefer400Networks'

%get the subject list,excluding those who have NAs
subjlist=readtable(fullfile(listdir,'n27_cohort_file_one_run_only_21019.csv'),'Delimiter',',')
subjlist=subjlist(:,1);
%% Z-score FC matrices
for n=1:height(subjlist)
    sub=char(subjlist.id0(n)) %look at this
    file=fullfile(datadir,strcat(sub,'/fcon/schaefer400/',sub,'_schaefer400_network.txt'));
    subfcmat=load(file);
    %make into adjacency matrix and save out
    size_vec=tril(ones(400,400),-1);
    adj_mat=size_vec;
    adj_mat(adj_mat==1)=subfcmat;
    subfcmat=adj_mat+adj_mat';
    outfile=fullfile(noz_outdir,strcat(sub,'_schaefer400_network.txt'));
    csvwrite(outfile, subfcmat);
    %elimate parcel 52 (parcel index 103), delete row 103
    %subfcmat=removerows(subfcmat, 'ind', [103]);
    %remove column 103
   % subfcmat(:,103)=[]; %never checked parcel coverage for this.
    %replace the diagonal of 1's with 0's
    for x=1:359
        subfcmat(x,x)=0;
    end
    %create an empty z-matrx
    zfc=[];
    for i=1:400
        %cycle through each column of the FC matrix and do a fisher r-to-z
        %for each value
        zfc(:,i)=fisherz(subfcmat(:,i));
    end
    outfile=fullfile(z_outdir, strcat(num2str(sub),'_Schaefer400subjectspace_znetwork.txt'));
    csvwrite(outfile, zfc);
end

%% Within and between network connectivity
%cluster mounted locally
datadir='~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/Schaefer400zNetworks'

for n=1:height(subjlist)
    sub=char(subjlist.id0(n)) %look at this
    file=fullfile(datadir,strcat(num2str(sub),'_Schaefer400subjectspace_znetwork.txt'))
    subfcmat = load(file);
    for x=1:359
        subfcmat(x,x)=0;
    end
%% SIGNED MATRICES 
%average network strength (the mean of all network weights in the matrix that are not equal to
%0))
avgweight(n,1)=mean(subfcmat(subfcmat~=0));

%Clustering coefficient sing Constantini & Perugini's generalization of the Zhang &
%Horvath formula (option3). This formula takes both positive &
%negative weights into account simultaneously, produces 1 value for each
%node. Then take the mean.
avgclustco_both(n,1)=mean(clustering_coef_wu_sign(subfcmat,3));

%Path length
%MUST CONVERT TO A CONNECTION_LENGTHS matrix
%%L and then a distance matrix D
connlenmat=weight_conversion(subfcmat, 'lengths');
distmat=distance_wei(connlenmat);
cpathleng(n,1)=charpath(distmat);

%Strength 
[Spos Sneg]=strengths_und_sign(subfcmat);
strengths_pos(n,:)=Spos; %nodal strength of positive weights
strengths_neg(n,:)=Sneg;    %nodal strength of negative weights

%% SIGNED NULL MODELS
% NULL ONE %
%take a given subject's network, randomize it, 
%this function preseves the strength and degree function in weighted networks
net=null_model_und_sign(subfcmat);

avgweight_null1(n,1)=mean(net(net~=0));

%Clustering coefficient sing Constantini & Perugini's generalization of the Zhang &
%Horvath formula (option3). This formula takes both positive &
%negative weights into account simultaneously, produces 1 value for each
%node. Then take the mean.
avgclustco_both_null1(n,1)=mean(clustering_coef_wu_sign(net,3));

%Path length
%MUST CONVERT TO A CONNECTION_LENGTHS matrix
%%L and then a distance matrix D
connlenmat=weight_conversion(net, 'lengths');
distmat=distance_wei(connlenmat);
cpathleng_null1(n,1)=charpath(distmat);

%Strength 
[Spos Sneg]=strengths_und_sign(net);
strengths_pos_null1(n,:)=Spos;
strengths_neg_null1(n,:)=Sneg;

% NULL TWO %
%take a given subject's network, randomize it, 
%this function DOES NOT preserve the strength function in weighted networks
%random, directed network with a specified number of fully connected modules linked together by evenly distributed
%remaining random connections. (num vertices, edges, modules)
%net=makeevenCIJ(2^16,400,8); %
net=randmio_und_signed(subfcmat, 30);

avgweight_null2(n,1)=mean(net(net~=0));

%Clustering coefficient sing Constantini & Perugini's generalization of the Zhang &
%Horvath formula (option3). This formula takes both positive &
%negative weights into account simultaneously, produces 1 value for each
%node. Then take the mean.
avgclustco_both_null2(n,1)=mean(clustering_coef_wu_sign(net,3));

%Path length
%MUST CONVERT TO A CONNECTION_LENGTHS matrix
%%L and then a distance matrix D
connlenmat=weight_conversion(net, 'lengths');
distmat=distance_wei(connlenmat);
cpathleng_null2(n,1)=charpath(distmat);

%Strength 
[Spos Sneg]=strengths_und_sign(net);
strengths_pos_null2(n,:)=Spos;
strengths_neg_null2(n,:)=Sneg;
%% THRESHOLDED MATRICES
thrsubfcmat=threshold_absolute(subfcmat, 0);
avgweight_thresh(n,1)=mean(thrsubfcmat(thrsubfcmat~=0));

%Path length
%MUST CONVERT TO A CONNECTION_LENGTHS matrix
%%L and then a distance matrix D
connlenmat=weight_conversion(thrsubfcmat, 'lengths');
distmat=distance_wei(connlenmat);
cpathleng_thresh(n,1)=charpath(distmat);

%% Thresholded null network models
% NULL ONE %
%take a given subject's network, randomize it, 
%this function preseves the strength and degree function in weighted networks
net=null_model_und_sign(thrsubfcmat);

avgweight_thresh_null1(n,1)=mean(net(net~=0));
%Path length
connlenmat=weight_conversion(net, 'lengths');
distmat=distance_wei(connlenmat);
cpathleng_thresh_null1(n,1)=charpath(distmat);
% NULL TWO %
%take a given subject's network, randomize it, 
net=randmio_und_signed(thrsubfcmat, 30);

avgweight_thresh_null2(n,1)=mean(net(net~=0));
%Path length
connlenmat=weight_conversion(net, 'lengths');
distmat=distance_wei(connlenmat);
cpathleng_thresh_null2(n,1)=charpath(distmat);

end

%% Export data
outdir= '~/Dropbox/'
outfile=dataset(char(subjlist.id0), avgweight, avgclustco_both, cpathleng, avgweight_null1, avgclustco_both_null1, cpathleng_null1,avgweight_null2, avgclustco_both_null2, cpathleng_null2, avgweight_thresh, cpathleng_thresh, avgweight_thresh_null1, cpathleng_thresh_null1, avgweight_thresh_null2, cpathleng_thresh_null2)
strengthoutfile=dataset(strengths_pos, strengths_neg)
strengthsnullone=dataset(strengths_pos_null1, strengths_neg_null1)
strengthsnulltwo=dataset(strengths_pos_null2, strengths_neg_null2)

export(outfile,'File',fullfile(outdir,'n27_one_run_only_net_meas_Schaefer400.csv'),'Delimiter',',')
export(strengthoutfile,'File',fullfile(outdir,'n27_one_run_only_strengths.csv'),'Delimiter',',')
export(strengthsnullone,'File',fullfile(outdir,'n27_one_run_only_strengths_null1.csv'),'Delimiter',',')
export(strengthsnulltwo,'File',fullfile(outdir,'n27_one_run_only_strengths_null2.csv'),'Delimiter',',')
%% SYSTEM SEGREGATION CODE-UNUSED

system_segreg=zeros(height(subjlist),1);
mean_within_sys=zeros(height(subjlist),1);
mean_between_sys=zeros(height(subjlist),1);


%average whole system segregation, from Micaela Chan 2018
[S, W, B] = segregation(subfcmat,yeo_nodes);
system_segreg(n,1)=S;
mean_within_sys(n,1)=W;
mean_between_sys(n,1)=B;

%Within and between connectivity, adapted from Micaela Chan 2018
Ci=yeo_nodes;
nCi = unique(Ci);
M=subfcmat;

for i = 1:length(nCi) % loop through communities
    for j = 1:length(nCi)
       Wi = Ci == nCi(i); % find index for within communitiy edges
       Bi = Ci == nCi(j); % find index for between communitiy edges to specific community
       
       Wv_temp = M(Wi,Wi); % extract within communitiy edges
       Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community
       
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
system_conn_vector = reshape(system_connectivity',[],1)';

system_conn(n,:)=system_conn_vector;

end
