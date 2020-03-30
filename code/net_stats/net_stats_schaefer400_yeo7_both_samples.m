% %Running on the cluster
% datadir=fullfile('/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
% listdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20/'
% z_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zNetworks'
% noz_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400Networks'
% z_avg_outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'
% addpath(genpath('/home/utooley/matlab/WSBM_v1.2'))
% addpath(genpath('/home/utooley/matlab/system_matrix_tools/'))
% addpath(genpath('/home/utooley/matlab/BCT/'))
% 
% %running with the cluster mounted locally
% %CHECK WHERE ON THE CLUSTER IS MOUNTED EXACTLY-THAT MAY MESS THIS UP
% datadir=fullfile('~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek')
% listdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20/'
% z_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zNetworks'
% noz_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400Networks'
% z_avg_outdir='~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/fc_matrices/site14site20_test_sample/Schaefer400zavgNetworks'
% 
% %get the subject list,excluding those who have NAs
% subjlist=readtable(fullfile(listdir,'n546_filtered_runs_site14site20_postprocess.csv')) %the first two runs here are those that were input into gwMRF
% %% Z-score FC matrices for the test sample
% for n=1:height(subjlist)
%     sub=char(subjlist.id(n)) %look at this
%     %sub=subjlist{n,:}; 
%     run1=char(subjlist.var1(n)) %%% FOR RUN 1 %%
%     file=fullfile(datadir,sub, run1, strcat('fcon/schaefer400x7/',sub,'_',run1,'_schaefer400x7_network.txt'));
%     try
%     %subfcmat=load(file{1});
%     subfcmat=load(file);
%     %make into adjacency matrix and save out
%     size_vec=tril(ones(400,400),-1);
%     adj_mat=size_vec;
%     adj_mat(adj_mat==1)=subfcmat;
%     subfcmat=adj_mat+adj_mat';
%     outfile=fullfile(noz_outdir,strcat(sub,'_',run1,'_schaefer400x7_network.txt'));
%     csvwrite(outfile, subfcmat);
%    % subfcmat(:,103)=[]; %never checked parcel coverage for this.
%     %replace the diagonal of 1's with 0's
%     for x=1:400
%         subfcmat(x,x)=0;
%     end
%     %create an empty z-matrx
%     zfc1=[];
%     for i=1:400
%         %cycle through each column of the FC matrix and do a fisher r-to-z
%         %for each value
%         zfc1(:,i)=fisherz(subfcmat(:,i));
%     end
%     outfile=fullfile(z_outdir, strcat(sub,'_', run1,'_Schaefer400x7_znetwork.txt'));
%     csvwrite(outfile, zfc1);
%     run2=char(subjlist.var2(n))  %%% FOR RUN 2 %%
%     file=fullfile(datadir,sub, run2, strcat('fcon/schaefer400x7/',sub,'_',run2,'_schaefer400x7_network.txt'));
%     subfcmat=load(file);
%     size_vec=tril(ones(400,400),-1); %make into adjacency matrix and save out
%     adj_mat=size_vec;
%     adj_mat(adj_mat==1)=subfcmat;
%     subfcmat=adj_mat+adj_mat';
%     outfile=fullfile(noz_outdir,strcat(sub,'_',run2,'_schaefer400x7_network.txt'));
%     csvwrite(outfile, subfcmat);
%     for x=1:400
%         subfcmat(x,x)=0;
%     end
%     %create an empty z-matrx
%     zfc2=[];
%     for i=1:400
%         %cycle through each column of the FC matrix and do a fisher r-to-z
%         %for each value
%         zfc2(:,i)=fisherz(subfcmat(:,i));
%     end
%     outfile=fullfile(z_outdir, strcat(sub,'_', run2,'_Schaefer400x7_znetwork.txt'));
%     csvwrite(outfile, zfc2);
%     %% AVERAGE THE TWO RUNS %%
%     aggregmat=(zfc2+zfc1)/2;
%     outfile=fullfile(z_avg_outdir, strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
%     csvwrite(outfile, aggregmat);
%     catch
%     fprintf('Cant read sub %s run %s, skipped. \n', sub);
%   end
% end
% 
% %% Run some partition quality metrics on the Yeo partition and estimate network statistics
% %read in the yeo partition
% yeo_nodes=dlmread('~/Desktop/cluster/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
% yeo_nodes=dlmread('/data/picsl/mackey_group/tools/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')
% subjlist=subjlist.id;
% %preallocate variables
% deviation_edge_weights_yeo=zeros(size(subjlist,1),1);
% modul_yeo=zeros(size(subjlist,1),1);
% 
% %load in the models for each subject
% for n=1:size(subjlist,1)
%     sub=subjlist{n,:}
%     try 
%         %% Load FC matrix
%         fcfile=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
%         %parcel 52 is already gone
%         subfcmat = load(fcfile);
%         %average edge weight for the subject
%         avgweight(n,1)=mean(subfcmat(subfcmat~=0));
%         %% Apply Yeo partition and calculate statistics with it
%         [S, W, B] = segregation(subfcmat,yeo_nodes);
%         system_segreg_yeo(n,1)=S;
%         mean_within_sys_yeo(n,1)=W;
%         mean_between_sys_yeo(n,1)=B;
%         %Connectivity between each of the 7 blocks
%         Ci=yeo_nodes;
%         nCi = unique(Ci);
%         M=subfcmat;
% 
%         for i = 1:length(nCi) % loop through communities
%             for j = 1:length(nCi)
%                Wi = Ci == nCi(i); % find index for within communitiy edges
%                Bi = Ci == nCi(j); % find index for between communitiy edges to specific community
% 
%                Wv_temp = M(Wi,Wi); % extract within communitiy edges
%                Bv_temp = M(Wi,Bi); % extract between communitiy edges to specific community
%                 
%                %calculate standard deviation of edges within blocks here
%                deviation_edge_weights_yeo(n,1)=deviation_edge_weights_yeo(n,1)+std(Wv_temp(Wv_temp~=0))^2; %Gu paper calculation, finished below
%                %Wv = [Wv_temp(logical(triu(ones(sum(Wi)),1)))'];  
%                Bv = [Bv_temp(:)'];
%                system_connectivity(i,j)=mean(Bv(Bv~=0));
%                %system_between(i,1)=mean(Bv);
%                %if i==j
%                %else
%             end
%         end
% 
%         %transpose these (or something) so they can be saved out on a subject basis.
%         system_connectivity;
%         system_connectivity_all(:,:,n)=system_connectivity;  %put this all into a matrix for everyone so that we can see the average
%         yeo_conn_vector = reshape(system_connectivity',[],1)';
%         yeo_connectivity(n,:)=yeo_conn_vector;
%         %finish the std deviation of edge weights calculation
%         deviation_edge_weights_yeo(n,1)=sqrt(deviation_edge_weights_yeo(n,1)/7);
%         %some net stats on the WSBM partition
%         %Participation coefficient average with Yeo partition
%         [Ppos Pneg]=participation_coef_sign(subfcmat, yeo_nodes);
%         %each row is a subject, with part coef for 359 nodes
%         sub_partcoef_pos_yeo(n,:)=mean(Ppos);
%         sub_partcoef_neg_yeo(n,:)=mean(Pneg);
%         %estimate modularity with this a-priori partition?
%          %% CHECK THAT THESE TWO Q's ARE THE SAME %%
%         [M Q]=modul_only(subfcmat, [], double(yeo_nodes), 'negative_asym');
%         %QFModul(double(yeo_nodes), subfcmat) %this gives the same answer as modul_only above when using on a pos-only matrix, but 
%         %we want the negative asymmteric weighting, so ignore for now.
%         % [M Q]=community_louvain(subfcmat, [], double(yeo_nodes), 'negative_asym'); %this optimizes modularity, which is not what we ant.
%         modul_yeo(n,1)=Q;
%     catch
%     fprintf('Cant read sub %s, skipped. \n', sub{1});
%     end
% end
% %save outfile
% outdir='/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/'
% outfile=dataset(subjlist, avgweight, system_segreg_yeo, mean_within_sys_yeo, mean_between_sys_yeo, deviation_edge_weights_yeo, sub_partcoef_pos_yeo, sub_partcoef_neg_yeo, modul_yeo, yeo_connectivity)
% export(outfile,'File',fullfile(outdir,'n546_test_sample_schaefer400_yeo_network_stats.csv'),'Delimiter',',')
% 
% %also save the mean system connectivity matrix from yeo
% mean_system_conn_mat=mean(system_connectivity_all,3)
% % header={'sys1', 'sys2', 'sys3','sys4','sys5','sys6','sys7'}
% % mean_system_conn_mat.Properties.VariableNames=header;
% save(fullfile(outdir, 'n546_mean_yeo_connectivity.mat'), 'mean_system_conn_mat')
% 

%% TRAINING SAMPLE
%This assumes you have already z-scored the FC matrices
listdir='/cbica/projects/spatial_topography/data/subjLists/release2/site16/parcellation'
z_avg_outdir='/cbica/projects/spatial_topography/data/imageData/fc_matrices/site16_training_sample/Schaefer400zavgNetworks'
outdir='/cbica/projects/spatial_topography/data/imageData/net_stats'
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/WSBM_v1.2'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/system_matrix_tools/'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/BCT/'))

subjlist=readtable(fullfile(listdir,'n670_filtered_runs_site16_postprocess.csv')) %the first two runs here are those that were input into WSBM before!
%read in the yeo partition
yeo_nodes=dlmread('/cbica/projects/spatial_topography/tools/parcellations/schaefer400/schaefer400x7CommunityAffiliation.1D.txt')

subjlist=subjlist.id;

%clear variables that might have overlapped with previous subject list
clear modul
clear avgweight
clear deviation_edge_weights_yeo
clear system_segreg_yeo
clear mean_within_sys_yeo
clear mean_between_sys_yeo
clear system_conn
clear sub_partcoef_pos_yeo
clear sub_partcoef_neg_yeo
clear modul_yeo
clear yeo_connectivity
clear system_connectivity_all

%preallocate variables
deviation_edge_weights_yeo=zeros(size(subjlist,1),1);
modul_yeo=zeros(size(subjlist,1),1);

%load in the FC matrix for each subject
for n=1:size(subjlist,1)
    sub=subjlist{n,:}
    try 
        %% LLLoad FC matrix
        fcfile=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
        %parcel 52 is already gone
        subfcmat = load(fcfile);
        %average edge weight for the subject
        avgweight(n,1)=mean(subfcmat(subfcmat~=0));
        %% Apply Yeo partition and calculate statistics with it
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
        system_connectivity_all(:,:,n)=system_connectivity;  %put this all into a matrix for everyone so that we can see the average
        yeo_conn_vector = reshape(system_connectivity',[],1)';
        yeo_connectivity(n,:)=yeo_conn_vector;
        %finish the std deviation of edge weights calculation
        deviation_edge_weights_yeo(n,1)=sqrt(deviation_edge_weights_yeo(n,1)/7);
        %some net stats on the WSBM partition
        %Participation coefficient average with Yeo partition
        [Ppos Pneg]=participation_coef_sign(subfcmat, yeo_nodes);
        %each row is a subject, with part coef for 359 nodes
        sub_partcoef_pos_yeo(n,:)=mean(Ppos);
        sub_partcoef_neg_yeo(n,:)=mean(Pneg);
        %estimate modularity with this a-priori partition?
         %% CHECK THAT THESE TWO Q's ARE THE SAME %%
        [M Q]=modul_only(subfcmat, [], double(yeo_nodes), 'negative_asym');
        %QFModul(double(yeo_nodes), subfcmat) %this gives the same answer as modul_only above when using on a pos-only matrix, but 
        %we want the negative asymmteric weighting, so ignore for now.
        % [M Q]=community_louvain(subfcmat, [], double(yeo_nodes), 'negative_asym'); %this optimizes modularity, which is not what we ant.
        modul_yeo(n,1)=Q;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end
%save outfile
outfile=dataset(subjlist, avgweight, system_segreg_yeo, mean_within_sys_yeo, mean_between_sys_yeo, deviation_edge_weights_yeo, sub_partcoef_pos_yeo, sub_partcoef_neg_yeo, modul_yeo, yeo_connectivity)
export(outfile,'File',fullfile(outdir,'n670_training_sample_schaefer400_yeo7_network_stats.csv'),'Delimiter',',')

%try to save outfile with correct header
header={'ID', 'avgweight', 'system_segreg_yeo', 'mean_within_sys_yeo', 'mean_between_sys_yeo', 'deviation_edge_weights_yeo', 'sub_partcoef_pos_yeo', 'sub_partcoef_neg_yeo', 'modul_yeo', 'sys1to1','sys1to2','sys1to3','sys1to4','sys1to5','sys1to6','sys1to7','sys2to1','sys2to2','sys2to3','sys2to4','sys2to5','sys2to6','sys2to7','sys3to1','sys3to2','sys3to3','sys3to4','sys3to5','sys3to6','sys3to7','sys4to1','sys4to2','sys4to3','sys4to4','sys4to5','sys4to6','sys4to7','sys5to1','sys5to2','sys5to3','sys5to4','sys5to5','sys5to6','sys5to7','sys6to1','sys6to2','sys6to3','sys6to4','sys6to5','sys6to6','sys6to7','sys7to1','sys7to2','sys7to3','sys7to4','sys7to5','sys7to6','sys7to7'}

outfile=table(subjlist, avgweight, system_segreg_yeo, mean_within_sys_yeo, mean_between_sys_yeo, deviation_edge_weights_yeo, sub_partcoef_pos_yeo, sub_partcoef_neg_yeo, modul_yeo, yeo_connectivity)
outfile2=splitvars(outfile, 'system_conn')
outfile2.Properties.VariableNames=header

save(fullfile(outdir, 'n670_training_sample_schaefer400_yeo7_network_stats.csv'), 'outfile2')
writetable(outfile2,fullfile(outdir,'n670_training_sample_schaefer400_yeo7_network_stats.csv'))

%outfile=table(char(unique_subjlist.id0), avgweight, modul, avgclustco_both, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg, mean_within_sys, mean_between_sys, system_conn)

%also save the mean system connectivity matrix from yeo
mean_system_conn_mat=mean(system_connectivity_all,3)
% header={'sys1', 'sys2', 'sys3','sys4','sys5','sys6','sys7'}
% mean_system_conn_mat.Properties.VariableNames=header;
save(fullfile(outdir, 'n670_mean_yeo7_connectivity.mat'), 'mean_system_conn_mat')
