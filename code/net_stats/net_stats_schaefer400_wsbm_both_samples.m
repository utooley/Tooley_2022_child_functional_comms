%% Loop through samples
outdir='/cbica/projects/spatial_topography/data/imageData/net_stats'
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/WSBM_v1.2'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/system_matrix_tools/'))
addpath(genpath('/cbica/projects/spatial_topography/tools/matlab/BCT/'))
%% For each sample
samples={'site16','site14site20'}
full_names={'site16_training_sample','site14site20_test_sample'}
subjectlist={'parcellation/n670_filtered_runs_site16_postprocess.csv','n544_filtered_runs_site14site20_postprocess.csv'}
for s=1:length(samples)
    sample=samples{s}
    full_name=full_names{s}
    subjlist_file=subjectlist{s}
    listdir=strcat('/cbica/projects/spatial_topography/data/subjLists/release2/', sample,'/')
    z_avg_outdir=strcat('/cbica/projects/spatial_topography/data/imageData/fc_matrices/',full_name,'/Schaefer400zavgNetworks')

    %load subject list
    subjlist=readtable(fullfile(listdir,subjlist_file)) %the first two runs here are those that were input into gwMRF
    subjlist=subjlist.id
%% Run some partition quality metrics on the WSBM consensus partition and estimate network statistics
%read in the WSBM partition
load('/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/n670_training_sample_consensus_partitions_yeorelabeled.mat')
consensus_nodes=consensus_iter_mode_yeorelabeled;

%clear variables that might have overlapped with previous subject list
clear modul
clear avgweight
clear deviation_edge_weights_wsbm
clear system_segreg_wsbm
clear mean_within_sys_wsbm
clear mean_between_sys_wsbm
clear system_conn
clear sub_partcoef_pos_wsbm
clear sub_partcoef_neg_wsbm
clear modul_wsbm
clear wsbm_connectivity
clear system_connectivity_all
%preallocate variables
deviation_edge_weights_wsbm=zeros(size(subjlist,1),1);
modul_wsbm=zeros(size(subjlist,1),1);

%load in the FC matrix for each subject
for n=1:size(subjlist,1)
    sub=subjlist{n,:}
    try 
        %% Load FC matrix
        fcfile=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
        %parcel 52 is already gone
        subfcmat = load(fcfile);
        %average edge weight for the subject
%         avgweight(n,1)=mean(subfcmat(subfcmat~=0));
%         %% Apply WSBM partition and calculate statistics with it
%         [S, W, B] = segregation(subfcmat,consensus_nodes);
%         system_segreg_wsbm(n,1)=S;
%         mean_within_sys_wsbm(n,1)=W;
%         mean_between_sys_wsbm(n,1)=B;
%         %Connectivity between each of the 7 blocks
        Ci=consensus_nodes;
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

        %transpose these (or something) so they can be saved out on a subject basis.
        system_connectivity;
        system_connectivity_all(:,:,n)=system_connectivity;  %put this all into a matrix for everyone so that we can see the average
        wsbm_conn_vector = reshape(system_connectivity',[],1)';
        wsbm_connectivity(n,:)=wsbm_conn_vector;
        %finish the std deviation of edge weights calculation
%         deviation_edge_weights_wsbm(n,1)=sqrt(deviation_edge_weights_wsbm(n,1)/7);
%         %some net stats on the WSBM partition
%         %Participation coefficient average with wsbm partition
%         [Ppos Pneg]=participation_coef_sign(subfcmat, consensus_nodes);
%         %each row is a subject, with part coef for 359 nodes
%         sub_partcoef_pos_wsbm(n,:)=mean(Ppos);
%         sub_partcoef_neg_wsbm(n,:)=mean(Pneg);
%         %estimate modularity with this a-priori partition?
%          %% CHECK THAT THESE TWO Q's ARE THE SAME %%
%         [M Q]=modul_only(subfcmat, [], double(consensus_nodes), 'negative_asym');
%         %QFModul(double(wsbm_nodes), subfcmat) %this gives the same answer as modul_only above when using on a pos-only matrix, but 
%         %we want the negative asymmteric weighting, so ignore for now.
%         % [M Q]=community_louvain(subfcmat, [], double(wsbm_nodes), 'negative_asym'); %this optimizes modularity, which is not what we ant.
%         modul_wsbm(n,1)=Q;
    catch
    fprintf('Cant read sub %s, skipped. \n', sub);
    end
end

%save outfile
% outfile=dataset(subjlist, avgweight, system_segreg_wsbm, mean_within_sys_wsbm, mean_between_sys_wsbm, deviation_edge_weights_wsbm, sub_partcoef_pos_wsbm, sub_partcoef_neg_wsbm, modul_wsbm, wsbm_connectivity)
% export(outfile,'File',fullfile(outdir,strcat('n',string(size(subjlist,1)),'_',full_name,'_schaefer400_wsbm_network_stats.csv')),'Delimiter',',')
% 
% %try to save outfile with correct header
% header={'ID', 'avgweight', 'system_segreg_wsbm', 'mean_within_sys_wsbm', 'mean_between_sys_wsbm', 'deviation_edge_weights_wsbm', 'sub_partcoef_pos_wsbm', 'sub_partcoef_neg_wsbm', 'modul_wsbm', 'sys1to1','sys1to2','sys1to3','sys1to4','sys1to5','sys1to6','sys1to7','sys2to1','sys2to2','sys2to3','sys2to4','sys2to5','sys2to6','sys2to7','sys3to1','sys3to2','sys3to3','sys3to4','sys3to5','sys3to6','sys3to7','sys4to1','sys4to2','sys4to3','sys4to4','sys4to5','sys4to6','sys4to7','sys5to1','sys5to2','sys5to3','sys5to4','sys5to5','sys5to6','sys5to7','sys6to1','sys6to2','sys6to3','sys6to4','sys6to5','sys6to6','sys6to7','sys7to1','sys7to2','sys7to3','sys7to4','sys7to5','sys7to6','sys7to7'}
% 
% outfile=table(subjlist, avgweight, system_segreg_wsbm, mean_within_sys_wsbm, mean_between_sys_wsbm, deviation_edge_weights_wsbm, sub_partcoef_pos_wsbm, sub_partcoef_neg_wsbm, modul_wsbm, wsbm_connectivity)
% outfile2=splitvars(outfile, 'wsbm_connectivity')
% outfile2.Properties.VariableNames=header
% 
% save(fullfile(outdir, strcat('n',string(size(subjlist,1)),'_',full_name,'_schaefer400_wsbm_network_stats.csv')), 'outfile2')
% writetable(outfile2,fullfile(outdir,strcat('n',string(size(subjlist,1)),'_',full_name,'_schaefer400_wsbm_network_stats.csv')))
% 
%outfile=table(char(unique_subjlist.id0), avgweight, modul, avgclustco_both, num_comms_modul, part_coef_pos, part_coef_neg, system_segreg, mean_within_sys, mean_between_sys, system_conn)

%also save the mean system connectivity matrix from wsbm
mean_system_conn_mat=mean(system_connectivity_all,3)
% header={'sys1', 'sys2', 'sys3','sys4','sys5','sys6','sys7'}
% mean_system_conn_mat.Properties.VariableNames=header;
save(fullfile(outdir, strcat('n',string(size(subjlist,1)),'_',full_name,'_mean_wsbm_connectivity.mat')), 'mean_system_conn_mat')

end
