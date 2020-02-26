function [temp_weights,temp_edges,temp_evidence,temp_labels] = wsbm_function(sub,num_comms, wsbm_dir, z_avg_outdir)
%WSBM FUNCTION FOR PARFOR LOOP
        outfile=fullfile(wsbm_dir,strcat(sub,'_k',num2str(num_comms),'_wsbm.mat'))
        if (exist(outfile)==2) %if it's already written don't do it again
            fprintf('Sub %s already exists. \n', sub);
            load(outfile)
            temp_weights=Model.Data.logHw;
            temp_edges=Model.Data.logHe;
            temp_evidence=Model.Para.LogEvidence;
            temp_labels=Labels;
        else
        file=fullfile(z_avg_outdir,strcat(sub,'_avg_Schaefer400x7_znetwork.txt'));
        %try
        subfcmat=load(file);
        [Labels Model]=wsbm(subfcmat, num_comms,'E_Distr','None', 'W_Distr', 'Normal', 'numTrials', 30, 'verbosity', 1, 'alpha', 0, 'parallel', 0);
        %logHw - Additive Log-likelihood constant for W_Distr
        temp_weights=Model.Data.logHw;
        %logHe - Additive Log-likelihood constant for E_Distr
        temp_edges=Model.Data.logHe;
        %LogEvidence - Marginal Log-likelihood (aka Log-Evidence), a model selection criterion 
        temp_evidence=Model.Para.LogEvidence;
        %a row for each subject's labels for the k=num_comms partition
        temp_labels=Labels;
        %Look at how many unique labels were output, since Rick says these might be
        %different
        %save the model
        outfile=fullfile(wsbm_dir,strcat(sub,'_k',num2str(num_comms),'_wsbm.mat'))
        save(outfile, 'Model', 'Labels')
        end
end

