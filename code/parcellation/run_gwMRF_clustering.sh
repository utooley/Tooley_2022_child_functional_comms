#!/bin/bash
set -euo pipefail
#$ -j y
#$ -o /cbica/home/tooleyu/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output/
#$ -l h_vmem=500.5G,s_vmem=500.3G
#$ -V
#$ -cwd

module unload matlab/R2014B
module load matlab/R2018A
# CHANGE SUBJECT INDICES FOR NUM OF SUBJECTS!
CBIG_code_dir='/cbica/home/tooleyu/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code'
filepaths='/cbica/home/tooleyu/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/parcellation/n131_site16_surfs_tworunsonly_nonan_subjects_CUBIC.csv'
outdir='/cbica/home/tooleyu/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/gwMRF_n131_no_nansubjs_site16_normalized_original_gamma/'
matlab -nodisplay -r "cd ${CBIG_code_dir};CBIG_gwMRF_build_data_and_perform_clustering('${filepaths}','${outdir}',1,130,200,200,150000,100,85,500000000,15); exit"

#input_fullpaths,output_path,start_idx,end_idx,num_left_cluster,num_right_cluster,smoothcost,num_iterations,num_runs,start_gamma,exponential)

#for just the clustering
#matlab -nodisplay -r "cd ${CBIG_code_dir};CBIG_gwMRF_graph_cut_clustering_split_newkappa_prod('${outdir}/mult_mat/lh_mult_matrix.mat',prams,'lh'); exit"
#prams is what I can't figure out how to input.

#for running in matlab on CUBIC
#cd /gpfs/fs001/cbica/home/tooleyu/tools/CBIG/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Code
#CBIG_gwMRF_build_data_and_perform_clustering(filepaths,outdir,1,130,200,200,150000,100,100,500000000,15)


#want to change these params to be suitable for actually running this.
#for 400-parcel parcellation, gradient weight c (smoothness) set to 150000
#Decay parameter k (exponential) set to 15
# Initial value of tau (gamma here), the smoothness parameter for spatial roundness, 5x 10^8
#num of iterations--not given in paper, default is 100 in the params code.
#number of runs-500 for comparison with other parcellations, 1000 for 400-parcels and above resolution
#CBIG_gwMRF_build_data_and_perform_clustering('/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/parcellation/example_input_fullpaths.csv','/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/gwMRF_cluster_test',1,2,50,50,5000,7,2,50000000,15)
