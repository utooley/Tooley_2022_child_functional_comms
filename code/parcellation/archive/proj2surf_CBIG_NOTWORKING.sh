#!/bin/bash
#$ -j y
#$ -l h_vmem=10.6G,s_vmem=10.5G
#$ -o /data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output/
#$ -q himem.q,all.q,basic.q,gpu.q

subject=${1}

$CBIG_CODE_DIR/stable_projects/preprocessing/CBIG_fMRI_Preproc2016/CBIG_preproc_native2fsaverage.csh \
-s ${subject} \
-d /data/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_36p_gsr_multruns \
-anat_s sub-${subject} -anat_d /data/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/freesurfer \
-bld 'run-01 run-02 run-03 run-04' \
-proj fsaverage6 \
-down fsaverage5 \
-sm 5 \
-BOLD_stem f \
-REG_stem mni152.register.dat \
#this registration is hard-coded in now.
#the BOLD_stem is unused, I coded it out.
#mni registration at /share/apps/freesurfer/6.0.0/average/mni152.register.dat

#this is the base projection command, checked that the script gave same results
#mri_vol2surf --mov ${subject}_img_sm4Std.nii.gz --mni152reg --hemi  lh --projfrac 0.5 --trgsubject fsaverage6 --o lh.test.nii.gz --reshape --interp trilinear
