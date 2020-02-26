#!/bin/bash
unset PYTHONPATH;
MACKEY_HOME=/data/picsl/mackey_group
project_dir=/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/
sublist_dir=/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists
tools_dir=${MACKEY_HOME}/tools/singularity


SNGL=/share/apps/singularity/2.5.1/bin/singularity
SIMG=${tools_dir}/xcpEngine.simg
FULL_COHORT=/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/n2_mult_runs_test.csv

$SNGL run --cleanenv -B /data:/mnt $SIMG \
  -c /mnt${FULL_COHORT} \
  -d /mnt${project_dir}/code/ppc/xcpEngine/fc-36p_scrub_edited_dropvols.dsn \
  -o /mnt/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_censored/ \
  -i $TMPDIR \
  -r /mnt/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/fmriprep/


#qsub -q himem.q -j y -o /data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output -l h_vmem=29.1G,s_vmem=29.0G /data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/ppc/xcpEngine/xcpEngine_cmd.sh
