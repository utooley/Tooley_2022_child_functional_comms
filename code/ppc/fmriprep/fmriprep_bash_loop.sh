#!/bin/sh

#Loop script for FMRIPREP
#unset PYTHONPATH;
URSULA_PROJ=/data/jux/mackey_group/Ursula/projects/in_progress
ERROR_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output
SUBLIST_DIR=${URSULA_PROJ}/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16
SUB_FILE=${SUBLIST_DIR}/n698_release2_site16_0.2mm_bids.txt
SCRIPTS_DIR=${URSULA_PROJ}/spatial_topography_parcellations_ABCD/code/ppc/fmriprep
output_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site16/derivatives/fmriprep
user=`whoami`

mkdir home/${user}/templateflow
#for sub in `find . -maxdepth 1 -mindepth 1 -type d -name "sub-*" | sed -e 's|.*/sub-||'`:
for sub in `cat ${SUB_FILE}`
do
echo $sub
sleep .5
# if [ -d ${output_dir}/sub-$sub ]; then
#   echo $sub
#   echo 'it exists'
# else
#   echo $sub
#   echo 'it doesnt exist'
  qsub -q all.q,basic.q,himem.q -o $ERROR_DIR ${SCRIPTS_DIR}/fmriprep_cmd.sh ${sub}
# fi
done
