##!/usr/bin/env bash
set -euo pipefail

if [ $# -eq 0 ]; then
echo "Usage: mk_full_paths_to_surfs.sh <subject_list> <output_dir>

This is DEPRECATED--no longer needed--filepaths to surfaces are now made in the R script
filter_subjlist_postxcp_deciding_pipeline.R"
exit
fi
# Set paths
data_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/surfaces
output_dir=${2}
subjlist=${1}
count=`wc -l < ${1}`

# Create a file to store example full paths
output_file="${output_dir}/n${count}_surfs_fullpaths_2runsonly.csv"
touch ${output_file}

for sub in `cat ${subjlist}`
do
  echo ${sub}
  echo -e ""  >>   ${output_file} #make a new line for the next subject
  for run in `find ${data_dir}/${sub}/surf/ -name "lh.fs6_sm4.8_${sub}_run-*" -type f| cut -d_ -f8 | head -n 2` #find only the first two runs
  do #echo all runs into one line
  echo -n "${data_dir}/${sub}/surf/lh.fs6_sm4.8_${sub}_${run}_task-rest_residualised_reshape.nii.gz " >>   ${output_file}
done
done


for a,b in $(cat ${subjlist}|cut -d "," -f 1,2);
do
echo $a
echo $b
done
