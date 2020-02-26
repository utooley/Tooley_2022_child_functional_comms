#delete half finished xcpEngine jobs because ran out of memory
set -euo pipefail
MACKEY_HOME=/data/picsl/mackey_group
project_dir=/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/
sublist_dir=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20
output_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek
FMRIPREP_DIR=/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/fmriprep
COHORT_FILE=${sublist_dir}/rerunning_site14site20_again_secondtime_cohort.csv
subjlist=${sublist_dir}/n597_site14site20_cohort.csv

echo id0,id1,img >> ${COHORT_FILE} #add header
for sub in `cat ${subjlist} | cut -d, -f1 | uniq`
do
  for run in {1..5};
    do
  #sub=sub-${sub} #check what your subject list is!
  #echo $sub
  #echo ${output_dir}/$sub/run-0${run}/qcfc/${sub}_run-0${run}_dvars-mean.txt
if [ -f ${output_dir}/${sub}/run-0${run}/qcfc/${sub}_run-0${run}_dvars-mean.txt ]; then
  echo $sub
  # echo 'fmriprep finished'
  #sleep 0.1
else
  echo $sub not finished
  #delete their run folder
  rm ${output_dir}/${sub}/run-0${run} -Rf
  #pull their input file and put it into a new cohort file
  cd ${FMRIPREP_DIR}
  find . -type f | grep "${sub}_task-rest_run-0${run}_space-T1w_desc-preproc_bold.nii.gz\$" | while read fname; do
  tmp=$(echo "$fname" | awk -F '_' '{print $3}') #this parses on underscores and pulls 'run-01'
  echo $sub,$tmp,${fname:2} >> ${COHORT_FILE}
done
  sleep 0.001
fi
done
done
