URSULA_PROJ=/data/jux/mackey_group/Ursula/projects/in_progress
ERROR_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output
SUBLIST_DIR=${URSULA_PROJ}/spatial_topography_parcellations_ABCD/data/subjLists
SUB_FILE=${SUBLIST_DIR}/n27_one_run_only_21019.txt #make sure this file has \n line endings
SCRIPTS_DIR=${URSULA_PROJ}/spatial_topography_parcellations_ABCD/code/ppc/fmriprep
FMRIPREP_DIR=/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/fmriprep

#NEED TO SEE WHICH SITE IT IS?

## REMOVE THE SUBJECTS WHO WEREN'T IN THE ORIGINAL n26 batch, rerun with 1.3.0 FMRIPREP
for sub in `find . -maxdepth 1 -mindepth 1 -type d -name "sub-*" | sed -e 's|.*/||'`:
do
echo $sub
if grep -Fxq "$sub" ${SUB_FILE}; then
  echo $sub
  echo 'it exists'
else
  echo $sub
  echo 'it doesnt exist'
  rm $sub -R
  # sleep 1
  # qsub -q all.q,basic.q,himem.q -j y -o $ERROR_DIR -l h_vmem=25.1G,s_vmem=25.0G ${SCRIPTS_DIR}/fmriprep_cmd2.sh ${sub}
fi
done

## COPY FROM /data/picsl/ to /data/jux/
output_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/fmriprep
fs_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/freesurfer
URSULA_PROJ=/data/jux/mackey_group/Ursula/projects/in_progress
ERROR_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output/fmriprep
SUBLIST_DIR=${URSULA_PROJ}/spatial_topography_parcellations_ABCD/data/subjLists/release2
SCRIPTS_DIR=${URSULA_PROJ}/spatial_topography_parcellations_ABCD/code/ppc/fmriprep
SUB_FILE=${SUBLIST_DIR}/site14site20/n611_release2_site14site20_0.2mm.txt
for sub in `cat ${SUB_FILE}`
#for sub in `find . -maxdepth 1 -mindepth 1 -type d -name "sub-*" | sed -e 's|.*/||'`:
do
  #sub=sub-${sub} #check what your subject list is!
#echo $sub
if [ -d ${output_dir}/$sub/func ]; then
  echo $sub
  # echo 'fmriprep finished'
  #sleep 0.1
else
  # echo $sub
  # echo 'func folder doesnt exist'
  # # echo ${output_dir}/${sub}
  # rm ${output_dir}/${sub} -R
  # echo deleted fmriprep folder
  # echo ${fs_dir}/${sub}
  # rm ${fs_dir}/${sub} -R
  # echo deleted freesurfer folder
  sleep 0.001
  # qsub -q all.q,basic.q,himem.q,gpu.q -j y -o $ERROR_DIR -l h_vmem=19.1G,s_vmem=19.0G ${SCRIPTS_DIR}/fmriprep_cmd2.sh ${sub}
fi
done

##RERUN THE REMAINDER with 1.2.6-1 FMRIPREP
for sub in `find . -maxdepth 1 -mindepth 1 -type d -name "sub-*" | sed -e 's|.*/sub-||'`
do
echo $sub
if [ -d ${FMRIPREP_DIR}/sub-$sub ];then
  echo 'FMRIPREP already ran'
else
  echo 'Submitting FMRIPREP again'
sleep 1
qsub -q all.q,basic.q,himem.q,gpu.q -j y -o $ERROR_DIR -l h_vmem=25.1G,s_vmem=25.0G ${SCRIPTS_DIR}/fmriprep_cmd2.sh ${sub}
fi
done
