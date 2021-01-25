#$ -j y
#$ -l h_vmem=19.1G,s_vmem=19.0G
#$ -o /data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output/fmriprep/
#$ -q himem.q,all.q,basic.q,gpu.q

unset PYTHONPATH;
sub=${1}
tools_dir=/data/picsl/mackey_group/tools/singularity
BIDS_folder=/data/jux/mackey_group/public_data/ABCD/bids_release2_site16
output_dir=${BIDS_folder}/derivatives/
echo ${sub}
echo ${BIDS_folder}
user=`whoami` #so that templateflow can go into the home dir of whoever is running.

export SINGULARITYENV_TEMPLATEFLOW_HOME=/home/${user}/templateflow
singularity run --cleanenv -B /home/${user}/templateflow:/home/${user}/templateflow,${BIDS_folder}:/mnt ${tools_dir}/fmriprep-1-4-1rc5.simg \
/mnt/ /mnt/derivatives/ participant \
--participant-label ${sub} \
--fs-license-file $HOME/license.txt \
--output-spaces MNI152NLin6Asym T1w fsaverage:den-10k \
--ignore slicetiming \
-w $TMPDIR \
