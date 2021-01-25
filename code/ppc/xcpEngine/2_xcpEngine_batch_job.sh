#!/bin/bash
set -euo pipefail
MACKEY_HOME=/data/picsl/mackey_group
project_dir=/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/
sublist_dir=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16
tools_dir=${MACKEY_HOME}/tools/singularity

FULL_COHORT=${sublist_dir}/n696_site16_cohort.csv
NJOBS=`wc -l < $FULL_COHORT`

if [[ ${NJOBS} == 0 ]]; then
    echo 'you dont have enough lines in your csv file'
    exit 0
fi
echo $NJOBS

#ALSO REMEMBER TO CHANGE YOUR .SGE_REQUEST FILE FOR MORE MEM, OTHERWISE THIS MAY NOT WORK
cat << EOF > xcpParallel.sh
#$ -j y
#$ -l h_vmem=75.6G,s_vmem=75.5G
#$ -o /data${project_dir}/output/job_output/xcpEngine
#$ -q himem.q,all.q,basic.q,gpu.q
#$ -t 1-${NJOBS}

# Adjust these so they work on your system
SNGL=/share/apps/singularity/2.5.1/bin/singularity
SIMG=${tools_dir}/xcpEngine-100219.simg
FULL_COHORT=${FULL_COHORT}

# Create a temp cohort file with 1 line
HEADER=\$(head -n 1 \$FULL_COHORT)
LINE_NUM=\$( expr \$SGE_TASK_ID + 1 )
LINE=\$(awk "NR==\$LINE_NUM" \$FULL_COHORT)
TEMP_COHORT=\${FULL_COHORT}.\${SGE_TASK_ID}.csv
echo \$HEADER > \$TEMP_COHORT
echo \$LINE >> \$TEMP_COHORT

\$SNGL run --cleanenv -B /data:/mnt \$SIMG \\
  -c /mnt\${TEMP_COHORT#/data} \\
  -d /mnt${project_dir}code/ppc/xcpEngine/fc_gsrwmcsf_scrub0.2mm_dropnonsteadystate.dsn \\
  -o /mnt/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols\\
  -r /mnt/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/fmriprep \\
  -i \$TMPDIR \\

EOF

qsub xcpParallel.sh
