#$ -j y
#$ -l h_vmem=75.6G,s_vmem=75.5G
#$ -o /data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD//output/job_output/xcpEngine
#$ -q himem.q,all.q,basic.q,gpu.q
#$ -t 1-20

# Adjust these so they work on your system
SNGL=/share/apps/singularity/2.5.1/bin/singularity
SIMG=/data/picsl/mackey_group/tools/singularity/xcpEngine-100219.simg
FULL_COHORT=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site14site20/rerunning_site14site20_again_secondtime_cohort.csv

# Create a temp cohort file with 1 line
HEADER=$(head -n 1 $FULL_COHORT)
LINE_NUM=$( expr $SGE_TASK_ID + 1 )
LINE=$(awk "NR==$LINE_NUM" $FULL_COHORT)
TEMP_COHORT=${FULL_COHORT}.${SGE_TASK_ID}.csv
echo $HEADER > $TEMP_COHORT
echo $LINE >> $TEMP_COHORT

$SNGL run --cleanenv -B /data:/mnt $SIMG \
  -c /mnt${TEMP_COHORT#/data} \
  -d /mnt/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/ppc/xcpEngine/fc_gsrwmcsf_scrub0.2mm_dropnonsteadystate.dsn \
  -o /mnt/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek\
  -r /mnt/jux/mackey_group/public_data/ABCD/bids_release2_site14site20/derivatives/fmriprep \
  -i $TMPDIR \

