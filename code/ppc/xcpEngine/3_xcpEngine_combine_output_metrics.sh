#command for running the whole cohort as a group to collate output
#!/bin/bash
MACKEY_HOME=/data/picsl/mackey_group
project_dir=/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/
sublist_dir=/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists
tools_dir=${MACKEY_HOME}/tools/singularity

#run this command to start the xco container interaactively
singularity shell --cleanenv -B /data:/mnt ${tools_dir}/xcpEngine-100219.simg \

tools_dir=${MACKEY_HOME}/tools/singularity

#then run this command
outputdir=/mnt/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols

${XCPEDIR}/utils/combineOutput \
   -p $outputdir  \
   -f "*quality.csv" \
   -o XCP_QAVARS.csv \

${XCPEDIR}/utils/combineOutput \
    -p $outputdir  \
    -f "*audit.csv" \
    -o XCP_AUDIT.csv \
