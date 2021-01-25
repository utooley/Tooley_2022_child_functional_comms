#!/bin/sh
SUBLIST_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16
COHORT_FILE=${SUBLIST_DIR}/n696_site16_cohort.csv

# THIS WORKS to pull all those that have completed fmriprep, or
#all those in a given list

echo id0,id1,img >> ${COHORT_FILE} #add header
for sub in `cat ${SUBLIST_DIR}/n696_site16_boldandt1_fmriprep.txt`
do
    echo $sub
    #find . -type f | grep "${sub}_*"
    #find . -iregex '.*preproc_bold.nii.gz$\|.*cgi$' -exec grep -il '${sub}_task' '{}' ';'
    #find . -type f | grep "sub-${sub}_task-rest_run-[0-9][0-9]_space-T1w_desc-preproc_bold.nii.gz\$"| while read fname; do
    find . -type f | grep "${sub}_task-rest_run-.._space-T1w_desc-preproc_bold.nii.gz\$" | while read fname; do
    tmp=$(echo "$fname" | awk -F '_' '{print $3}') #this parses on underscores and pulls 'run-01'
    echo $sub,$tmp,${fname:2} >> ${COHORT_FILE}
done;
done;
