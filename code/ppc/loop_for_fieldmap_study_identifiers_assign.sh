# #/bin/bash

cd /data/jux/mackey_group/public_data/ABCD/release2_site20/niftis/
SCRIPTS_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/ppc

for sub in `ls`
do
  sub=${sub:4} #take off the sub- part
  echo ${sub}
  python ${SCRIPTS_DIR}/assign_fieldmaps_to_IntendedFor_field.py ${sub}
done


### loop for fixing study identifiers when heudiconv breaks###
################################
# SCRIPTS_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/ppc
# #list of subjs to fix
# declare -a arr=(NDARINVAU2AGPUR NDARINVC0W9C4R6 NDARINVM7FJXJRR NDARINVNER12NXX NDARINVNJRMP62L NDARINVPR3T5MAK NDARINVV4KUL4N2 NDARINVDMMAKV5Y NDARINVV7NEVHLK NDARINVZLP46GRP)
# cd /data/picsl/mackey_group/public_data/ABCD/release2_site14site16_dicoms/
#
# ## now loop through the above array
# for i in "${arr[@]}"
# do
#    echo "$i"
#    cd /data/jux/mackey_group/public_data/ABCD/release2_site14site16_dicoms/sub-${i}
#    python ${SCRIPTS_DIR}/fix_conflicting_study_identifiers.py
# done
