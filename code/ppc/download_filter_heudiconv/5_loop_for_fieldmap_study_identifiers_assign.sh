# #/bin/bash

cd /data/jux/mackey_group/public_data/ABCD/release2_site16/niftis/
SCRIPTS_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/ppc

for sub in `ls` #bad practice, don't do this in the future
do
  sub=${sub:4} #take off the sub- part
  echo ${sub}
  python ${SCRIPTS_DIR}/5_assign_fieldmaps_to_IntendedFor_field.py ${sub}
done
