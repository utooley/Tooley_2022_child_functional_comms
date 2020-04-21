#!/bin/bash
subj_list=/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/n129_one_site_0.2mm_nobogus
scripts_dir=/data/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/parcellation/

for sub in `cat $subj_list`
do
  #sub=${sub:4} #take off the sub- part
  #echo ${sub}
  qsub ${scripts_dir}/proj2surf_CBIG.sh ${sub}
done
