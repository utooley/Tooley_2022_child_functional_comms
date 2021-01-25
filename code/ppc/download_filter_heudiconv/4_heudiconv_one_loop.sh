#/bin/bash
set -euo pipefail

# obtain scan and session labels
scans=/data/jux/mackey_group/public_data/ABCD/release2_site16_dicoms/*
SCRIPTS_DIR=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/ppc
#for a custom set of subjects
#declare -a arr=(NDARINVAU2AGPUR NDARINVC0W9C4R6 NDARINVM7FJXJRR NDARINVNER12NXX NDARINVNJRMP62L NDARINVPR3T5MAK NDARINVV4KUL4N2 NDARINVDMMAKV5Y NDARINVV7NEVHLK NDARINVZLP46GRP)

for sc in $scans; #if running everyone in the folder
for subID in "${arr[@]}";
	do subID=$(echo $sc|cut -d'/' -f8 |cut -c 5-); 	#if running everyone in the folder
	do echo ${subID}
  qsub ${SCRIPTS_DIR}/4_heudiconv_cmd.sh ${subID}

done
