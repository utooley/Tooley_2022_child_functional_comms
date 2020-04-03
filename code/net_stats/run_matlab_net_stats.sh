#!/bin/bash
set -euo pipefail
#$ -j y
#$ -o /cbica/projects/spatial_topography/output/job_output/
#$ -l h_vmem=200.5G,s_vmem=200.3G
#$ -V
#$ -cwd

module unload matlab/R2014B
module load matlab/R2018A
subject=${1}

code_dir='/cbica/projects/spatial_topography/code/net_stats'
matlab -nodisplay -r "cd ${code_dir}, net_stats_surfaces_yeodev_subject_script('${subject}'); exit"
