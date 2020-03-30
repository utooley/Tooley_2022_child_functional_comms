#!/bin/bash
set -euo pipefail
#$ -j y
#$ -o /cbica/projects/spatial_topography/output/job_output/
#$ -l h_vmem=15.5G,s_vmem=15.3G
#$ -V
#$ -cwd

module unload matlab/R2014B
module load matlab/R2018A

code_dir='/cbica/projects/spatial_topography/code/net_stats'
matlab -nodisplay -r "cd ${code_dir}, run('net_stats_schaefer400_yeo7_both_samples.m'); exit"
