export SUBJECTS_DIR=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/freesurfer
data_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols
out_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site16/derivatives/surfaces
scripts_dir=/cbica/projects/spatial_topography/code/proj2surf

subjlist=/cbica/projects/spatial_topography/data/subjLists/release2/site16/n670_filtered_runs_site16_postprocess.csv

for subject in `cat ${subjlist} | cut -d, -f2 | uniq`
do
  echo $subject
  qsub -j y -l h_vmem=10.6G,s_vmem=10.5G ${scripts_dir}/proj2surf.sh ${subject}
done
