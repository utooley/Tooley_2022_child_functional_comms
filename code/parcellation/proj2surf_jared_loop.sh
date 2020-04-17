export SUBJECTS_DIR=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/freesurfer
data_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek
out_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/surfaces
scripts_dir=/cbica/projects/spatial_topography/code/parcellation

#make a subject list with only the good runs (motion) and good subjects (coreg)
#for now, just trying with everyone who was run through the original xcp pipeline with 36 p
subjlist=/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/n546_filtered_runs_site14site20_postprocess.csv

for subject in `cat ${subjlist} | cut -d, -f2 | uniq`
do
  echo $subject
  qsub -j y -l h_vmem=10.6G,s_vmem=10.5G ${scripts_dir}/proj2surf_jaredmethod_site14site20_CUBIC.sh ${subject}
done
