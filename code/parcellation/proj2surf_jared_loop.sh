export SUBJECTS_DIR=/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/freesurfer
data_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek
out_dir=/data/jux/mackey_group/public_data/ABCD/bids_release2_site16/derivatives/surfaces
scripts_dir=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/code/parcellation

#make a subject list with only the good runs (motion) and good subjects (coreg)
#for now, just trying with everyone who was run through the original xcp pipeline with 36 p
subjlist=/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/site16/n670_filtered_runs_site16_postprocess.csv

for subject in `cat ${subjlist} | cut -d, -f2 | uniq`
do
  echo $subject
  qsub -j y -l h_vmem=10.6G,s_vmem=10.5G ${scripts_dir}/proj2surf_jaredmethod.sh ${subject}
done
