#!/bin/bash
#$ -j y
#$ -l h_vmem=10.6G,s_vmem=10.5G
#$ -o /cbica/projects/spatial_topography/output/job_output/
set -euo pipefail

export SUBJECTS_DIR=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/freesurfer
data_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/xcpEngine_gsrwmcsf_scrub0.2mm_dropvols_marek
out_dir=/cbica/projects/spatial_topography/public_data/ABCD/bids_release2_site14site20/derivatives/surfaces

subject=${1}
fs_subject=${subject}
fwhm=4.8 #smoothing kernel, 2x voxel size, voxel size is 2.4
last_run=$(find ${data_dir}/${subject} -name "run-*" -type d| cut -d- -f3 |sort -n | tail -n1)
last_run=${last_run:1}
mkdir ${out_dir}/${subject}
mkdir ${out_dir}/${subject}/coreg
mkdir ${out_dir}/${subject}/surf

for sesh in $(seq 1 $last_run); #create the registration from the bold subject space (refvolume below) to freesurfer surfaces
do
if [ -f ${out_dir}/${subject}/surf/lh.fs6_sm${fwhm}_${subject}_run-0${sesh}_task-rest_uncensor_trunc_reshape.nii.gz ]; then
  echo $sub $sesh 'done already'
  # if it exists already, skip it
else
  if [ -f ${data_dir}/${subject}/run-0${sesh}/regress/${subject}_run-0${sesh}_uncensored.nii.gz ]; then
    img=${data_dir}/${subject}/run-0${sesh}/regress/${subject}_run-0${sesh}_uncensored
  else
    img=${data_dir}/${subject}/run-0${sesh}/regress/${subject}_run-0${sesh}_residualised
  fi
  #truncate functional runs to 370 volumes
  fslroi ${img}.nii.gz ${data_dir}/${subject}/run-0${sesh}/regress/${subject}_run-0${sesh}_uncensored_truncated.nii.gz 0 370
  echo session is $sesh
  bbregister --s ${fs_subject} \
      --mov ${data_dir}/${subject}/run-0${sesh}/prestats/*referenceVolumeBrain.nii.gz \
      --reg ${out_dir}/${subject}/coreg/${subject}_run-0${sesh}_fs_epi2struct.dat \
      --init-fsl --bold
#to check run
#tkregisterfv --mov ${data_dir}/${subject}/run-0${sesh}/prestats/*referenceVolumeBrain.nii.gz --reg ${out_dir}/${subject}/coreg/${subject}_run-0${sesh}_fs_epi2struct.dat --surfs &

for hem in lh rh
do
  echo hem is $hem
  img=${data_dir}/${subject}/run-0${sesh}/regress/${subject}_run-0${sesh}_uncensored_truncated.nii.gz #image in subject/input space from xcpEngine
  mri_vol2surf --mov ${img} \
              --reg ${out_dir}/${subject}/coreg/${subject}_run-0${sesh}_fs_epi2struct.dat \
              --hemi ${hem} \
              --o ${out_dir}/${subject}/surf/${hem}.${subject}_run-0${sesh}_task-rest_uncensor_trunc_reshape.nii.gz \
              --projfrac 0.5 \
              --interp trilinear \
              --reshape \
              #Yeo script uses --reshape and no smoothing with surf-fwhm here, smooths later in surf2surf
done

#this looks good, check this by running freeview --recon ${sub} and overlaying the functional image on the lh.inflated of the subject.

for hem in lh rh
do
# you can also use --trgsubject ico and --icoorder instead of trgsubject, Yeo code uses --cortex
mri_surf2surf --srcsubject ${fs_subject} \
            --srcsurfval ${out_dir}/${subject}/surf/${hem}.${subject}_run-0${sesh}_task-rest_uncensor_trunc_reshape.nii.gz \
            --trgsubject fsaverage6 \
            --trgsurfval ${out_dir}/${subject}/surf/${hem}.fs6_sm${fwhm}_${subject}_run-0${sesh}_task-rest_uncensor_trunc_reshape.nii.gz \
            --hemi ${hem} \
            --fwhm-trg ${fwhm} \
            --cortex \
            --reshape \
#Check this by running freeview -f $SUBJECTS_DIR/fsaverage6/surf/lh.inflated:overlay={hem}_fs6_sm${fwhm}_${subject}_run-0${sesh}_task-rest_residualised.nii.gz

#This is not needed if you want to stay in fsaverage6 space.
# mri_surf2surf --srcsubject fsaverage6 \
#             --srcsurfval ${out_dir}/${subject}/surf/${hem}_fs6_sm${fwhm}_${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
#             --trgsubject fsaverage5 \
#             --trgsurfval ${out_dir}/${subject}/surf/${hem}_fs5_sm${fwhm}_${subject}_run-0${sesh}_task-rest_residualised_fromfs6.nii.gz \
#             --hemi ${hem} \
#             --nsmooth-in 1

done
fi
done
