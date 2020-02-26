#Before was trying proj2surf_jaredmethod.sh, tried these, they don't work!

###################
### NOT WORKING ###
####################
#alternative to the above, try doing both in one step and projecting directly to fsaverage6, then smoothing and downsampling to fsaverage5
# THIS DOES NOT WORK! BECAUSE THE REGISTRATION FILE DOES NOT MATCH THE TARGET SUBJECT
img=${data_dir}/${subject}/run-0${sesh}/regress/${subject}_run-0${sesh}_residualised.nii.gz #image in subject/input space from xcpEngine
mri_vol2surf --mov ${img} \
            --reg ${out_dir}/${subject}/coreg/${subject}_run-0${sesh}_fs_epi2struct.dat \
            --hemi ${hem} \
            --o ${out_dir}/${subject}/surf/${hem}.fsaverage5.sm${fwhm}.${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
            --projfrac 0.5 \
            --interp trilinear \
            --reshape --surf-fwhm ${fwhm} \
            --trgsubject fsaverage5 \ #or try going to fsaverage6 and then downsampling

####################
####### An alternative after vol2surf ####
#######################
#OR project to 6, smooth there in a separate command, and then downsample. This gives slightly different (2nd decimal level) results
#than doing the first two in one step.
mri_surf2surf --srcsubject ${fs_subject} \
            --srcsurfval ${out_dir}/${subject}/surf/${hem}.${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
            --trgsubject fsaverage6 \
            --trgsurfval ${out_dir}/${subject}/surf/${hem}.fs6first.${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
            --hemi ${hem} \
            --reshape

mri_surf2surf --s fsaverage6 \
            --sval ${out_dir}/${subject}/surf/${hem}.fs6first.${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
            --tval ${out_dir}/${subject}/surf/${hem}.fs6.smafter${fwhm}.${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
            --hemi ${hem} \
            --reshape \
            --fwhm-trg ${fwhm} \
            --cortex \

mri_surf2surf --srcsubject fsaverage6 \
            --sval ${out_dir}/${subject}/surf/${hem}.fs6.smafter${fwhm}.${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
            --trgsubject  fsaverage5 \
            --tval ${out_dir}/${subject}/surf/${hem}.fs5.smafter${fwhm}.${subject}_run-0${sesh}_task-rest_residualised.nii.gz \
            --nsmooth-in 1 \
            --hemi ${hem} \
            --reshape \
