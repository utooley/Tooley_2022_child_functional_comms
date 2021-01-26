#Following instructions here:https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP_5_8.pdf
resample_dir=${CBIG_CODE_DIR}/data/templates/surface/standard_mesh_atlases_20170508/resample_fsaverage

#LEFT
current_sphere=${resample_dir}/fsaverage6_std_sphere.L.41k_fsavg_L.surf.gii
new_sphere=${resample_dir}/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii
current_area=${resample_dir}/fsaverage6.L.midthickness_va_avg.41k_fsavg_L.shape.gii
new_area=${resample_dir}/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii
infile=lh.wsbm.fsaverage6.annot
outfile=wsbmpartition.lh.32k_fs_LR.label.gii

echo $infile

#convert to gifti
mris_convert --annot ${infile} ${SUBJECTS_DIR}/fsaverage6/surf/lh.white ${infile}.gii

#map label data from fsaverage group to fs_LR
wb_command -label-resample ${infile}.gii ${current_sphere} ${new_sphere} ADAP_BARY_AREA ${outfile} -area-metrics ${current_area} ${new_area}

#RIGHT
current_sphere=${resample_dir}/fsaverage6_std_sphere.R.41k_fsavg_R.surf.gii
new_sphere=${resample_dir}/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii
current_area=${resample_dir}/fsaverage6.R.midthickness_va_avg.41k_fsavg_R.shape.gii
new_area=${resample_dir}/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii
infile=rh.wsbm.fsaverage6.annot
outfile=wsbmpartition.rh.32k_fs_LR.label.gii

echo $infile

#convert to gifti
mris_convert --annot ${infile} ${SUBJECTS_DIR}/fsaverage6/surf/rh.white ${infile}.gii

#map label data from fsaverage group to fs_LR
wb_command -label-resample ${infile}.gii ${current_sphere} ${new_sphere} ADAP_BARY_AREA ${outfile} -area-metrics ${current_area} ${new_area}
