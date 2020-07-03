resample_dir=${CBIG_CODE_DIR}/data/templates/surface/standard_mesh_atlases_20170508/resample_fsaverage
task_maps_dir=/cbica/projects/spatial_topography/data/imageData/task_maps

for file in `find ${task_maps_dir} -iname "*_cohen_c1.dscalar.nii"  -exec basename {} \;`
do

cifti_in=${file}
#cifti_in=MID_antic_large_loss_vs_neutral_dat_cohen_c1.dscalar.nii
metric_out_left=${task_maps_dir}/lh.${cifti_in%.dscalar.nii}.func.gii
metric_out_right=${task_maps_dir}/rh.${cifti_in%.dscalar.nii}.func.gii

current_sphere_left=${resample_dir}/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii
current_sphere_right=${resample_dir}/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii
new_sphere_left=${resample_dir}/fsaverage6_std_sphere.L.41k_fsavg_L.surf.gii
new_sphere_right=${resample_dir}/fsaverage6_std_sphere.R.41k_fsavg_R.surf.gii
final_out_left=${task_maps_dir}/ABCD.${cifti_in%dat_cohen_c1.dscalar.nii}.lh.41k_fsavg_L.func.gii
final_out_right=${task_maps_dir}/ABCD.${cifti_in%dat_cohen_c1.dscalar.nii}.rh.41k_fsavg_R.func.gii
current_area_left=${resample_dir}/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii
current_area_right=${resample_dir}/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii
new_area_left=${resample_dir}/fsaverage6.L.midthickness_va_avg.41k_fsavg_L.shape.gii
new_area_right=${resample_dir}/fsaverage6.R.midthickness_va_avg.41k_fsavg_R.shape.gii


#wb_command -cifti-separate ${task_maps_dir}/${cifti_in} COLUMN -metric CORTEX_LEFT ${metric_out_left} -metric CORTEX_RIGHT ${metric_out_right}

#wb_command -metric-resample ${metric_out_left} ${current_sphere_left} ${new_sphere_left} ADAP_BARY_AREA ${final_out_left} -area-metrics ${current_area_left} ${new_area_left}

wb_command -metric-resample ${metric_out_right} ${current_sphere_right} ${new_sphere_right} ADAP_BARY_AREA ${final_out_right} -area-metrics ${current_area_right} ${new_area_right}

done
