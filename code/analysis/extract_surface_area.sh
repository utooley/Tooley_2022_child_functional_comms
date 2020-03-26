
### FOR YEO7 AND YEO-DEVELOPMENTAL ###

declare -a array=(lh
rh)

#Yeo7
for hemi in "${array[@]}";
do
annot_file=/cbica/projects/spatial_topography/tools/parcellations/yeo7_from_freesurfer/fsaverage6/label/${hemi}.Yeo2011_7Networks_N1000.annot

#get surface area
mri_segstats --annot fsaverage6 lh $annot_file --i $SUBJECTS_DIR/fsaverage6/surf/${hemi}.white.avg.area.mgh --accumulate --sum /cbica/projects/spatial_topography/tools/parcellations/yeo7_from_freesurfer/fsaverage6/${hemi}.7nets.surfarea.stats
done


#Yeo dev

for hemi in "${array[@]}"
do
annot_file=/cbica/projects/spatial_topography/data/imageData/yeo_clustering_networks/yeo7_n670_2runsonly_1000tries/${hemi}.yeodev.fsaverage6.annot
#get surface area
mri_segstats --annot fsaverage6 lh $annot_file --i $SUBJECTS_DIR/fsaverage6/surf/${hemi}.white.avg.area.mgh --accumulate --sum freesurfer/${hemi}.surfarea.stats
done

### FOR WSBM AND YEO7 ###
#first, relabel the Schaefer fsaverage6 annotation to have WSBM labels with visualization.m, then calculate surface area.

declare -a array=(lh
rh)

for hemi in "${array[@]}"
do
annot_file=/cbica/projects/spatial_topography/data/imageData/wsbm/site16_training_sample/brains/${hemi}.wsbm.consensus.fsaverage6.annot

#get surface area
mri_segstats --annot fsaverage6 lh $annot_file --i $SUBJECTS_DIR/fsaverage6/surf/${hemi}.white.avg.area.mgh --accumulate --sum ${hemi}.surfarea.stats
done
