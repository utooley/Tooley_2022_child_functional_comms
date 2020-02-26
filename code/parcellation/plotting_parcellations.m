CBIG_DrawSurfaceMaps(results.lh_label, results.rh_label, 'fsaverage6', 'inflated')

#this only works in fsaverage5 or fsaverage space.
CBIG_SaveParcellationToFreesurferAnnotation('/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/gwMRF_cluster_test/clustering/results_seed1.mat', '/data/jux/mackey_group/Ursula/lh.out', '/data/jux/mackey_group/Ursula/rh.out')

#an example CT structure
[vertices, label, example] = read_annotation('/data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/gwMRF_cluster_test/clustering/lh.Schaefer2018_400Parcels_7Networks_order.annot')

#however many parcels you have in each hemisphere. numEntries must be this + 1.
my_ct = ones(50+1, 3);
my_ct(2:end, :) = round(rand(50, 3)*255);
for i = 1:ct.numEntries
  r = my_ct(i, 1);
  g = my_ct(i, 2);
  b = my_ct(i, 3);
  ct.table(i, :) = [r g b 0 r + g*2^8 + b*2^16];
end

lh_labels_annot = lh_labels;
for i=0:(ct.numEntries-1)
    lh_labels_annot(lh_labels == i) = ct.table(i+1, 5);
end

ct.struct_names=num2cell(1:51) #this doesn't quite work to get names in, work on this.
ct.orig_tab='Mycolortable'

Write_Brain_Annotation('lh.annot', [0:40961], lh_labels_annot, ct)



#other figure commands
abs_path_to_lh_annot_file = fullfile('/Applications/freesurfer','subjects', 'fsaverage6', 'label', 'lh.aparc.annot')
abs_path_to_rh_annot_file = fullfile(getenv('FREESURFER_HOME'), 'subjects', 'fsaverage6', 'label', 'rh.aparc.annot');
ref_medial_wall_label = 0;

CBIG_DrawSurfaceDataAsAnnotation(results.lh_label, results.rh_label, 'fsaverage6', abs_path_to_lh_annot_file, abs_path_to_rh_annot_file, ref_medial_wall_label, '~/Downloads', 'testcluster')
