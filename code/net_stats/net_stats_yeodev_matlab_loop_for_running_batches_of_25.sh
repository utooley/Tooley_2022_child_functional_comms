#Loop to submit matlab jobs for running full corr mats with yeo-dev
code_dir=/cbica/projects/spatial_topography/code/net_stats
subject_list=/cbica/projects/spatial_topography/data/subjLists/release2/site16/n670_subjects_only_filtered_runs_site16_postprocess.txt

x=26 #increment this by 1 each time people stop running. Run in batches of 25

mapfile -t myArray < $subject_list
begin=$(expr ${x} \* 25)

for i in "${myArray[@]:${begin}:25}";
do
echo $i
qsub ${code_dir}/run_matlab_net_stats.sh ${i}
done


## LOOP TO AGGREGATE SUBJECT DATA AND SUBMIT NEW JOBS FOR THEM

code_dir=/cbica/projects/spatial_topography/code/net_stats
subject_list=/cbica/projects/spatial_topography/data/subjLists/release2/site16/n670_subjects_only_filtered_runs_site16_postprocess.txt
data_dir=/cbica/projects/spatial_topography/data/imageData/net_stats/site16_subjectwise

mapfile -t myArray < $subject_list

for i in "${myArray[@]}";
do
echo $i
if [ -e ${data_dir}/${i}_site16_fsaverage6_yeodev_network_stats.csv ]; then
  cat ${data_dir}/${i}_site16_fsaverage6_yeodev_network_stats.csv | head -2 | tail -1 >> ${data_dir}/INCOMPLETE_n670_training_sample_fsaverage6_yeodev_network_stats.csv
else
  echo ${i} 'doesnt exist'
  qsub ${code_dir}/run_matlab_net_stats.sh ${i}
fi
done
