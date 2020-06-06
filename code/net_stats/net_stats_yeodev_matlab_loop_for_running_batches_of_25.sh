#Loop to submit matlab jobs for running full corr mats with yeo-dev
code_dir=/cbica/projects/spatial_topography/code/net_stats
subject_list=/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/n544_subjects_only_filtered_runs_site14site20_postprocess.txt

x=23 #increment this by 1 each time people stop running. Run in batches of 25, stop at 22

# FOR SITE14SITE20 subject list
for subject in `cat ${subject_list} | cut -d, -f2 | uniq | head -$(expr ${x} \* 25)| tail -25`
do
echo $subject
qsub ${code_dir}/run_matlab_net_stats.sh ${subject}
done

# FOR SITE 16 subject list
mapfile -t myArray < $subject_list
begin=$(expr ${x} \* 25)

for i in "${myArray[@]:${begin}:25}";
do
echo $i
qsub ${code_dir}/run_matlab_net_stats.sh ${i}
done

## LOOP TO AGGREGATE SUBJECT DATA AND SUBMIT NEW JOBS FOR THEM

code_dir=/cbica/projects/spatial_topography/code/net_stats
subject_list=/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/n544_subjects_only_filtered_runs_site14site20_postprocess.txt
data_dir=/cbica/projects/spatial_topography/data/imageData/net_stats/site14site20_subjectwise

mapfile -t myArray < $subject_list

for i in "${myArray[@]}";
do
echo $i
if [ -e ${data_dir}/${i}_site16_fsaverage6_yeodev_network_stats.csv ]; then
  cat ${data_dir}/${i}_site16_fsaverage6_yeodev_network_stats.csv | head -2 | tail -1 >> ${data_dir}/n544_test_sample_fsaverage6_yeodev_network_stats.csv
else
  echo ${i} 'doesnt exist'
  qsub ${code_dir}/run_matlab_net_stats.sh ${i}
fi
done
