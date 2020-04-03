#Loop to submit matlab jobs
code_dir=/cbica/projects/spatial_topography/code/net_stats
subject_list=/cbica/projects/spatial_topography/data/subjLists/release2/site16/n670_subjects_only_filtered_runs_site16_postprocess.txt

x=0 #increment this by 1 each time people stop running.

mapfile -t myArray < $subject_list
begin=$(expr ${x} \* 25)
end=$(expr 25 + ${x} \* 25)

for i in "${myArray[@]:$begin:$end}";
do
echo $i
#qsub ${code_dir}/run_matlab_net_stats.sh ${subject}
done
