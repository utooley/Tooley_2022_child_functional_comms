#$ -l h_vmem=6.1G,s_vmem=6.0G
#$ -j y
#$ -q all.q,basic.q,himem.q,gpu.q
#$ -e /data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/output/job_output
source activate aws_cli
for line in `cat /data/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/release2/AWS_release2_site16_0.2mm_rest_T1s.txt`;
do
echo ${line}
aws s3 cp ${line} /data/jux/mackey_group/public_data/ABCD/release2_site16_dicoms
done


cd /data/jux/mackey_group/public_data/ABCD/release2_site16_dicoms
for line in `ls NDAR*`;
do
echo ${line}
tar --warning=no-timestamp -xvzf ${line}
done
