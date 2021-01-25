#$ -j y
#$ -l h_vmem=10.1G,s_vmem=10.0G
#$ -o /data/picsl/mackey_group/CBPD/output/qsub_output
#$ -q himem.q,all.q,basic.q,gpu.q

subID=${1}

singularity run -B /data/jux/mackey_group/public_data/ABCD/release2_site16_dicoms:/mnt --cleanenv /data/picsl/mackey_group/tools/singularity/heudiconv0.5.4.simg -d /mnt/sub-{subject}/ses-{session}/*/*/sub-{subject}*.dcm -o /mnt/niftis -f /mnt/4_heuristic_ABCD.py -s ${subID} -ss baselineYear1Arm1 -c dcm2niix -b;
