#$ -j y
#$ -l h_vmem=10.1G,s_vmem=10.0G
#$ -o /data/picsl/mackey_group/CBPD/output/qsub_output
#$ -q himem.q,all.q,basic.q,gpu.q

subID=${1}

singularity run -B /data/picsl/mackey_group:/mnt --cleanenv /data/picsl/mackey_group/tools/singularity/heudiconv0.5.4.simg -d /mnt/BPD/dicoms/{subject}/*.dcm -o /mnt/CBPD/CBPD_bids -f /mnt/CBPD/bids_ppc_scripts/heudiconv/heuristic.py -s ${subID} -ss 01 -c dcm2niix -b --minmeta;
