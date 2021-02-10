# Partitions

![image](../cover_fig.png)
Both the developmental clustering partition and the developmental WSBM are provided in `fsaverage6`, `MNI152`, and `fs_LR 32k` space. The adult clustering partition was estimated by Yeo _et al._ in their [2011 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3174820/), and is available from the CBIG group [here](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Yeo2011_fcMRI_clustering).

## fsaverage6
This folder contains partitions in the `fsaverage6` space. The structure of this folder follows that of a preprocessed Freesurfer subject. In particular, fsaverage6/label/ contains the partition files, and fsaverage6/surf contains the fsaverage6 surfaces needed to visualize the partition.

To visualize, make sure Freesurfer is set up and configured (http://surfer.nmr.mgh.harvard.edu/fswiki/SetupConfiguration), and `cd` to this folder.

Run the following command:

`freeview -f fsaverage6/surf/rh.inflated:annot=fsaverage6/label/rh.clustering.fsaverage6.annot`

## MNI152

The developmental clustering partition provided was projected to MNI152 1mm space from fsaverage6 (using the RF-ANTS method published by Wu *et al.* (2018) in *Human Brain Mapping*). The WSBM partition is provided as a relabeling of the Schaefer-400 parcellation in MNI152 1mm space (available from the CBIG group [here](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal)).


The WSBM is also provided as a vector of assignments for the Schaefer 400-region parcellation, in `wsbm_schaefer400_assignments.txt`.

## fsLR

Partitions are provided in fsLR 32k space, transformed from fsaverage6 space following instructions from the HCP group [here](https://wiki.humanconnectome.org/download/attachments/63078513/Resampling-FreeSurfer-HCP_5_8.pdf).

