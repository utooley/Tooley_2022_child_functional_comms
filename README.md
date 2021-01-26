
# Tooley et al., 2021

This repo contains code accompanying the manuscript [Functional brain network community structure in childhood: Unfinished territories and fuzzy boundaries](https://www.biorxiv.org/content/10.1101/2021.01.21.427677v1), as well as data imported directly into the manuscript.

We provide two freely available partitions (assignments of regions to communities, also sometimes called networks), that were estimated using data from children ages 9-11. For more details, see the [preprint](https://www.biorxiv.org/content/10.1101/2021.01.21.427677v1)!


## Code

Preprocessing was conducted using fMRIprep ([link](https://fmriprep.org/en/stable/)) and xcpEngine ([link](https://xcpengine.readthedocs.io/)). 

Code used in this repo depends heavily on that of the Computational Brain Imaging Group (found [here](https://github.com/ThomasYeoLab/CBIG)), which must be set up to run functions in the `code/yeo_networks` folder. The weighted stochastic block model (WSBM, `code/wsbm`) was run using code from Aicher et al. (2015), found [here](https://aaronclauset.github.io/wsbm/).

## Partitions

![image](cover_fig.png)

Both the developmental clustering partition and the developmental WSBM are provided in `fsaverage6`, `MNI152`, and `fs_LR 32k` space in the `partitions` folder. The WSBM partition is also provided as a vector of assignments for the Schaefer 400-region parcellation, available from the CBIG group [here](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal).

## Manuscript

Statistics and data that are directly imported into the manuscript live in this folder.
