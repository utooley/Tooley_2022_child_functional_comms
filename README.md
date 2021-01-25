
## Tooley et al., 2021

This repo contains code accompanying the manuscript *Functional brain network community structure in childhood: Unfinished territories and fuzzy boundaries*, as well as data imported directly into the manuscript. We provide two freely available partitions (assignments of regions to communities, also sometimes called networks), that were estimated using data from children ages 9-11.

For more details, see the preprint [here](https://www.biorxiv.org/content/10.1101/2021.01.21.427677v1)!


## Code

Preprocessing was conducted using fMRIprep (link) and xcpEngine (link. Code used in this repo depends heavily on that of the Computational Brain Imaging Group (https://github.com/ThomasYeoLab/CBIG), which must be set up to run functions in the `code/yeo_networks` folder. The weighted stochastic block model (WSBM) was employed using code from Aicher et al. (2015), found [here](https://aaronclauset.github.io/wsbm/).

## Partitions

(Image here of brains)

Both the developmental clustering partition and the developmental WSBM are provided in xxx and xxx space in `brains`. The WSBM is also provided as a vector of assignments for the Schaefer xx parcellation, available from the CBIG group [here](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal).

## Manuscript

Statistics and data that are directly imported into the manuscript live in this folder.