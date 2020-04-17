
subjlist_dir=/cbica/projects/spatial_topography/data/subjLists/release2/site16/yeo_networks
input=${subjlist_dir}/n670_site16_fullpaths_surf2surf_rh.txt

subjlist_dir=~/Downloads
input=~/Downloads/n670_site16_fullpaths_surf2surf_lh.txt
while read -r line
do
echo $line
sub=$(echo $line | cut -d / -f 10)
filename=${subjlist_dir}/subjects/${sub}_fullpaths_surf2surf_lh.txt
echo $line > $filename
done < $input


## FOR site14site20

subjlist_dir=/cbica/projects/spatial_topography/data/subjLists/release2/site14site20/
input=${subjlist_dir}/n546_filtered_runs_site14site20_postprocess.csv

## WHEN SUBBING LH FOR RH, MAKE SURE IS CASE SENSITIVE!
subjlist_dir=~/Downloads
input=~/Downloads/n546_filtered_runs_site14site20_fullpaths_surf2surf_rh.txt
while read -r line
do
echo $line
# newline=$(echo $line | cut -d , -f 7)
# newline2=$(echo $line | cut -d , -f 8)
sub=$(echo $line | cut -d / -f 10)
echo $sub
filename=${subjlist_dir}/subjects_site14/${sub}_fullpaths_surf2surf_rh.txt
echo $line > $filename
done < $input
