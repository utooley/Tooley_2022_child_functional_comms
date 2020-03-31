
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
