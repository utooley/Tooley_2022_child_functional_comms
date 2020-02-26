#run this in the /surfaces directory, and then check for mincost values that are too high!
for sub in `ls`
do
echo ${sub}\n >> mincost_site16.csv
cat ${sub}/**/*.mincost >> mincost_site16.csv
done
