cd /data/picsl/mackey_group/public_data/ABCD/site14_2/submission_16889
for line in `ls`;
do
echo ${line}
#tar --warning=no-timestamp -xvzf ${line}
tar -xvzf ${line}
done
