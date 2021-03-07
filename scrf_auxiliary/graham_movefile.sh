#!/bin/bash
#SBATCH --job-name=movescrf
#SBATCH --time=0-12:00:00
#SBATCH --mem=1G
#SBATCH --account=rpp-kshook


path_0=/home/gut428/scratch/SCRF/rndnum_tmean
path_tar=/home/gut428/projects/rpp-kshook/gut428/EMGLB/SCRF/rndnum_tmean

y1=1950
y2=2019
e1=1
e2=20
m1=1
m2=12

for ((e=e1;e<=e2;e++))
do
for ((y=y1;y<=y2;y++))
do
for ((m=m1;m<=m2;m++))
do
file=${path_0}/scrf_$((y*100+m))_$(printf %03d $e).nc
echo move - $file - to - $path_tar
start_time="$(date -u +%s)"
mv $file $path_tar
end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed for process"
done
done
done
