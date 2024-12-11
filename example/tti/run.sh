#!/bin/bash
set -e 

# clean directory
rm -rf out 
mkdir -p out 

sourcedir=../../

# run sem
for i in `seq 0 99`;
do 
    phi=`echo $i | awk '{print 0 + 360 / 99 * $1}'`
    $sourcedir/bin/surftti modelhti.txt $phi 0.01 0.01 1

    mv out/swd.txt out/swd.$i.txt 
done
# $sourcedir/bin/surftti modelhti.txt 45 0.01 0.01 1

