#!/bin/bash
set -e 

# clean directory
rm -rf out 
mkdir -p out 

sourcedir=../../

# # run sem
# for i in `seq 0 99`;
# do 
#     phi=`echo $i | awk '{print 0 + 360 / 99 * $1}'`
#     $sourcedir/bin/surftti modelhti.txt $phi 0.01 0.01 1

#     mv out/swd.txt out/swd.$i.txt 
# done
time $sourcedir/bin/surftti model.txt 0 0.01 0.5 100

# bin2h5
echo "converting to hdf5..."
python $sourcedir/scripts/binary2h5.py out/database.bin out/swd.txt out/kernels.h5


# run benckmark
\rm -f  *.so bench_cps.py 
ln -s $sourcedir/lib/cps* . 
\cp $sourcedir/scripts/bench_cps.py .
python bench_cps.py 2
mv out/swd.cps.txt out/swd.cps.rayl.txt 
python bench_cps.py 1
mv out/swd.cps.txt out/swd.cps.love.txt 

# plot dispersion curves
python plot_disp.py 

# plot eigenfunctions
python plot_kernels.py out/kernels.h5 50 0