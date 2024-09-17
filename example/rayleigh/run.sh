#!/bin/bash
set -e 

# clean directory
rm -rf out 
mkdir -p out 

sourcedir=../../

# run sem
time $sourcedir/bin/surf24 model.txt 2

# run benckmark
python bench_cps.py

# plot 
python $sourcedir/plot/plot_disp.py
#python plot_kernel.py 119.txt
