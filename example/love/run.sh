#!/bin/bash
set -e 

# clean directory
rm -rf out 
mkdir -p out 

sourcedir=../../

# run sem
$sourcedir/bin/surf24 model.txt 1

# run benckmark
python bench_cps.py

# plot 
python $sourcedir/plot/plot_disp.py
python plot_kernel.py 119.txt

