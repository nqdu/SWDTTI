#!/bin/bash
set -e 

# clean directory
rm -rf out 
mkdir -p out 

sourcedir=../../

# run sem
time $sourcedir/bin/surftti modelhti.txt
