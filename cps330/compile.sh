# macos
CC="g++-14 -O3 -march=native -shared"
FC="gfortran -O3 -march=native -shared "
include=`python -m pybind11 --includes`
FFLAGS="-ffixed-line-length-none "
outname="../libsurf"`python3.11-config --extension-suffix`
set -x 
$FC $FFLAGS -g -c surfdisp96.f -o surfdisp96.o
#$CC -g -c surfdisp.cpp -o surfdisp.o 
$CC -g -c main.cpp -o main.o  $include
$CC *.o -o $outname -lgfortran  -undefined dynamic_lookup
rm *.o 


# linux
# CC="g++ -O3 -march=native -shared -fPIC"
# FC="gfortran -O3 -march=native -shared -fPIC"
# include=`python -m pybind11 --includes`
# FFLAGS="-ffixed-line-length-none "
# outname="../libsurf"`python3-config --extension-suffix`
# set -x 
# $FC $FFLAGS -g -c surfdisp96.f -o surfdisp96.o
# #$CC -g -c surfdisp.cpp -o surfdisp.o 
# $CC -g -c main.cpp -o main.o  $include
# $CC *.o -o $outname -lgfortran
# rm *.o