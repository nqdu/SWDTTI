cc = g++ -O3 -march=native
cc1 = gfortran
prom = SurfDisp

f77flags= -ffixed-line-length-none -ffloat-store -W  -fbounds-check
EIGEN = -I/mnt/d/linuxapp/eigen-3.3.9
cppsrc = $(shell find src/*.cpp)
objcpp = $(cppsrc:%.cpp=%.o) 

$(prom): $(objcpp)
	$(cc) -o $(prom) $(objf77) $(objcpp) -lm $(EIGEN)

%.o: %.cpp
	$(cc) -c  $(EIGEN) $< -o $@ 

clean:
	rm src/*.o

