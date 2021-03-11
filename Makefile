cc = g++ -O3 -march=native
prom = SurfDisp
EIGEN = -I/mnt/d/linuxapp/eigen-3.3.9
cppsrc = $(shell find src/*.cpp)
objcpp = $(cppsrc:%.cpp=%.o) 

$(prom): $(objcpp)
	$(cc) -o $(prom) $(objcpp) -lm $(EIGEN)

%.o: %.cpp
	$(cc) -c  $(EIGEN) $< -o $@ 

clean:
	rm src/*.o

