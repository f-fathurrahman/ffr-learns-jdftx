CXX = g++

#CXX_FLAGS = -O3 -Wall -Wno-unused-result -ftemplate-depth-512 -std=c++0x

# No optimization
CXX_FLAGS = -Wall -Wno-unused-result -ftemplate-depth-512 -std=c++0x

CXX_DEFINES = -DMKL_PROVIDES_BLAS -DMPI_ENABLED -DTHREADED_BLAS

CXX_INCLUDES = -I../ \
-I/home/efefer/mysoftwares/include -I/home/efefer/intel/mkl/include \
-I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent \
-I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent/include \
-I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi

LIBS_EXT =

include make.sources

include make.objs

# Libraries
lib: $(OBJS)
	cp -v *.o ../../objs

clean:
	rm -rv *.o

%.o : %.cpp
	$(CXX) $(CXX_INCLUDES) $(CXX_FLAGS) $(CXX_DEFINES) -c -o $(*F).o $<



