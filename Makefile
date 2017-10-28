#include make.inc
CXX = g++
CXX_FLAGS =   -Wall -Wno-unused-result -ftemplate-depth-512    -std=c++0x
CXX_DEFINES = -DMKL_PROVIDES_BLAS -DMPI_ENABLED -DTHREADED_BLAS
CXX_INCLUDES = -I/home/efefer/mysoftwares/include -I/home/efefer/intel/mkl/include -I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent -I/usr/lib/openmpi/include/openmpi/opal/mca/event/libevent2021/libevent/include -I/usr/lib/openmpi/include -I/usr/lib/openmpi/include/openmpi -I./jdftx

LIBS_EXT =

include make.sources

include make.objs

# Libraries
lib: $(OBJS)
	ar rcs libmain.a objs/*.o

clean:
	rm -rv objs/*/*.o

include make.rules
