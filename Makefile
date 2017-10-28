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
	ar rcs libjdftx.a objs/*/*.o

jdftx: lib jdftx.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_DEFINES) $(CXX_INCLUDES) -Wl,-rpath  -Wl,/usr/lib/openmpi/lib  -Wl,--enable-new-dtags jdftx.cpp  -o jdftx.x  -Wl,--whole-archive libjdftx.a -Wl,--no-whole-archive /usr/lib/openmpi/lib/libmpi_cxx.so /usr/lib/openmpi/lib/libmpi.so /home/efefer/mysoftwares/lib/libgsl.a /usr/lib/x86_64-linux-gnu/libfftw3_threads.a /home/efefer/mysoftwares/lib/libfftw3.a /home/efefer/intel/mkl/lib/intel64/libmkl_intel_lp64.so /home/efefer/intel/mkl/lib/intel64/libmkl_gnu_thread.so /home/efefer/intel/mkl/lib/intel64/libmkl_core.so -lgomp -lpthread -Wl,-rpath,/usr/lib/openmpi/lib

clean:
	rm -rv objs/*/*.o

include make.rules
