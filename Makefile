#include make.inc
CXX = g++

#CXX_FLAGS = -O3 -Wall -Wno-unused-result -ftemplate-depth-512 -std=c++0x
CXX_FLAGS = -O0 -g -Wall -Wno-unused-result -ftemplate-depth-512 -std=c++0x

CXX_DEFINES = -DMPI_ENABLED

CXX_INCLUDES = \
-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi \
-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent \
-I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include \
-I/usr/lib/x86_64-linux-gnu/openmpi/include \
-I/home/efefer/WORKS/JDFTX/jdftx-1.4.2/build \
 -I./jdftx

LIBS_EXT =

include make.sources

include make.objs

# Libraries
lib: $(OBJS)
	ar rcs libjdftx.a objs/*/*.o

jdftx: lib jdftx.cpp
	$(CXX) $(CXX_FLAGS) $(CXX_DEFINES) $(CXX_INCLUDES) -Wl,-rpath -Wl,/usr/lib/openmpi/lib -Wl,--enable-new-dtags -pthread jdtfx.cpp -o jdftx libjdftx.so -Wl,-rpath,/usr/lib/openmpi/lib /usr/lib/openmpi/lib/libmpi_cxx.so /usr/lib/openmpi/lib/libmpi.so -lgsl -lfftw3_threads -lfftw3 -lcblas -lopenblas -lpthread

clean:
	rm -rv objs/*/*.o

include make.rules

