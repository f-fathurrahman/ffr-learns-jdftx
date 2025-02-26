MAKE = make
AR = ar

CXX = g++
CXX_INC = -I./jdftx
CXX_OPTS = -Wall -O0 -g -ftemplate-depth-512 -std=c++0x
CXX_OPTS_DYN = $(CXX_OPTS) -fPIC

LIBS = -Wl,-rpath,/home/efefer/WORKS/my_github_repos/ffr-learns-jdftx: libjdftx.so \
-lgsl -lfftw3_threads -lfftw3 -lgslcblas -llapack -lblas -lpthread


CXX_SRC = \
MyElecGradient.cpp \
my_elecEnergyAndGrad.cpp my_linminQuad.cpp \
MyElecMinimizer.cpp my_elecFluidMinimize.cpp \
write_info.cpp \
my_calcDensity.cpp my_applyHamiltonian.cpp \
my_ElecMinimizer_step.cpp \
my_ElecMinimizer_minimize.cpp \
my_ElecMinimizer_compute.cpp \
initialize_Haux.cpp \
simple_minimize.cpp \
export_variables.cpp

OBJ = $(CXX_SRC:.cpp=.o)

%.o : %.cpp
	$(CXX) $(CXX_INC) $(CXX_OPTS) -c -o $(*F).o $<


# Targets
lib: $(OBJ)
	ar rcs libjdftx_debug.a *.o

clean:
	rm -rf *.o libjdftx_debug.a *.x


