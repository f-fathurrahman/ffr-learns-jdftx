prefix=`basename $1 .cpp`

FLAGS="-Wall -O0 -g -ftemplate-depth-512 -std=c++0x"

INCLUDE="-I./jdftx"

#LIBS="libjdftx.a -lgsl -lfftw3_threads -lfftw3 -lcblas -lopenblas -lpthread"
#LIBS="-L ./ -ljdftx -lgsl -lfftw3_threads -lfftw3 -lcblas -lopenblas -lpthread"

LIBS="-Wl,-rpath,/home/efefer/WORKS/my_github_repos/ffr-learns-jdftx: libjdftx.so -lgsl -lfftw3_threads -lfftw3 -lgslcblas -llapack -lblas -lpthread"

g++ $INCLUDE $FLAGS -c $1
