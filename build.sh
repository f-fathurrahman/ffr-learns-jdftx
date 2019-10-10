prefix=`basename $1 .cpp`

FLAGS="-Wall -O0 -g -ftemplate-depth-512 -std=c++0x"

INCLUDE="-I./jdftx"

#LIBS="libjdftx.a -lgsl -lfftw3_threads -lfftw3 -lcblas -lopenblas -lpthread"
LIBS="-L ./ -ljdftx -lgsl -lfftw3_threads -lfftw3 -lcblas -lopenblas -lpthread"

mpic++ $INCLUDE $FLAGS $1 -o $prefix.x $LIBS
