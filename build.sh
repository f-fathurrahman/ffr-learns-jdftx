prefix=`basename $1 .cpp`

FLAGS="-Wall -O0 -g -ftemplate-depth-512 -std=c++0x"

INCLUDE="-I./jdftx"

LIBS="libjdftx.a"

g++ $INCLUDE $FLAGS $1 -o $prefix.x
