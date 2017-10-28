prefix=`basename $1 .cpp`

FLAGS="-Wall -O3 -ftemplate-depth-512 -std=c++0x"

INCLUDE="-I./jdftx"

LIBS=

g++ $INCLUDE $FLAGS $1 -o $prefix.x
