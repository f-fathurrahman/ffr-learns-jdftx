JDFTX_HOME=/home/efefer/WORKS/JDFTX/jdftx-1.4.2

dirs="wannier phonon fluid electronic core commands"

for d in $dirs
do
  cp -v $JDFTX_HOME/jdftx/$d/*.cpp jdftx/$d
  cp -v $JDFTX_HOME/jdftx/$d/*.h jdftx/$d
done

