dirs="wannier phonon fluid electronic core commands"
for d in $dirs
do
  cp -v ../jdftx-git/jdftx/$d/*.h jdftx/$d
done
