dirs="wannier phonon fluid electronic core commands"
for d in $dirs
do
  cp -v ../jdftx-git/jdftx/$d/*.cpp jdftx/$d
done
