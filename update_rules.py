from __future__ import print_function

import subprocess
import os

cmd = 'ls jdftx/*/*.cpp'
lines = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout.read()

# remove various unwanted files
lines = lines.replace('jdftx/wannier/main.cpp','')
lines = lines.replace('jdftx/phonon/main.cpp','')

Nlines = len(lines)

f = open('make.sources','w')
f.write('SRC = \\\n')
ilines = 0
for l in lines.split():
    ilines = ilines + 1
    f.write(l)
    if(ilines != Nlines):
        f.write(' \\\n')
    else:
        f.write('\n')
f.close()

f = open('make.objs','w')
f.write('OBJS = \\\n')
ilines = 0
for l in lines.split():
    ilines = ilines + 1
    objfile = l.replace('jdftx/','').replace('.cpp','.o')
    f.write('objs/' + objfile)
    if(ilines != Nlines):
        f.write(' \\\n')
    else:
        f.write('\n')
f.close()


f = open('make.rules','w')
for l in lines.split():
    objfile = l.replace('jdftx/','').replace('.cpp','.o')
    f.write('objs/' + objfile + ':\n')
    f.write('\t$(CXX) -c $(CXX_INCLUDES) $(CXX_FLAGS) $(CXX_DEFINES) ' + l + ' -o objs/' + objfile + '\n\n')
f.close()

#print(retval)
#print(type(retval))
