from __future__ import print_function

import subprocess
import os

cmd = 'ls jdftx/*/*.cpp'
lines = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE).stdout.read()

# remove various unwanted files
lines = lines.replace('jdftx/electronic/Blip.cpp','')
lines = lines.replace('jdftx/electronic/matrix.cpp','')
lines = lines.replace('jdftx/electronic/operators.cpp','')
lines = lines.replace('jdftx/electronic/RadialFunction.cpp','')
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


# a dumb way to figure out dependencies by analyzing include statements
# does not work for preprocessing
def find_deps(filename):
    #
    f = open(filename)
    #
    include_files = []
    while True:
        #
        line = f.readline()
        if not line:
            break
        #
        # Start processing
        #
        # heuristic for this project to determined whether the dependencies is
        # internal or external (system's *.h files)
        if ("#include" in line) & ("/" in line) & (not "gsl" in line) & (not "sys/" in line):
            #
            # special case for this project, append jdtfx/
            #
            l = "jdftx/" + line.split()[1].replace("<","").replace(">","")
            #
            include_files.append(l)
    #
    f.close()
    return include_files

def find_all_deps(filename):
    inc_one = find_deps(filename)
    #
    inc_all = []
    inc_all.append(inc_one)
    #
    for f in inc_one:
        inc_all.append( find_deps(f) )
    #
    # Flatten the list
    #
    flat_list = []
    for sublist in inc_all:
        for l in sublist:
            flat_list.append(l)
    return flat_list


f = open('make.rules','w')
for l in lines.split():
    objfile = l.replace('jdftx/','').replace('.cpp','.o')
    #
    f.write('objs/' + objfile + ': ' + l + ' ')
    # find deps
    print('Finding deps for file ' + l)
    deps_list = find_all_deps(l)
    #
    for fil in deps_list:
        f.write(fil + ' ')
    f.write('\n')
    #
    # compilation command
    #
    f.write('\t$(CXX) -c $(CXX_INCLUDES) $(CXX_FLAGS) $(CXX_DEFINES) ' + l + ' -o objs/' + objfile + '\n\n')
f.close()


"""
f = open('make.rules','w')
for l in lines.split():
    objfile = l.replace('jdftx/','').replace('.cpp','.o')
    f.write('objs/' + objfile + ': ' + l + '\n')
    f.write('\t$(CXX) -c $(CXX_INCLUDES) $(CXX_FLAGS) $(CXX_DEFINES) ' + l + ' -o objs/' + objfile + '\n\n')
f.close()
"""


#print(retval)
#print(type(retval))
