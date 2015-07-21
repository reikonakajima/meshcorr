#
# SConstruct
# 
import os
import sys

home_dir = os.getcwd()


# CPPPATH: sets the include path for C++
env = Environment()
env.Append(CPPPATH=['/vol/software/software/tools/tmv/tmv0.72/x86_64/include'])
env.Append(CPPPATH=['/opt/local/include'])
env.Append(CPPPATH=[os.path.join(home_dir,'src'),
		    os.path.join(home_dir,'src','utilities'),
		    os.path.join(home_dir,'src','cosmology'),
		    ])
env.Append(LIBS=['tmv', 'blas', 'CCfits'])
env.Append(LIBPATH = ['/vol/software/software/tools/tmv/tmv0.72/x86_64/lib', '/usr/local/lib'])
env.Append(CPPPATH=['-stdlib=libstdc++'])

def ReadFileList(fname):
    """
    This reads a list of whitespace separated values from the input file fname
    and stores it as a list.  We will make this part of the environment so
    other SConscripts can use it
    """
    try:
        files=open(fname).read().split()
    except:
        print 'Could not open file:',fname
	sys.exit(45)
    files = [f.strip() for f in files]
    return files

# subdirectory SConscript files can use this function
env['_ReadFileList'] = ReadFileList



# build the sub-objects
Export("env")
SConscript("src/SConscript")

# append options specific to the head node
env.Append(CPPPATH=['utilities',])

# specify the sub-objects
sub_objects = '''
	    src/utilities/.obj/StringStuff.o
	    src/.obj/GalaxyObjects.o
	    src/.obj/GAMAObjects.o
	    '''.split()

# build the main programs

#test = env.Program(target='GGLensTest', source=sub_objects+['src/GGLensDriver.cpp',])
#starhalo = env.Program(target='gglens_starhalo', source=sub_objects+['src/gglens_starhalo.cpp',])
#env.Install('bin', [test, starhalo])
#env.Alias('install', 'bin')

# (eventually, build a library code)
