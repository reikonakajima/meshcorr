#!/usr/bin/env python

import subprocess as S
import numpy as np
import pylab as P
import sys
sys.path.append('/Users/reiko/2code/python/mylib')

from mpl_toolkits.mplot3d import Axes3D

"""
xyz.py
  plot x,y,z coordinates of GAMA (debug)

"""

datafname = 'comovingcoord3d.out'
(ra,dec,redshift,x,y,z) = np.loadtxt(datafname, unpack=True)

fig = P.figure()
ax = fig.add_subplot(111,projection='3d')

ax.scatter(x,y,z, c='r', marker='.')

for ii in xrange(0,360,1):
    ax.view_init(elev=10, azim=ii)
    P.savefig("movie%s"%ii+".png")



notetext = " "
P.figtext(0.52, 0.925, notetext, fontsize='xx-small', multialignment='left', weight='bold',)

annotatestr = S.Popen(['pwd'], stdout=S.PIPE).communicate()[0]
annotatestr = annotatestr + S.Popen(['date'], stdout=S.PIPE).communicate()[0]
annotatestr = annotatestr + sys.argv[0] + '::' + datafname
P.figtext(0.1, 0.93, annotatestr, fontsize='xx-small', multialignment='left',)

figfname = sys.argv[0].replace('.py', '.png')
P.savefig(figfname)
P.clf()

