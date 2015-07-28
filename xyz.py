#!/usr/bin/env python

import subprocess as S
import numpy as np
import pylab as P
import sys
sys.path.append('/Users/reiko/2code/python/mylib')

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

"""
xyz.py
  plot x,y,z coordinates of GAMA (debug)

"""

fig = P.figure()
ax = fig.add_subplot(111,projection='3d')

datafname0 = 'comovingcoord3d.out'
(ra0,dec0,redshift0,x0,y0,z0) = np.loadtxt(datafname0, unpack=True)

datafname = 'test.out'
(ra,dec,redshift,x,y,z) = np.loadtxt(datafname, unpack=True)

ax.scatter(x0,y0,z0, c='r', marker='.')
ax.scatter(x,y,z, c='b', marker='x')

# parametric plot
xx = -300.
yy = 400.
zz = -90.
epsilon = 0.1
theta = np.linspace(np.pi/2.-epsilon, np.pi/2.+epsilon, 20)
phi = np.arctan2(yy, xx)
x1 = np.cos(phi) * np.cos(theta) * 500.
y1 = np.sin(phi) * np.cos(theta) * 500.
z1 = np.sin(theta) * 500.
ax.plot(x1,y1,z1, 'b-')

for ii in xrange(0,360,1):
    ax.view_init(elev=10., azim=ii)
    ax.set_xlim3d(-800,0)
    ax.set_ylim3d(0,900)
    ax.set_zlim3d(0,20)
    P.savefig("xyz%04d"%ii + ".png")


shcommand = "ffmpeg -i xyz%04d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4"
S.call(shcommand, shell=True)


"""
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360, interval=20, blit=True)
anim.save('xyz.mp4', fps=30,) # extra_args=['-vcodec', 'libx264',])
"""
"""
notetext = " "
P.figtext(0.52, 0.925, notetext, fontsize='xx-small', multialignment='left', weight='bold',)

annotatestr = S.Popen(['pwd'], stdout=S.PIPE).communicate()[0]
annotatestr = annotatestr + S.Popen(['date'], stdout=S.PIPE).communicate()[0]
annotatestr = annotatestr + sys.argv[0] + '::' + datafname
P.figtext(0.1, 0.93, annotatestr, fontsize='xx-small', multialignment='left',)

figfname = sys.argv[0].replace('.py', '.png')
P.savefig(figfname)
P.clf()
"""
