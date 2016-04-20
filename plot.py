#!/usr/bin/env python3

from pylab import *
import sys

def smooth(y, box_pts=3):
    box = ones(box_pts)/box_pts
    y_smooth = convolve(y, box, mode='same')
    return y_smooth

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'data.out'

a = loadtxt('data.out')
t = a[:,0]
itp = a[:,1]
args = argsort(t)
t = t[args]
itp = itp[args]

wtd = gradient(gradient(itp, t[1]-t[0]), t[1]-t[0])
swtd = smooth(wtd)

# plot(t, wtd)
plot(t, wtd)
show()
