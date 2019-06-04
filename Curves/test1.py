import numpy as np
from math import cos, sin, pi, radians, sqrt
from scipy.special import fresnel
import matplotlib.pyplot as plt


def spiral_interp_centre(distance, x, y, hdg, length, curvEnd):
    '''Interpolate for a spiral centred on the origin'''
    # s doesn't seem to be needed...
    theta = hdg                    # Angle of the start of the curve
    Ltot = length                  # Length of curve
    Rend = 1 / curvEnd             # Radius of curvature at end of spiral

    # Rescale, compute and unscale
    a = 1 / sqrt(2 * Ltot * Rend)  # Scale factor
    distance_scaled = distance * a # Distance along normalised spiral
    deltay_scaled, deltax_scaled = fresnel(distance_scaled)
    deltax = deltax_scaled / a
    deltay = deltay_scaled / a

    # deltax and deltay give coordinates for theta=0
    deltax_rot = deltax * cos(theta) - deltay * sin(theta)
    deltay_rot = deltax * sin(theta) + deltay * cos(theta)

    # Spiral is relative to the starting coordinates
    xcoord = x + deltax_rot
    ycoord = y + deltay_rot

    return xcoord, ycoord

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# This version
xs = []
ys = []
for n in range(-100, 100+1):
    x, y = spiral_interp_centre(n, 50, 100, radians(56), 40, 1/20.)
    xs.append(x)
    ys.append(y)
ax.plot(xs, ys)

'''
# Your version
from Curves.test2 import spiralInterpolation
stations = spiralInterpolation(1, 77, 50, 100, radians(56), 40, 0, 1/20.)
ax.plot(stations[:,0], stations[:,1])

ax.legend(['My spiral', 'Your spiral'])
fig.savefig('spiral.png')
'''
plt.show()
