"""
MIT License

Copyright (c) 2019 ming

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
===============================
@File  : test2.py
@Author: ming.ustb@outlook.com
@Date  : 19-5-27
@GitHub: https://github.com/yangmingustb/PTPSim
"""

import numpy as np
from math import cos, sin, pi, radians
from scipy.special import fresnel
import matplotlib.pyplot as plt


def spiralInterpolation(resolution, s, x, y, hdg, length, curvStart, curvEnd):
    points = np.zeros((int(length/resolution), 1))
    points = [i*resolution for i in range(len(points))]
    xx = np.zeros_like(points)
    yy = np.zeros_like(points)
    hh = np.zeros_like(points)
    if curvStart == 0 and curvEnd > 0:
        print("Case 1: curvStart == 0 and curvEnd > 0")
        radius = np.abs(1/curvEnd)
        A_sq = radius*length
        ss, cc = fresnel(np.square(points)/(2*A_sq*np.sqrt(np.pi/2)))
        xx = points*cc
        yy = points*ss
        hh = np.square(points)*2*radius*length
        xx, yy, hh = rotate(xx, yy, hh, hdg)
        xx, yy = translate(xx, yy, x, y)
        xx = np.insert(xx, 0, x, axis=0)
        yy = np.insert(yy, 0, y, axis=0)
        hh = np.insert(hh, 0, hdg, axis=0)

    elif curvStart == 0 and curvEnd < 0:
        print("Case 2: curvStart == 0 and curvEnd < 0")
        radius = np.abs(1/curvEnd)
        A_sq = radius*length
        ss, cc = fresnel(np.square(points)/(2*A_sq*np.sqrt(np.pi/2)))
        xx = points*cc
        yy = points*ss*-1
        hh = np.square(points)*2*radius*length
        xx, yy, hh = rotate(xx, yy, hh, hdg)
        xx, yy = translate(xx, yy, x, y)
        xx = np.insert(xx, 0, x, axis=0)
        yy = np.insert(yy, 0, y, axis=0)
        hh = np.insert(hh, 0, hdg, axis=0)

    elif curvEnd == 0 and curvStart > 0:
        print("Case 3: curvEnd == 0 and curvStart > 0")

    elif curvEnd == 0 and curvStart < 0:
        print("Case 4: curvEnd == 0 and curvStart < 0")

    else:
        print("The curvature parameters differ from the 4 predefined cases. Change curvStart and/or curvEnd")

    n_stations = int(length/resolution) + 1
    stations = np.zeros((n_stations, 3))
    for i in range(len(xx)):
        stations[i][0] = xx[i]
        stations[i][1] = yy[i]
        stations[i][2] = hh[i]

    return stations

def rotate(x, y, h, angle):
    # This function rotates the x and y vectors around zero
    xx = np.zeros_like(x)
    yy = np.zeros_like(y)
    hh = np.zeros_like(h)
    for i in range(len(x)):
        xx[i] = x[i]*cos(angle) - y[i]*sin(angle)
        yy[i] = x[i]*sin(angle) + y[i]*cos(angle)
        hh[i] = h[i] + angle
    return xx, yy, hh

def translate(x, y, x_delta, y_delta):
    # This function translates the x and y vectors with the delta values
    xx = np.zeros_like(x)
    yy = np.zeros_like(y)
    for i in range(len(x)):
        xx[i] = x[i] + x_delta
        yy[i] = y[i] + y_delta
    return xx, yy

stations = spiralInterpolation(1, 77, 50, 100, radians(56), 40, 0, 1/20)

x = []
y = []
h = []

for station in stations:
    x.append(station[0])
    y.append(station[1])
    h.append(station[2])

plt.figure(figsize=(20,13))
plt.plot(x, y, '.')
plt.grid(True)
plt.axis('equal')
plt.show()

def get_heading_components(x, y, h, length=1):
    xa = np.zeros_like(x)
    ya = np.zeros_like(y)
    for i in range(len(x)):
        xa[i] = length*cos(h[i])
        ya[i] = length*sin(h[i])
    return xa, ya

xa, ya = get_heading_components(x, y, h)
plt.figure(figsize=(20,13))
plt.quiver(x, y, xa, ya, width=0.005)
plt.grid(True)
plt.axis('equal')
plt.show()