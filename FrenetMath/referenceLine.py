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
@File  : referenceLine.py
@Author: ming.ustb@outlook.com
@Date  : 19-5-29
@GitHub: https://github.com/yangmingustb/PTPSim
"""


import math
from scipy.linalg import solve
import numpy as np
import matplotlib.pyplot as plt


showLane=True

LaneWidth = 3.75  # [m]
left = 1


def wayPointDistribution(rx, ry, ryaw, s):
	"""

	:param rx:
	:param ry:
	:param ryaw:
	:param s:
	:return: generate the efficients of the reference line
	"""

	x_list = []
	y_list = []
	theta_list = []
	s_list = []
	for i in range(len(rx)):
		if 20 * i > (len(rx) - 1):
			break
		x_list.append(rx[20 * i])
		y_list.append(ry[20 * i])
		theta_list.append(ryaw[20 * i])
		s_list.append(s[20 * i])

	x_list.append(rx[-1])
	y_list.append(ry[-1])
	theta_list.append(ryaw[-1])
	s_list.append(s[-1])

	efficients = []
	for i in range(len(x_list) - 1):
		x_init = [x_list[i], y_list[i], theta_list[i]]
		x_goal = [x_list[i + 1], y_list[i + 1], theta_list[i + 1]]
		s0 = s_list[i]
		sf = s_list[i + 1]

		b, c = arcLengthCurve(x_init, x_goal, s0, sf, left)
		tmp_efficients = [s0]
		tmp_efficients.append(b)
		tmp_efficients.append(c)
		efficients.append(tmp_efficients)

	return efficients


def arcLengthCurve(x_init, x_goal, s0, sf, rightLeft):
	"""
	参考线弧长参数化
	:param x_init:
	:param x_goal:
	:param s0:
	:return:
	"""
	x0 = x_init[0]
	y0 = x_init[1]
	theta0 = x_init[2]

	xg = x_goal[0]
	yg = x_goal[1]
	thetag = x_goal[2]

	B = np.array([[1, s0, s0 ** 2, s0 ** 3],
	              [1, sf, sf ** 2, sf ** 3],
	              [0, 1, 2 * s0, 3 * s0 ** 2],
	              [0, 1, 2 * sf, 3 * sf ** 2]])

	b_b = np.array([x0, xg, np.cos(theta0), np.cos(thetag)])

	b_vector = solve(B, b_b)
	# print(a)
	# 最右侧车道线
	s = np.arange(s0, sf, 0.01)
	fun_x = b_vector[0] + b_vector[1] * s + b_vector[2] * s ** 2 + b_vector[3] * s ** 3
	d_funx = b_vector[1] + 2 * b_vector[2] * s + 3 * b_vector[3] * s ** 2

	c_b = np.array([y0, yg, np.sin(theta0), np.sin(thetag)])
	c_vector = solve(B, c_b)
	# print(a)
	fun_y = c_vector[0] + c_vector[1] * s + c_vector[2] * s ** 2 + c_vector[3] * s ** 3
	d_funy = c_vector[1] + 2 * c_vector[2] * s + 3 * c_vector[3] * s ** 2

	theta = []
	for i in range(len(d_funx)):
		tmp_theta = math.atan2((d_funy[i]),(d_funx[i]))
		theta.append(tmp_theta)
	# print('len_dx',len(d_funx))
	# print('len_theta', len(theta))
	# print('theta[0]:',theta[0])

	# 第二车道线
	fun_x2 = []
	fun_y2 = []
	for i in range(len(fun_x)):
		x2 = fun_x[i] +rightLeft* LaneWidth * math.cos(theta[i] + math.pi / 2.0)
		y2 = fun_y[i] +rightLeft* LaneWidth * math.sin(theta[i] + math.pi / 2.0)
		fun_x2.append(x2)
		fun_y2.append(y2)

	# 第三车道线
	fun_x3 = []
	fun_y3 = []
	for i in range(len(fun_x)):
		x3 = fun_x[i] + rightLeft*2 * LaneWidth * math.cos(theta[i] + math.pi / 2.0)
		y3 = fun_y[i] + rightLeft*2 * LaneWidth * math.sin(theta[i] + math.pi / 2.0)
		fun_x3.append(x3)
		fun_y3.append(y3)

	# refLine
	refx = []
	refy = []
	for i in range(len(fun_x)):
		x = fun_x[i] + rightLeft * 0.5 * LaneWidth * math.cos(theta[i] + math.pi / 2.0)
		y = fun_y[i] + rightLeft * 0.5 * LaneWidth * math.sin(theta[i] + math.pi / 2.0)
		refx.append(x)
		refy.append(y)
	# print('refPosition:',[refx[0], refy[0]])
	if showLane:
		plt.plot(fun_x, fun_y, c='k', lw=0.6, alpha=1)
		plt.plot(fun_x2, fun_y2, c='k', lw=0.4, alpha=1, ls='--')
		plt.plot(fun_x3, fun_y3, c='k', lw=0.6, alpha=1, ls='-')
		plt.plot(refx, refy, c='green', lw=0.3, alpha=0.8, ls='--')

	return b_vector, c_vector


if __name__ == '__main__':
	pass
