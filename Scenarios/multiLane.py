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
# @File  : multiLane.py
# @Author: ming.ustb@outlook.com
# @Date  : 19-5-18
# @GitHub: https://github.com/yangmingustb/PTPSim
==============================
这里以最右侧车道线为frenet坐标系
"""

from scipy.linalg import solve
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import FrenetMath.calHeading as heading
from scipy import integrate

p0 = [0.0, 0.0, 70.0 * math.pi / 180.0]
p1 = [50.0, 100.0, 40.0 * math.pi / 180.0]
p2 = [100.0, 120.0, 0.0 * math.pi / 180.0]
p3 = [150.0, 100.0, -30.0 * math.pi / 180.0]
p = [p0, p1, p2, p3]
segment_num = len(p) - 1  # 路径分为几段
LaneWidth = 3.75  # [m]


def wayPoint(x_init, x_goal):
	"""
	道路点离散化
	:param x_init:
	:param x_goal:
	:return:输出的步长为1m
	"""
	x0 = x_init[0]
	y0 = x_init[1]
	theta0 = x_init[2]

	xg = x_goal[0]
	yg = x_goal[1]
	thetag = x_goal[2]

	A = np.array([[1, x0, x0 ** 2, x0 ** 3],
	              [1, xg, xg ** 2, xg ** 3],
	              [0, 1, 2 * x0, 3 * x0 ** 2],
	              [0, 1, 2 * xg, 3 * xg ** 2]])

	b = np.array([y0, yg, np.tan(theta0), np.tan(thetag)])

	a = solve(A, b)
	#   print(a)
	x = np.arange(x0, xg, 0.1)
	y = a[0] + a[1] * x + a[2] * x ** 2 + a[3] * x ** 3
	# plt.plot(x, y, c='black', lw=0.5)
	theta = []
	for i in range(len(x) - 1):
		state0 = [x[i], y[i]]
		state1 = [x[i + 1], y[i + 1]]
		tmp_theta = calHeading(state0, state1)
		theta.append(tmp_theta)
	theta.append(thetag)

	x_list = []
	y_list = []
	theta_list = []
	for i in range(len(x)):
		if 20 * i > (len(x) - 1):
			break
		x_list.append(x[20 * i])
		y_list.append(y[20 * i])
		theta_list.append(theta[20 * i])
	x_list.append(x[-1])
	y_list.append(y[-1])
	theta_list.append(thetag)
	return x_list, y_list, theta_list


def calHeading(state0, state1):
	x0 = state0[0]
	y0 = state0[1]
	x1 = state1[0]
	y1 = state1[1]

	theta = math.atan((y1 - y0) / (x1 - x0))
	return theta


def arcLengthCurve(x_init, x_goal, s0):
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

	A = np.array([[1, x0, x0 ** 2, x0 ** 3],
	              [1, xg, xg ** 2, xg ** 3],
	              [0, 1, 2 * x0, 3 * x0 ** 2],
	              [0, 1, 2 * xg, 3 * xg ** 2]])

	b = np.array([y0, yg, np.tan(theta0), np.tan(thetag)])

	a = solve(A, b)
	#   print(a)

	# x = np.arange(x0, xg, 0.1)
	# y = a[0] + a[1] * x + a[2] * x ** 2 + a[3] * x ** 3
	# plt.plot(x, y, c='r', lw=0.5)

	quad_fun = lambda x: math.sqrt((a[1] + 2 * a[2] * x + 3 * a[3] * x ** 2) ** 2 + 1)

	# 积分
	result = integrate.quad(quad_fun, x0, xg)
	sf = result[0]
	sf = s0 + sf
	# print(sf)

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
	plt.plot(fun_x, fun_y, c='k', lw=0.6, alpha=1)
	theta = []
	for i in range(len(d_funx)):
		tmp_theta = math.atan((d_funy[i]) / (d_funx[i]))
		theta.append(tmp_theta)
	# print('len_dx',len(d_funx))
	# print('len_theta', len(theta))
	# print('theta[0]:',theta[0])

	# 第二车道线
	fun_x2 = []
	fun_y2 = []
	for i in range(len(fun_x)):
		x2 = fun_x[i] + LaneWidth * math.cos(theta[i] + math.pi / 2.0)
		y2 = fun_y[i] + LaneWidth * math.sin(theta[i] + math.pi / 2.0)
		fun_x2.append(x2)
		fun_y2.append(y2)
	plt.plot(fun_x2, fun_y2, c='k', lw=0.4, alpha=1, ls='--')

	# 第三车道线
	fun_x3 = []
	fun_y3 = []
	for i in range(len(fun_x)):
		x3 = fun_x[i] + 2 * LaneWidth * math.cos(theta[i] + math.pi / 2.0)
		y3 = fun_y[i] + 2 * LaneWidth * math.sin(theta[i] + math.pi / 2.0)
		fun_x3.append(x3)
		fun_y3.append(y3)
	plt.plot(fun_x3, fun_y3, c='k', lw=0.4, alpha=1, ls='--')

	# 第四车道线
	fun_x4 = []
	fun_y4 = []
	for i in range(len(fun_x)):
		x4 = fun_x[i] + 3 * LaneWidth * math.cos(theta[i] + math.pi / 2.0)
		y4 = fun_y[i] + 3 * LaneWidth * math.sin(theta[i] + math.pi / 2.0)
		fun_x4.append(x4)
		fun_y4.append(y4)
	plt.plot(fun_x4, fun_y4, c='k', lw=0.6, alpha=1)

	# refLine
	refx = []
	refy = []
	for i in range(len(fun_x)):
		x = fun_x[i] + 1.5 * LaneWidth * math.cos(theta[i] + math.pi / 2.0)
		y = fun_y[i] + 1.5 * LaneWidth * math.sin(theta[i] + math.pi / 2.0)
		refx.append(x)
		refy.append(y)
	# print('refPosition:',[refx[0], refy[0]])
	plt.plot(refx, refy, c='green', lw=0.3, alpha=0.8)

	return sf, b_vector, c_vector


def curvePath():
	s0 = 0.0
	for i in range(segment_num):
		x, y, theta = wayPoint(p[i], p[i + 1])
		for i in range(len(x) - 1):
			start = [x[i], y[i], theta[i]]
			end = [x[i + 1], y[i + 1], theta[i + 1]]
			s0, b, c = arcLengthCurve(start, end, s0)


# ==================================================

def frenetToXY(s, rho, thetaRho, efficients):
	xr, yr, thetar = findEfficients(s, efficients)
	x = xr + rho * math.cos(thetar + math.pi / 2.0)
	y = yr + rho * math.sin(thetar + math.pi / 2.0)
	theta = thetar + thetaRho
	return x, y, theta


def findEfficients(s, efficients):

	s_id = 0
	for i in range(len(efficients)):
		if (s > efficients[i][0]):
			s_id = i
			continue
		else:
			break

	s0, b_vector, c_vector = efficients[s_id]

	x = b_vector[0] + b_vector[1] * s + b_vector[2] * s ** 2 + b_vector[3] * s ** 3
	d_x = b_vector[1] + 2 * b_vector[2] * s + 3 * b_vector[3] * s ** 2

	y = c_vector[0] + c_vector[1] * s + c_vector[2] * s ** 2 + c_vector[3] * s ** 3
	d_y = c_vector[1] + 2 * c_vector[2] * s + 3 * c_vector[3] * s ** 2

	theta = math.atan(d_y / d_x)
	return x,y,theta


def saveEfficients():
	s0 = 0.0
	efficients = []

	for i in range(segment_num):
		x, y, theta = wayPoint(p[i], p[i + 1])
		for i in range(len(x) - 1):
			start = [x[i], y[i], theta[i]]
			end = [x[i + 1], y[i + 1], theta[i + 1]]
			tmp_efficients = [s0]
			s0, b, c = calEfficients(start, end, s0)
			tmp_efficients.append(b)
			tmp_efficients.append(c)
			efficients.append(tmp_efficients)
	return efficients


def calEfficients(x_init, x_goal, s0):
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

	A = np.array([[1, x0, x0 ** 2, x0 ** 3],
	              [1, xg, xg ** 2, xg ** 3],
	              [0, 1, 2 * x0, 3 * x0 ** 2],
	              [0, 1, 2 * xg, 3 * xg ** 2]])

	b = np.array([y0, yg, np.tan(theta0), np.tan(thetag)])

	a = solve(A, b)
	#   print(a)

	quad_fun = lambda x: math.sqrt((a[1] + 2 * a[2] * x + 3 * a[3] * x ** 2) ** 2 + 1)

	# 积分
	result = integrate.quad(quad_fun, x0, xg)
	sf = result[0]
	sf = s0 + sf
	# print(sf)

	B = np.array([[1, s0, s0 ** 2, s0 ** 3],
	              [1, sf, sf ** 2, sf ** 3],
	              [0, 1, 2 * s0, 3 * s0 ** 2],
	              [0, 1, 2 * sf, 3 * sf ** 2]])

	b_b = np.array([x0, xg, np.cos(theta0), np.cos(thetag)])

	b_vector = solve(B, b_b)
	# print(a)
	# 最右侧车道线
	c_b = np.array([y0, yg, np.sin(theta0), np.sin(thetag)])
	c_vector = solve(B, c_b)
	return sf, b_vector, c_vector


if __name__ == '__main__':
	plt.figure(figsize=(3.5, 3.5 * 0.9))  # 单位英寸， 3.5
	# plt.subplot(111, axisbg='#FFDAB9')  # 在图标1中创建子图
	# plt.subplot(111, facecolor='black')  # 在图标1中创建子图
	# plt.style.use('seaborn-dark')
	plt.axes([0.15, 0.15, 0.75, 0.75])
	plt.axis("equal")
	plt.grid(linestyle="--", linewidth=0.2, alpha=1)

	# x2, y2, sf = Polynomial(p1, p2, sf)

	# x3, y3, sf = Polynomial(p2, p3, sf)
	font1 = {'family': 'Times New Roman',
	         'weight': 'normal',
	         'size': 10,
	         }
	plt.xlabel('x (m)', font1)
	plt.ylabel('y (m)', font1)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)

	curvePath()

	plt.savefig('/home/ming/桌面/PTPSim/SimGraph/multiLane.svg')
	plt.show()

	efficients = saveEfficients()
	x, y, theta = frenetToXY(p0[0], p0[1], p0[2],  efficients)
	print('x,y:', x,y)
