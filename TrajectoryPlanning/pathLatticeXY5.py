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
@File  : pathLatticeXY5.py
@Author: ming.ustb@outlook.com
@Date  : 19-6-12
@GitHub: https://github.com/yangmingustb/PTPSim
scenario:overtaking
"""


import numpy as np
import Curves.Cubic as cubic
import matplotlib.pyplot as plt
import math
import Curves.cubic_spline as cubicSpline
import FrenetMath.FrenetToCartesian as ftc
import model.simModel as car
# 注意，这里有两个cubic，一个是连接两个状态，一个是用样条连接所有waypoints构成行驶环境


showSamplingPath = True
show_obstacle = True
showVehicleStart = True
longitudinal_num = 4
lateral_num = 9  # 横向采样个数
longitudinal_step = 20.0
lateral_step = 0.5

lane_width = 3.75
obstacleHeading = 0.0


def parameters(s0, offset, refline):
	# s0 = 400.0    # 全局变量，全局调整
	refLineRho = lane_width * 0.5
	laneChaneRefLine = lane_width * 1.5
	start_SRho = [s0, refline-offset, 0.0 * math.pi / 180.0]
	static_obs = [[s0 + 20, refLineRho - 0.3], [s0 + 30, refLineRho + 0.5], [s0 + 85, laneChaneRefLine - 0.3]]
	return start_SRho, static_obs


def sampling(x_row, y_column, lateral_step, longitudinal_step, refline, start_SRho):
	"""

	:param x_row: s采样个数
	:param y_column: lateral采样个数
	:param lateral_step: 采样步长
	:param longitudinal_step: 纵向采样步长
	:return:s-rho坐标系端点
	"""

	end_set = np.empty(shape=[x_row, y_column, 3])

	for i in range(x_row):
		x_i = (i + 1) * longitudinal_step + start_SRho[0]
		for j in range(y_column):
			y_i = (j - lane_width) * lateral_step + refline
			target_point = [x_i, y_i, 0.0 * math.pi / 180.0]

			end_set[i, j] = np.array(target_point)

	return end_set


def generate_lattice(efficients, refline, start_SRho):
	end_set = sampling(longitudinal_num, lateral_num, lateral_step, longitudinal_step, refline, start_SRho)
	#   print(end_set)
	end_size = end_set.shape
	# print("end_size:")
	# print(end_size)

	# 生成车辆起点到第一列采样点的图
	for i in range(end_size[1]):
		s, rho, thetaRho = cubic.Polynomial(start_SRho, end_set[0, i])
		x = []
		y = []
		# print(s)
		for j in range(len(s)):
			tmpX, tmpY, tmpTheta = ftc.frenetToXY(s[j], rho[j], thetaRho[j], efficients)
			x.append(tmpX)
			y.append(tmpY)
		# plt.scatter(end_set[0, i][0], end_set[0, i][1], color='b', s=2, alpha=0.8)
		plt.plot(x, y, c='b', linewidth=0.2, alpha=1.0)

	# 采样点之间的图
	for i in range(end_size[0] - 1):
		for j in range(end_size[1]):
			#   print([i, j])
			for q in range(end_size[1]):
				#   mptg.test_optimize_trajectory(end_set[1, 0], end_set[0, 1])
				s, rho, thetaRho = cubic.Polynomial(end_set[i, q], end_set[i + 1, j])
				x = []
				y = []
				# print(s)
				for q in range(len(s)):
					tmpX, tmpY, tmpTheta = ftc.frenetToXY(s[q], rho[q], thetaRho[q], efficients)
					x.append(tmpX)
					y.append(tmpY)
				# plt.scatter(end_set[i + 1, j][0], end_set[i + 1, j][1], color='b', s=2, alpha=0.8)
				plt.plot(x, y, c='b', linewidth=0.1, alpha=0.80)
	return None


def plotGraph():
	# plt.style.use('ggplot')
	plt.figure(figsize=(3.5, 3.5 * 0.62))  # 单位英寸， 3.5
	plt.axes([0.2, 0.2, 0.7, 0.7])
	plt.axis("equal")
	plt.grid(linestyle="--", linewidth=0.5, alpha=1)

	# plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 如果要显示中文字体，则在此处设为：SimHei
	# plt.rcParams['axes.unicode_minus'] = False  # 显示负号

	# # show multilane
	# multiLane.curvePath()

	# # 计算多车道环境的弧长参数曲线的系数
	# efficients = multiLane.saveEfficients()
	# 计算回环环境的弧长参数曲线的系数
	efficients = cubicSpline.saveEfficients()

	# show sampling path
	s0 = 400
	offset = 0.3
	refLineRho = lane_width * 0.5
	laneChaneRefLine = lane_width * 1.5
	refline = refLineRho
	start_SRho, static_obs=parameters(s0, offset, refline)
	generate_lattice(efficients, refLineRho, start_SRho)
	generate_lattice(efficients, laneChaneRefLine, start_SRho)

	if show_obstacle:
		c = ['lightseagreen', 'orange', 'green']
		for i in range(len(static_obs)):
			tmpx, tmpy, tmptheta = ftc.frenetToXY(static_obs[i][0], static_obs[i][1], obstacleHeading,
			                                      efficients)
			car.simVehicle([tmpx, tmpy], tmptheta, c[i], 1)

	if showVehicleStart:
		tmpx, tmpy, tmptheta = ftc.frenetToXY(start_SRho[0], start_SRho[1], obstacleHeading, efficients)
		car.simVehicle([tmpx, tmpy], tmptheta, 'b', 1)

	s0 = 480
	offset = 0.0
	refline = laneChaneRefLine
	start_SRho, static_obs=parameters(s0, offset, refline)
	generate_lattice(efficients, refLineRho, start_SRho)
	generate_lattice(efficients, laneChaneRefLine, start_SRho)

	font1 = {'family': 'Times New Roman',
	         'weight': 'normal',
	         'size': 10,
	         }
	plt.xlabel("x (m)", font1)
	plt.ylabel('y (m)', font1)

	# 设置坐标刻度值的大小以及刻度值的字体
	plt.tick_params(labelsize=10)
	# x = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0])
	# x = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
	# y = np.array([-2.0, -1, 0.0, 1.0, 2.0])
	# y = np.array([-2.0, -1, 0.0, 1.0, 2.0])

	# xgroup_labels = ['0.0', '20.0', '40.0', '60.0', '80.0', '100.0']  # x轴刻度的标识
	# ygroup_labels = ['-2.0', '-1.0', '0.0', '1.0', '2.0']  # y轴刻度的标识

	# plt.xticks(x, xgroup_labels, fontproperties='Times New Roman', fontsize=10)  # 默认字体大小为10
	# plt.yticks(y, ygroup_labels, fontproperties='Times New Roman', fontsize=10)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)
	plt.xlim(300, 500)
	plt.ylim(60, 190)
	plt.savefig('../SimGraph/pathLatticeXY5_061203.tiff', dpi=600)
	plt.show()


def ToPathPlanner(efficients):
	# show sampling path
	# show sampling path
	s0 = 400
	offset = 0.3
	refLineRho = lane_width * 0.5
	laneChaneRefLine = lane_width * 1.5
	refline = refLineRho
	start_SRho, static_obs=parameters(s0, offset, refline)
	generate_lattice(efficients, refLineRho, start_SRho)
	generate_lattice(efficients, laneChaneRefLine, start_SRho)

	if showVehicleStart:
		tmpx, tmpy, tmptheta = ftc.frenetToXY(start_SRho[0], start_SRho[1], obstacleHeading, efficients)
		car.simVehicle([tmpx, tmpy], tmptheta, 'b', 1)
	if show_obstacle:
		c = ['lightseagreen', 'orange', 'green']
		for i in range(len(static_obs)):
			tmpx, tmpy, tmptheta = ftc.frenetToXY(static_obs[i][0], static_obs[i][1], obstacleHeading,
			                                      efficients)
			car.simVehicle([tmpx, tmpy], tmptheta, c[i], 1)
	s0 = 480
	offset = 0.0
	refline = laneChaneRefLine
	start_SRho, static_obs=parameters(s0, offset, refline)
	generate_lattice(efficients, refLineRho, start_SRho)
	generate_lattice(efficients, laneChaneRefLine, start_SRho)
	return None


if __name__ == '__main__':

	plotGraph()

