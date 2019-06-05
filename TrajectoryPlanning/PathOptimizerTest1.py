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
@File  : PathOptimizerTest1.py
@Author: ming.ustb@outlook.com
@Date  : 19-6-5
@GitHub: https://github.com/yangmingustb/PTPSim
"""


from casadi import *
import numpy as np
import math
import matplotlib.pyplot as plt
import model.simModel as car
import TrajectoryPlanning.PathPlanner as pp
import model.simModel as car


show_animation = 1
s_max = 100.0
reso_s = 1.0
n_s = int(s_max / reso_s)
r_circle = 1.0  # 车体包络圆
d_circle = 2.0  # 圆心间距
obs_inflation = 1.0
safe_distance = obs_inflation + r_circle
kappa_max = 0.187
w_d = 1.0
w_dd = 10.0
w_ddd = 50.0
w_ref = 0.1

lane_width = 3.75
# refLineRho = lane_width*0.5
# 障碍物表示
# static_obs = [[20, refLineRho-2], [40, refLineRho - 1], [70, refLineRho-2]]  # 障碍物的frenet坐标,
rho_max = lane_width  # 行驶通道的左右边界
rho_min = -lane_width
refLineRho = -lane_width*0.5
static_obs = [[20, -2], [40, 2], [70, -2]]
rho_r = [refLineRho for i in range(n_s + 3)]
# 初值
rho_init_guess = refLineRho * np.ones((1, n_s + 3))
rho_min_list = [rho_min for i in range(n_s+3)]
rho_max_list = [rho_max for i in range(n_s+3)]

# 起点状态
x0 = 0
y0 = refLineRho
s0 = 0
rho0 = refLineRho
theta0 = 0
kappa0 = 0
x_init = [x0, y0, s0, rho0, theta0, kappa0]

# 里程s的序列
s = [i for i in np.arange(reso_s, (s_max + 4 * reso_s), reso_s)]
# print("len(s)", len(s))


def RefLine(s):
	xr = s
	yr = 0.0
	thetar = 0.0 * math.pi / 180.0

	return xr, yr, thetar


def frenetToXY(s, rho):
	xr, yr, thetar = RefLine(s)
	x = xr + rho * math.cos(thetar + math.pi / 2.0)
	y = yr + rho * math.sin(thetar + math.pi / 2.0)

	return x, y


def boundValue():
	# 等式约束
	xr0, yr0, thetar0 = RefLine(s0)
	rho_1 = reso_s * (math.tan(theta0 - thetar0)) + rho0

	x1, y1 = frenetToXY(s0 + reso_s, rho_1)
	s2 = s0 + 2 * reso_s
	xr2, yr2, thetar2 = RefLine(s2)
	m = kappa0 * ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** (3 / 2)
	n = (x1 - x0) * (yr2 - 2 * y1 + y0) + (y1 - y0) * (xr2 - 2 * x1 + x0)
	q = (x1 - x0) * math.sin(thetar2 + math.pi / 2.0) - (y1 - y0) * math.cos(thetar2 + math.pi / 2.0)
	rho_2 = (m - n) / q
	return rho_1, rho_2


# Objective
def objective(rho):
	f_d = 0
	f_dd = 0
	f_ddd = 0
	f_ref = 0
	for i in range(n_s):
		f_d = f_d + (rho[i + 1] - rho[i]) ** 2 / (reso_s ** 2)
		f_dd = f_dd + (rho[i + 2] - 2 * rho[i + 1] + rho[i]) ** 2 / (reso_s ** 4)
		f_ddd = f_ddd + (rho[i + 3] - 3 * rho[i + 2] + 3 * rho[i + 1] - rho[i]) ** 2 / reso_s ** 6
		f_ref = f_ref + (rho[i] - rho_r[i]) ** 2

	# print('=========================================')
	# print(f_d)
	# print(f_dd)
	# print(f_ddd)
	# print(f_ref)

	f_d = reso_s * w_d * f_d
	f_dd = reso_s * w_dd * f_dd
	f_ddd = reso_s * w_ddd * f_ddd
	f_ref = reso_s * w_ref * f_ref
	f = f_d + f_dd + f_ddd + f_ref
	return f


def inequality_cons(rho, s, optimization):
	# 曲率约束
	for i in range(n_s):
		x0, y0 = frenetToXY(s[i], rho[i])
		x1, y1 = frenetToXY(s[i + 1], rho[i + 1])
		x2, y2 = frenetToXY(s[i + 2], rho[i + 2])
		k1 = (x1 - x0) * (y2 - 2 * y1 + y0)
		k2 = (y1 - y0) * (x2 - 2 * x1 + x0)
		k3 = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** (3.0 / 2.0)
		kappa1 = -k1 + k2 + k3 * kappa_max
		kappa2 = k1 - k2 + k3 * kappa_max
		optimization.subject_to(kappa1 >= 0)
		optimization.subject_to(kappa2 >= 0)

	# 障碍物约束
	for i in range(n_s + 3):
		# the first obs
		if ((s[i] - static_obs[0][0]) ** 2 < 9):

			c1_obs = (s[i] - static_obs[0][0]) ** 2 + (rho[i] - static_obs[0][1]) ** 2
			c2_obs = (s[i] + d_circle - static_obs[0][0]) ** 2 + (rho[i] - static_obs[0][1]) ** 2
			c3_obs = (s[i] + 2 * d_circle - static_obs[0][0]) ** 2 + (rho[i] - static_obs[0][1]) ** 2
			optimization.subject_to(c1_obs >= safe_distance ** 2)
			optimization.subject_to(c2_obs >= safe_distance ** 2)
			optimization.subject_to(c3_obs >= safe_distance ** 2)
		# the second obs
		if ((s[i] - static_obs[1][0]) ** 2 < 9):
			c1_obs2 = (s[i] - static_obs[1][0]) ** 2 + (rho[i] - static_obs[1][1]) ** 2
			c2_obs2 = (s[i] + d_circle - static_obs[1][0]) ** 2 + (rho[i] - static_obs[1][1]) ** 2
			c3_obs2 = (s[i] + 2 * d_circle - static_obs[1][0]) ** 2 + (rho[i] - static_obs[1][1]) ** 2

			optimization.subject_to(c1_obs2 >= safe_distance ** 2)
			optimization.subject_to(c2_obs2 >= safe_distance ** 2)
			optimization.subject_to(c3_obs2 >= safe_distance ** 2)

		# the second obs
		if ((s[i] - static_obs[2][0]) ** 2 < 9):
			c1_obs3 = (s[i] - static_obs[2][0]) ** 2 + (rho[i] - static_obs[2][1]) ** 2
			c2_obs3 = (s[i] + d_circle - static_obs[2][0]) ** 2 + (rho[i] - static_obs[2][1]) ** 2
			c3_obs3 = (s[i] + 2 * d_circle - static_obs[2][0]) ** 2 + (rho[i] - static_obs[2][1]) ** 2

			optimization.subject_to(c1_obs3 >= safe_distance ** 2)
			optimization.subject_to(c2_obs3 >= safe_distance ** 2)
			optimization.subject_to(c3_obs3 >= safe_distance ** 2)


def allConstraints(rho, s, optimization):
	# 曲率和避障约束
	inequality_cons(rho, s, optimization)
	rho1, rho2 = boundValue()
	# 边值条件约束
	optimization.subject_to(rho[0] == rho1)
	optimization.subject_to(rho[1] == rho2)

	# 取值范围约束
	optimization.subject_to(optimization.bounded(rho_min_list, rho, rho_max_list))


def ipoptSolver():
	# construct a optimization problem
	optimization = casadi.Opti()

	# column vector, decision variables
	rho = optimization.variable(n_s + 3)  # number of variables:n_s+3
	# print("rho_shape[0]:", rho.shape[0])

	f = objective(rho)
	optimization.minimize(f)

	allConstraints(rho, s, optimization)

	p_opts = {"expand": True}
	s_opts = {"max_iter": 100}
	optimization.solver("ipopt", p_opts, s_opts)
	# 设置初值
	optimization.set_initial(rho, rho_init_guess)
	solution = optimization.solve()
	rho_vector = solution.value(rho)
	# print("rho_vector:", rho_vector)
	return  rho_vector


def plotGraph():
	# plot graph
	plt.figure(figsize=(3.5, 3.0))  # 单位英寸， 3.5
	font1 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 10}
	plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 如果要显示中文字体，则在此处设为：SimHei

	p1 = [0.2, 0.2, 0.7, 0.6]
	plt.axes(p1)

	rho_vector = ipoptSolver()
	# plot trajectory
	plt.plot(s, rho_vector, c='cyan', linestyle="-", linewidth=1, alpha=1, label='Refined path')

	# 调用pathplanner的环境
	pp.plotGraph()


	car.simVehicle([x_init[2], x_init[3]], 0 * math.pi / 180, 'g', 1)

	# plot obstacles
	car.simVehicle([static_obs[0][0], static_obs[0][1]], 0 * math.pi / 180, 'r', 1)
	car.simVehicle([static_obs[1][0], static_obs[1][1]], 0 * math.pi / 180, 'r', 1)
	car.simVehicle([static_obs[2][0], static_obs[2][1]], 0 * math.pi / 180, 'r', 1)

	plt.grid(linestyle="--", linewidth=0.5, alpha=1)
	plt.title('x-y Graph', font1)
	plt.xlabel('x (m)', font1)
	plt.ylabel('y (m)', font1)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)
	# plt.xlim(-1, 110)
	# plt.ylim(-4, 4)
	plt.legend(loc=0)  # 图例位置自动
	plt.savefig('../SimGraph/pathOptimizationPath060414.svg')

	plt.figure(figsize=(3.5, 3.5))  # 单位英寸， 3.5
	font1 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 10}
	plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 如果要显示中文字体，则在此处设为：SimHei

	p2 = [0.2, 0.2, 0.7, 0.6]
	plt.axes(p2)
	kappa_list = trajectory_kappa(rho_vector, s)
	plt.plot(s[0:n_s + 1], kappa_list, c= 'cyan', linestyle="-", linewidth=0.8, alpha=1, label='Refined Path')

	# plot rough path kappa
	pp.plotkappaGraph()

	y = [kappa_max for i in range(n_s + 1)]
	y2 = [-kappa_max for i in range(n_s + 1)]
	plt.plot(s[0:n_s + 1], y, c= 'r', linestyle="--", linewidth=0.5, alpha=1, label='Maximum Curvature Limit')
	plt.plot(s[0:n_s + 1], y2, c= 'orange', linestyle="--", linewidth=0.5, alpha=1, label='Minimum Curvature Limit')
	plt.title('Curvature Profile', font1)
	plt.grid(linestyle="--", linewidth=0.5, alpha=1)
	plt.xlabel('s (m)', font1)
	plt.ylabel('kappa (1/m)', font1)
	plt.xlim(-1, 110)
	plt.ylim(-0.2, 0.2)
	plt.legend(loc=0)  # 图例位置自动
	plt.savefig('../SimGraph/pathOptimizationKappa060413.svg')
	plt.show()


def trajectory_kappa(rho, s):
	tmp_kappa = []

	for i in range(len(rho)-2):
		x0, x1, x2 = s[i], s[i+1], s[i+2]
		y0, y1, y2 = rho[i], rho[i+1], rho[i+2]
		k1 = (x1 - x0)*(y2-2*y1+y0)
		k2 = (y1-y0)*(x2-2*x1+x0)
		k3 = ((x1-x0)**2+(y1-y0)**2)**(3.0/2.0)

		if (k3 == 0.0):
			tmp_kappa.append(0)

		else:
			tmp_kappa.append((k1 - k2) / k3)
	return tmp_kappa


if __name__ == '__main__':
	plotGraph()

