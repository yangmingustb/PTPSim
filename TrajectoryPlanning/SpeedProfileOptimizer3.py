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
====================
speed profile optimizer
通过conda安装cvxopt,cvxpy
calculate quadratic programming problems
======================
min     j=1/2(x.T)Hx+qx
s.t.    Gx<h
        Ax=b
=======================

"""

import pprint
from cvxopt import matrix, solvers
import math
import numpy
import matplotlib.pyplot as plt


max_t = 8  # 8s
max_s = 100  # 100m
reso_t = 0.1
reso_s = 2
max_v = 20
max_acc = 2.5  # liu chang liu. 2017 IV. speed profile planning
ref_v = 40.0/3.6  # 前方没有障碍物时的参考速度
#   ref_v_static = 0  #   前方有静态障碍物时的速度
max_jerk = 5
n_t = int(max_t / reso_t)
lane_width = 3.75
refLineRho = lane_width*0.5
# 起点状态
s0 = 100.0
s_end = s0+max_s
t0 = 12.4
start_SRho = [s0, refLineRho, 0.0 * math.pi / 180.0]
start = [s0, t0]
init_v = 30/3.6  # init_speed
init_acc = 0  # init_acceleration
init_jerk = 0  # init_jerk
init_s = s0
init_t = t0

r = 1  # 碰撞检测圆半径
d = 2  # 圆心间距
# obs = [[60, 0]]
obs = []
obs_threshold = 2  # 静态障碍物膨胀阈值
safeDistance = r + obs_threshold
safeSpeed = 30/3.6      # 安全速度
overtakingSpeed = 60/3.6
flag_obs = 1  # 障碍物标志位，0前方无障碍物，1前方静态障碍物，2前方动态障碍物
rough_st = []
rough_v = []
guidance_v = ref_v * numpy.ones((1, n_t))
for i in range(n_t):
	if i == int(1/max_s * n_t):
		# print(i)
		guidance_v[0][i] = overtakingSpeed
		guidance_v[0][i+1] = overtakingSpeed
		guidance_v[0][i+2] = overtakingSpeed
		guidance_v[0][i+3] = overtakingSpeed
		guidance_v[0][i+4] = overtakingSpeed
	elif i == int(40/max_s * n_t):
		guidance_v[0][i] = ref_v
		guidance_v[0][i+1] = ref_v
		guidance_v[0][i+2] = ref_v
		guidance_v[0][i+3] = ref_v
		guidance_v[0][i+4] = ref_v
	elif i == int(70/max_s * n_t):
		guidance_v[0][i] = safeSpeed
		guidance_v[0][i+1] = safeSpeed
		guidance_v[0][i+2] = safeSpeed
		guidance_v[0][i+3] = safeSpeed
		guidance_v[0][i+4] = safeSpeed
	else:continue
print("----------------")
print(guidance_v)
w_vel = 100
w_acc = 1
w_jerk = 10
showAnimation = 1


def cal_H_vel():
	H_vel = numpy.eye(n_t + 1)

	for i in range(n_t + 1):
		if i > 0 and i < n_t:
			H_vel[i][i] = 2
		if i > 0 and i < n_t + 1:
			H_vel[i - 1][i] = -1
			H_vel[i][i - 1] = -1

	#   print("H_vel size:", numpy.shape(H_vel), H_vel)

	return H_vel


def cal_q_vel():
	"""
	暂时使用默认参考速度进行优化，等调试成熟再使用粗糙速度来优化
	:return:
	"""
	q_vel = numpy.zeros((1, n_t + 1))

	if flag_obs == 0:
		q_vel[0][0] = -ref_v
		q_vel[0][n_t] = ref_v

	if flag_obs == 1:
		for i in range(n_t + 1):
			if i < 1:
				q_vel[0][i] = -guidance_v[0][i]
			elif i >= n_t:

				q_vel[0][i] = guidance_v[0][i - 1]

			else:
				q_vel[0][i] = guidance_v[0][i - 1] - guidance_v[0][i]
	#   print('q_vel:', numpy.shape(q_vel), q_vel)

	return q_vel


def cal_H_acc():
	H_acc = numpy.eye(n_t + 2)

	for i in range(n_t + 2):
		if i == 1:
			H_acc[i][i] = 5
			H_acc[i - 1][i] = - 2
			H_acc[i][i - 1] = - 2
		if i > 1 and i < n_t:
			H_acc[i][i] = 6
			H_acc[i - 1][i] = - 4
			H_acc[i - 2][i] = 1
			H_acc[i][i - 1] = - 4
			H_acc[i][i - 2] = 1

		if i == n_t:
			H_acc[i][i] = 5
			H_acc[i - 1][i] = - 4
			H_acc[i - 2][i] = 1
			H_acc[i][i - 1] = - 4
			H_acc[i][i - 2] = 1

		if i == n_t + 1:
			H_acc[i - 1][i] = - 2
			H_acc[i - 2][i] = 1
			H_acc[i][i - 1] = - 2
			H_acc[i][i - 2] = 1

	#   print("H_acc size:", numpy.shape(H_acc), "\n", H_acc)

	return H_acc


def cal_H_jerk():
	H_jerk = numpy.eye(n_t + 3)

	for i in range(n_t + 3):
		if i == 1:
			H_jerk[i][i] = 10
			H_jerk[i - 1][i] = - 3
			H_jerk[i][i - 1] = - 3
		if i == 2:
			H_jerk[i][i] = 19
			H_jerk[i - 1][i] = - 12
			H_jerk[i - 2][i] = 3
			H_jerk[i][i - 1] = - 12
			H_jerk[i][i - 2] = 3

		if i > 2 and i < n_t:
			H_jerk[i][i] = 20
			H_jerk[i - 1][i] = - 15
			H_jerk[i - 2][i] = 6
			H_jerk[i - 3][i] = - 1

			H_jerk[i][i - 1] = - 15
			H_jerk[i][i - 2] = 6
			H_jerk[i][i - 3] = - 1

		if i == n_t:
			H_jerk[i][i] = 19
			H_jerk[i - 1][i] = - 15
			H_jerk[i - 2][i] = 6
			H_jerk[i - 3][i] = - 1

			H_jerk[i][i - 1] = - 15
			H_jerk[i][i - 2] = 6
			H_jerk[i][i - 3] = - 1

		if i == n_t + 1:
			H_jerk[i][i] = 10
			H_jerk[i - 1][i] = - 12
			H_jerk[i - 2][i] = 6
			H_jerk[i - 3][i] = - 1
			H_jerk[i][i - 1] = - 12
			H_jerk[i][i - 2] = 6
			H_jerk[i][i - 3] = - 1

		if i == n_t + 2:
			H_jerk[i - 1][i] = - 3
			H_jerk[i - 2][i] = 3
			H_jerk[i - 3][i] = - 1
			H_jerk[i][i - 1] = - 3
			H_jerk[i][i - 2] = 3
			H_jerk[i][i - 3] = - 1

	#   print("H_jerk size:", numpy.shape(H_jerk), "\n", H_jerk)
	#   print(H_jerk[82])
	return H_jerk


def cal_H_total():
	H_vel = cal_H_vel()
	tmp_y = numpy.zeros((n_t + 1, 2))
	H_vel = numpy.c_[H_vel, tmp_y]
	tmp_x = numpy.zeros((2, n_t + 3))
	H_vel = numpy.r_[H_vel, tmp_x]

	H_acc = cal_H_acc()
	tmp_y = numpy.zeros((n_t + 2, 1))
	H_acc = numpy.c_[H_acc, tmp_y]
	tmp_x = numpy.zeros((1, n_t + 3))
	H_acc = numpy.r_[H_acc, tmp_x]

	H_jerk = cal_H_jerk()
	#   print(H_vel)
	#   print(H_acc)
	#   print(H_jerk)

	H = 2 * w_vel / reso_t * H_vel + 2 * w_acc / (reso_t ** 3) * H_acc + 2 * w_jerk / (reso_t ** 5) * H_jerk
	return H


def cal_q():
	q_vel = cal_q_vel()
	tmp_y = numpy.zeros((1, 2))
	tmp_q = numpy.c_[q_vel, tmp_y]
	q = -2 * w_vel * tmp_q

	#   print(q)

	return q


def constraints_acc():
	h_acc = reso_t ** 2 * max_acc * numpy.ones((n_t, 1))
	g_acc = numpy.zeros((n_t, n_t + 3))
	for i in range(n_t):
		g_acc[i][i] = 1
		g_acc[i][i + 1] = -2
		g_acc[i][i + 2] = 1

	return g_acc, h_acc


def constraints_jerk():
	h_jerk = reso_t ** 3 * max_jerk * numpy.ones((n_t, 1))
	g_jerk = numpy.zeros((n_t, n_t + 3))
	for i in range(n_t):
		g_jerk[i][i] = -1
		g_jerk[i][i + 1] = 3
		g_jerk[i][i + 2] = -3
		g_jerk[i][i + 3] = 1

	return g_jerk, h_jerk


def constraints_speed():
	h_vel = reso_t * max_v * numpy.ones((n_t, 1))
	g_vel = numpy.zeros((n_t, n_t + 3))
	for i in range(n_t):
		g_vel[i][i] = -1
		g_vel[i][i + 1] = 1

	return g_vel, h_vel


def constraints_monotonicity():
	"""
	单调性约束
	:return: 维度, (n_t+2)*(n_t+3)
	"""
	h_inc = numpy.zeros((n_t + 2, 1))
	g_inc = numpy.zeros((n_t + 2, n_t + 3))
	for i in range(n_t + 2):
		g_inc[i][i] = 1
		g_inc[i][i + 1] = -1

	return g_inc, h_inc


def constraints_range():
	h_max = s_end * numpy.ones((n_t, 1))
	h_init = init_s * numpy.ones((n_t, 1))
	g_range = numpy.zeros((n_t, n_t + 3))
	for i in range(n_t):
		g_range[i][i] = 1

	return g_range, h_max, h_init


def static_obs_constraints():
	h_obs1 = (obs[0][0] - safeDistance - 2 * d) * numpy.ones((n_t + 3, 1))
	h_obs2 = (-obs[0][0] - safeDistance) * numpy.ones((n_t + 3, 1))
	g_obs1 = numpy.eye(n_t + 3)
	g_obs2 = -numpy.eye(n_t + 3)

	return g_obs1, g_obs2, h_obs1, h_obs2


def constraints_equation():
	b = numpy.zeros((6, 1))
	b[0][0] = reso_t * init_v + init_s
	b[1][0] = init_acc * reso_t ** 2 + init_s + 2 * reso_t * init_v
	b[2][0] = init_jerk * reso_t ** 3 + 3 * reso_t * init_v + 3 * init_acc * reso_t ** 2 + init_s

	a = numpy.zeros((6, n_t + 3))
	a[0][0] = 1
	a[1][1] = 1
	a[2][2] = 1

	a[3][n_t - 2] = 1
	a[3][n_t - 1] = -2
	a[3][n_t] = 1

	a[4][n_t - 2] = 2
	a[4][n_t - 1] = -3
	a[4][n_t + 1] = 1

	a[5][n_t - 2] = 3
	a[5][n_t - 1] = -4
	a[5][n_t + 2] = 1

	return a, b


def all_constraints():
	g_acc, h_acc = constraints_acc()
	g_jerk, h_jerk = constraints_jerk()
	g_vel, h_vel = constraints_speed()
	# g_inc, h_inc = constraints_monotonicity()
	g_max, h_max, h_init = constraints_range()

	g = numpy.r_[g_acc, -g_acc]
	g = numpy.r_[g, g_jerk]
	g = numpy.r_[g, -g_jerk]
	g = numpy.r_[g, g_vel]
	# g = numpy.r_[g, g_inc]
	g = numpy.r_[g, g_max]
	g = numpy.r_[g, -g_max]

	h = numpy.r_[h_acc, h_acc]
	h = numpy.r_[h, h_jerk]
	h = numpy.r_[h, h_jerk]
	h = numpy.r_[h, h_vel]
	# h = numpy.r_[h, h_inc]
	h = numpy.r_[h, h_max]
	h = numpy.r_[h, h_init]

	if len(obs):
		g_obs1, g_obs2, h_obs1, h_obs2 = static_obs_constraints()
		g = numpy.r_[g, g_obs1]
		# g = numpy.r_[g, g_obs2]
		h = numpy.r_[h, h_obs1]
		# h = numpy.r_[h, h_obs2]

	a, b = constraints_equation()
	return g, h, a, b


def QP_solver():
	H = cal_H_total()
	q = cal_q()
	q = q.T
	g, h, a, b = all_constraints()

	P = matrix(H)
	q = matrix(q)
	G = matrix(g)
	h = matrix(h)
	A = matrix(a)
	b = matrix(b)
	result = solvers.qp(P, q, G, h, A, b)
	s_set = result['x']

	return s_set


def plotTrajectory():
	s_set = QP_solver()
	s_set = numpy.array(s_set).T

	#   print(type(s_set))
	shape = numpy.shape(s_set)
	t_tmp = [init_t]
	s_tmp = [init_s]

	for i in range(shape[1]):
		t_tmp.append(init_t + reso_t * (i + 1))
		s_tmp.append(s_set[0][i])
	return t_tmp, s_tmp


def plot_figure():
	s_set = QP_solver()
	s_set = numpy.array(s_set).T

	#   print(type(s_set))
	shape = numpy.shape(s_set)
	t_tmp = [init_t]
	s_tmp = [init_s]

	for i in range(shape[1]):
		t_tmp.append(init_t + reso_t * (i + 1))
		s_tmp.append(s_set[0][i])

	print('s_tmp:', s_tmp)

	font1 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 10}
	plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 如果要显示中文字体，则在此处设为：SimHei

	# plt.style.use('ggplot')
	plt.figure(figsize=(3.5, 4.2))  # 单位英寸， 3.5

	p1 = [0.18, 0.85, 0.72, 0.1]
	p2 = [0.18, 0.60, 0.72, 0.08]
	p3 = [0.18, 0.35, 0.72, 0.05]
	p4 = [0.18, 0.1, 0.72, 0.05]
	ax1 = plt.axes(p1)  # 2行2列,第一个图
	ax2 = plt.axes(p2)  # 2行2列,第一个图
	ax3 = plt.axes(p3)  # 2行2列,第一个图
	ax4 = plt.axes(p4)  # 2行2列,第一个图

	ax1.plot(t_tmp, s_tmp, lw=1, label='s-t Graph', c='b', alpha=1)
	ax1.grid(linestyle="--", linewidth=0.5, alpha=1)
	ax1.set_title('s-t Graph', font1)
	ax1.set_xlabel('Time (s)', font1)
	ax1.set_ylabel('s (m)', font1)
	# plt.xticks(fontproperties='Times New Roman', fontsize=10)
	# plt.yticks(fontproperties='Times New Roman', fontsize=10)
	# ax1.set_xlim(-1, 9)
	# ax1.set_ylim(-1, 110)
	#   plt.grid(True)
	#   plt.title('A simple plot')

	v_tmp = []
	for i in range(len(s_tmp) - 1):
		v = (s_tmp[i + 1] - s_tmp[i]) / reso_t
		v_tmp.append(v)

	ax2.plot(t_tmp[0: (len(t_tmp) - 1)], v_tmp, lw=1, label='Speed Profile', c='b', alpha=0.8)
	ax2.set_title('Speed Profile', font1)
	ax2.grid(linestyle="--", linewidth=0.5, alpha=1)
	# plt.legend(loc=0)  # 图例位置自动
	# ax2.axis('tight')
	ax2.set_xlabel('Time (s)', font1)
	ax2.set_ylabel('v (m/s)', font1)
	# ax2.set_xlim(-1, 9)
	# ax2.set_ylim(0, 15)

	acc_tmp = []
	for i in range(len(s_tmp) - 2):
		acc = (s_tmp[i + 2] - 2 * s_tmp[i + 1] + s_tmp[i]) / reso_t ** 2
		acc_tmp.append(acc)

	ax3.plot(t_tmp[0: (len(t_tmp) - 2)], acc_tmp, lw=1, label='Acceleration Profile', alpha=0.8, c='b')
	ax3.set_title('Acceleration Profile', font1)
	ax3.grid(linestyle="--", linewidth=0.5, alpha=1)

	# ax3.legend(loc=0)  # 图例位置自动
	# ax3.axis('tight')
	ax3.set_xlabel('Time (s)', font1)
	ax3.set_ylabel(r'a (m/s$\mathrm{^2}$)', font1)
	# ax3.set_xlim(-1, 9)
	# ax3.set_ylim(-5, 5)

	jerk_tmp = []
	for i in range(len(s_tmp) - 3):
		jerk = (s_tmp[i + 3] - 3 * s_tmp[i + 2] + 3 * s_tmp[i + 1] - s_tmp[i]) / reso_t ** 3
		jerk_tmp.append(jerk)

	ax4.plot(t_tmp[0: (len(t_tmp) - 3)], jerk_tmp, lw=1, label='Jerk Profile', c='b', alpha=0.8)
	ax4.set_title('Jerk Profile', font1)
	ax4.grid(linestyle="--", linewidth=0.5, alpha=1)

	# plt.legend(loc=0)  # 图例位置自动
	# ax4.axis('tight')
	ax4.set_xlabel('Time (s)', font1)
	ax4.set_ylabel(r'jerk (m/s$\mathrm{^3}$)', font1)
	# ax4.set_xlim(-1, 9)
	# ax4.set_ylim(-10, 10)
	plt.savefig('../SimGraph/speedOptimization3_061101.svg')
	plt.savefig('../SimGraph/speedOptimization3_061101.tiff', dpi=600)
	#   plt.grid(True)
	plt.show()


if __name__ == '__main__':
	if showAnimation:
		plot_figure()
