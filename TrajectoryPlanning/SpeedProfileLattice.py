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
"""

import matplotlib.pyplot as plt
import TrajectoryPlanning.SpeedProfileGenerator as speed_generator
import numpy as np


max_t = 8   # 8s
max_s = 100	    # 采样高度，100m
max_v = 20
reso_s = 4
reso_t = 1

inc_s = max_v*reso_t    # s方向最大增量
num_s = int(inc_s/reso_s+1)     # s方向子节点数量
num_t = int(max_t / reso_t)

start = [0, 0]
init_v = 0      # init_speed
init_acc = 0    # init_acceleration
init_jerk = 0   # init_jerk
init_s = 0
init_t = 0

max_acc = 2.5      # liu chang liu. 2017 IV. speed profile planning
ref_speed = 12

static_obs1 = [3, 30]
static_obs2 = [4, 30]
static_obs3 = [5, 30]


def child_node(node):
	children = []
	for i in range(num_s):
		child = [node[0] + reso_t, node[1] + i * reso_s]
		if child[0] > max_t:
			return children
		if child[1] > max_s:
			continue
		else:
			children.append(child)
			#   print("child:", child)

	# print("children's length:", len(children))
	return children


def num_node():
	"""
	:return: lattice nodes total number,包括第一个节点，即起始节点也算进去了
	"""
	total = 0
	num_list = []
	for i in range(num_t+1):
		tmp_num = i*(num_s-1) + 1
		if (tmp_num-1)*reso_s <= max_s:
			total = total + tmp_num
			num_list.append(tmp_num)
		else:
			total = total + int(max_s / reso_s) + 1
			num_list.append(int(max_s / reso_s) + 1)

	return total, num_list


def lattice():
	total_node, num_list = num_node()
	print("total:", total_node)
	for k in range(num_t):
		for i in range(num_list[k]):
			tmp_start = [k*reso_t, i * reso_s]
			children = child_node(tmp_start)
			for j in range(len(children)):
				speed_generator.speedprofile(tmp_start, children[j])


if __name__ == '__main__':
	plt.figure(figsize=(3.5, 3.5*0.8))  # 单位英寸， 3.5
	plt.axes([0.18, 0.16, 0.75, 0.76])
	plt.grid(linestyle="--", linewidth=0.3, alpha=1)
	
	lattice()

	font1 = {'family': 'Times New Roman', 'weight': 'normal', 'size': 10}
	plt.xlabel("Time (s)", font1)
	plt.ylabel("s (m)", font1)
	
	# 设置坐标刻度值的大小以及刻度值的字体
	plt.tick_params(labelsize=10)
	x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
	# x = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
	y = np.array([0.0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
	# y = np.array([-2.0, -1, 0.0, 1.0, 2.0])
	
	xgroup_labels = ['0.0', '', '2.0', '', '4.0', '', '6.0', '', '8.0']  # x轴刻度的标识
	ygroup_labels = ['0.0', '', '20.0', '', '40.0', '', '60.0', '', '80.0', '', '100.0']  # y轴刻度的标识
	
	plt.xticks(x, xgroup_labels, fontproperties='Times New Roman', fontsize=10)  # 默认字体大小为10
	plt.yticks(y, ygroup_labels, fontproperties='Times New Roman', fontsize=10)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)

	plt.xlim(-0.5, 8.5)
	plt.ylim(-10, 110)
	
	# plt.legend(loc='upper right')
	# plt.legend()
	# 将文件保存至文件中并且画出图
	plt.savefig('../SimGraph/speedLattice01.svg')
	plt.show()


