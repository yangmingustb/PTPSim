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
=========================
计算ST图中边的权重,速度
"""

import math
import numpy

max_t = 8  # 8s
max_s = 100  # 采样里程，100m
reso_t = 0.5
reso_s = 2
max_v = 20
inc_s = max_v * reso_t  # s方向最大增量
n_s = int(inc_s / reso_s + 1)  # s方向子节点数量
max_ns = int(max_s / reso_s) + 1  # 每一列最多的采样节点的个数

start = [0, 0]
init_v = 0  # init_speed
init_acc = 0  # init_acceleration
init_jerk = 0  # init_jerk
init_s = 0
init_t = 0

max_acc = 2.5  # liu chang liu. 2017 IV. speed profile planning
ref_speed = 12

# static_obs = numpy.array([[3, 3.5, 4, 4.5, 5, 5, 5.5, 6, 4, 4.5], [30, 30, 30, 30, 30.0, 60, 60, 60, 60, 60]])
static_obs=[]
r_circle = 1.0
d_circle = 2.0
obs_inflation = 1.0
# static_obs4 = [6, 30]

alpha1 = 0.50
alpha2 = 0.1
alpha3 = 0.1
alpha4 = 150
alpha5 = 0


def reference_speed(node, closedset):
	pind = node.pind
	father_node = closedset[pind]

	v = (node.y - father_node.y) / reso_t
	g = (v - ref_speed) ** 2
	return g


def acc_weight(node, closedset):
	pind = node.pind
	father_node = closedset[pind]
	v = (node.y - father_node.y) / reso_t

	pind2 = father_node.pind
	if pind2 >= 0:
		grandpa_node = closedset[pind2]

		acc = (node.y - 2 * father_node.y + grandpa_node.y) / reso_t ** 2
	else:
		acc = (v - init_v) / reso_t

	g = acc ** 2
	return g


def jerk_weight(node, closedset):
	pind = node.pind
	father_node = closedset[pind]
	v = (node.y - father_node.y) / reso_t

	pind2 = father_node.pind
	if pind2 >= 0:
		grandpa_node = closedset[pind2]

		acc = (node.y - 2 * father_node.y + grandpa_node.y) / reso_t ** 2

		pind3 = grandpa_node.pind
		if pind3 >= 0:
			grandgrandpa_node = closedset[pind3]
			jerk = (node.y - 3 * father_node.y + 3 * grandpa_node.y - grandgrandpa_node.y) / reso_t ** 3

		else:
			v2 = (father_node.y - init_s) / reso_t

			father_a = (v2 - init_v) / reso_t
			jerk = (acc - father_a) / reso_t

	else:
		acc = (v - init_v) / reso_t
		jerk = (acc - init_acc) / reso_t
	g = jerk ** 2
	return g


def static_obstacles_risk(node):
	d = math.inf
	if static_obs:
		for i in range(static_obs.shape[1]):
			dt = node.x - static_obs[0][i]
			ds = node.y - static_obs[1][i]
			tmp_d = math.sqrt(dt ** 2 + ds ** 2)
			if tmp_d < d:
				d = tmp_d

		# 如果距离等于0,需要加一个小值使得代价值为一个较大的数，而不是无穷大
		sigma = 0.00001
		g = 1 / (d + sigma)
		return g
	else:return 0


def total_cost(node, closedset):
	cost1 = reference_speed(node, closedset)
	cost2 = acc_weight(node, closedset)
	cost3 = jerk_weight(node, closedset)
	cost4 = static_obstacles_risk(node)
	cost5 = 0

	print("cost1:", cost1)
	print("cost2:", cost2)
	print("cost3:", cost3)
	print("cost4:", cost4)

	#   alpha1 = 550

	cost = alpha1 * cost1 + alpha2 * cost2 + alpha3 * cost3 + alpha4 * cost4 + alpha5 * cost5
	print("cost:", cost)
	print("==============================================")
	return cost


if __name__ == '__main__':
	pass
