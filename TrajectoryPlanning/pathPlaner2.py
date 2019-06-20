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
@File  : pathPlaner2.py
@Author: ming.ustb@outlook.com
@Date  : 19-5-30
@GitHub: https://github.com/yangmingustb/PTPSim
x-y graph
"""


import matplotlib.pyplot as plt
import math
import TrajectoryPlanning.PathCost as CostFunction
import Curves.Cubic as cubic
import model.simModel as car
import Scenarios.multiLane as mlane
import Curves.cubic_spline as cubicSpline
import FrenetMath.FrenetToCartesian as ftc
import TrajectoryPlanning.pathLatticeXY2 as pathlattice


show_vehicle_animation = False
show_sample_point = False
show_lane = True
show_obstacle = True
show_rough_path = True
showPathLattice = True

s_max = 100.0
longi_num = 5
lateral_num = 9  # 横向采样个数
longi_step = 20.0
latera_step = 0.5
lane_width = 3.75
numberUnderRefLine = lateral_num//2    # sampling number under reference line
s0 = 0.0
s_end=s0+s_max
refLineRho = lane_width*0.5
start_SRho = [s0, refLineRho, 0.0 * math.pi / 180.0]
# start_XY = [-5.285, 1.924, 70.0 * math.pi / 180.0]

obstacle = [[20, refLineRho - 1], [40, refLineRho + 2], [70, refLineRho + 2]]  # 障碍物的frenet坐标,
obstacleHeading=0.0 * math.pi / 180.0

last_column_id = [lateral_num * longi_num, lateral_num * (longi_num-1) + 1]  # 最后一列的编号


class Node:

	def __init__(self, x, y, cost, pind):
		"""

		:param x:
		:param y:
		:param cost: 累计最小代价
		:param pind: 指针，指向父节点
		"""
		self.x = x
		self.y = y
		self.cost = cost
		self.pind = pind

	def __str__(self):
		return str(self.x) + "," + str(self.y) + "," + str(self.cost) + "," + str(self.pind)


def dijkstra_planning(start, longitudinal_step, lateral_step, lateral_number, efficients):
	"""

	:param start:
	:param dyn_obs:
	:param longitudinal_step:
	:param lateral_step:
	:param longitudinal_number:
	:param lateral_number:
	:return:
	"""

	nstart = Node(start[0], start[1], 0.0, -1)

	child_node = get_children_node(longitudinal_step, lateral_step)

	open_set, closed_set = dict(), dict()  # build empty dict
	open_set[0] = nstart

	while True:
		if not open_set:
			break

		c_id = min(open_set, key=lambda o: open_set[o].cost)
		current = open_set[c_id]
		print("current", current)

		# show graph
		if show_sample_point:
			plt.plot(current.x, current.y, "oc")

		# Remove the item from the open set
		del open_set[c_id]
		# Add it to the closed set
		closed_set[c_id] = current

		# expand search grid based on motion model
		for i in range(len(child_node)):
			node = Node(current.x + child_node[i][0], child_node[i][1], current.cost, c_id)
			cost = CostFunction.total_cost(current, node, obstacle, refLineRho, efficients)
			node.cost = node.cost + cost
			n_id = calc_index(node, nstart, longitudinal_step, lateral_step, lateral_number)

			if node.x > s_end:
				break

			if n_id in closed_set:
				continue
			# Otherwise if it is already in the open set
			elif n_id in open_set:
				if open_set[n_id].cost > node.cost:
					open_set[n_id].cost = node.cost
					open_set[n_id].pind = c_id
			else:
				open_set[n_id] = node

	return closed_set


def calc_index(node, nstart, longitudinal_step, lateral_step, latera_num):
	"""

	:param node: 子节点
	:param nstart: 起点
	:param longitudinal_step:
	:param lateral_step:
	:param latera_num:横向采样点个数
	:return:
	"""
	id = (node.y - (refLineRho-numberUnderRefLine*lateral_step)) / lateral_step + ((node.x - nstart.x) / longitudinal_step - 1) * latera_num + 1
	return id


def get_children_node(longitudinal_step, lateral_step):
	"""

	:param longitudinal_step:
	:param lateral_number:
	:return:
	"""
	motion = []
	for i in range(lateral_num):
		tmp_motion = [longitudinal_step, (i - numberUnderRefLine) * lateral_step + refLineRho]
		motion.append(tmp_motion)

	return motion


def determine_goal(closed_set):
	"""

	:param closed_set:
	:return: 找到最后一列采样点里面的目标点
	"""
	tem_dic = dict()
	print('len(closedset):', len(closed_set))
	print('closedSet:', closed_set)
	for i in range(last_column_id[1], last_column_id[0] + 1):
		print(closed_set[i])
		tem_dic[i] = closed_set[i]

	c_id = min(tem_dic, key=lambda o: tem_dic[o].cost)
	goal = tem_dic[c_id]

	return goal


def calc_final_path(ngoal, closedset):
	# generate final course
	rx, ry = [ngoal.x], [ngoal.y]
	pind = ngoal.pind
	while pind != -1:
		n = closedset[pind]
		rx.append(n.x)
		ry.append(n.y)
		pind = n.pind

	return rx, ry


def plotGraph():
	print(__file__ + " start!!")
	plt.figure(figsize=(3.5, 3.0))  # 单位英寸， 3.5
	p1 = [0.15, 0.15, 0.80, 0.8]
	plt.axes(p1)
	plt.axis("equal")
	plt.grid(linestyle="--", linewidth=0.5, alpha=1)

	# # 计算多车道环境的弧长参数曲线的系数
	# efficients = multiLane.saveEfficients()

	# 计算回环环境的弧长参数曲线的系数
	efficients = cubicSpline.saveEfficients()

	# if show_lane:
	# 	# show multilane
	# 	mlane.curvePath()

	if show_obstacle:
		for i in range(len(obstacle)):
			x, y, theta = ftc.frenetToXY(obstacle[i][0], obstacle[i][1], obstacleHeading, efficients)
			car.simVehicle([x, y], theta, 'r', 0.8)

	#	sampling_point = PathLattice.sampling(longitudinal_num, lateral_num, latera_step, longitudinal_step)
	#	用来显示lattice图像
	#   PathLattice.generate_lattice()

	closed_set = dijkstra_planning(start_SRho, longi_step, latera_step, lateral_num, efficients)
	#   print('----------------------------------')
	#   print('closed_set:', len(closed_set))

	goal = determine_goal(closed_set)
	rx, ry = calc_final_path(goal, closed_set)
	#	print("rx, ry: %s" % rx, ry)

	tmp_s = []
	tmp_rho = []
	tmp_thetaRho = []
	for i in range(len(rx) - 1):
		point_s = [rx[-(i + 1)], ry[-(i + 1)], 0.0 * math.pi / 180.0]
		point_e = [rx[-(i + 2)], ry[-(i + 2)], 0.0 * math.pi / 180.0]
		s, rho, thetaRho = cubic.Polynomial(point_s, point_e)
		tmp_s.extend(s)
		tmp_rho.extend(rho)
		tmp_thetaRho.extend(thetaRho)

	if showPathLattice:
		pathlattice.ToPathPlanner2(efficients)
	x = []
	y = []
	theta = []
	if show_rough_path:
		# print(s)
		for j in range(len(tmp_s)):
			tmpX, tmpY, tmpTheta = ftc.frenetToXY(tmp_s[j], tmp_rho[j], tmp_thetaRho[j], efficients)
			x.append(tmpX)
			y.append(tmpY)
			theta.append(tmpTheta)
		# plt.scatter(end_set[0, i][0], end_set[0, i][1], color='b', s=2, alpha=0.8)
		plt.plot(x, y, 'magenta', linewidth=0.5, alpha=1)

	if show_vehicle_animation:
		x0, y0, theta0 = ftc.frenetToXY(start_SRho[0], start_SRho[1], start_SRho[2], efficients)
		car.simVehicle([x0, y0], theta0, 'g', 1)

		time_stamp = 0
		for j in range(0, len(tmp_s), 50):
			time_stamp = tmp_s[j] / s_max
			car.simVehicle([x[j], y[j]], theta[j], 'b', time_stamp)
			plt.pause(0.01)

	font1 = {'family': 'Times New Roman',
	         'weight': 'normal',
	         'size': 10,
	         }
	plt.xlabel("x (m)", font1)
	plt.ylabel('y (m)', font1)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)
	plt.xlim(-10, 30)
	plt.ylim(-5, 80)
	# plt.savefig('/home/ming/桌面/PTPSim/SimGraph/pathPlannerXY_5_30_8_48.svg')
	plt.savefig('../SimGraph/pathPlanner2_061101.tiff', dpi=600)
	plt.show()


def ToPathOptimizer2(efficients):
	closed_set = dijkstra_planning(start_SRho, longi_step, latera_step, lateral_num, efficients)
	#   print('----------------------------------')
	#   print('closed_set:', len(closed_set))

	goal = determine_goal(closed_set)
	rx, ry = calc_final_path(goal, closed_set)
	#	print("rx, ry: %s" % rx, ry)

	tmp_s = []
	tmp_rho = []
	tmp_thetaRho = []
	for i in range(len(rx) - 1):
		point_s = [rx[-(i + 1)], ry[-(i + 1)], 0.0 * math.pi / 180.0]
		point_e = [rx[-(i + 2)], ry[-(i + 2)], 0.0 * math.pi / 180.0]
		s, rho, thetaRho = cubic.Polynomial(point_s, point_e)
		tmp_s.extend(s)
		tmp_rho.extend(rho)
		tmp_thetaRho.extend(thetaRho)
	x = []
	y = []
	theta = []
	if show_rough_path:
		# print(s)
		for j in range(len(tmp_s)):
			tmpX, tmpY, tmpTheta = ftc.frenetToXY(tmp_s[j], tmp_rho[j], tmp_thetaRho[j], efficients)
			x.append(tmpX)
			y.append(tmpY)
			theta.append(tmpTheta)
		# plt.scatter(end_set[0, i][0], end_set[0, i][1], color='b', s=2, alpha=0.8)
		plt.plot(x, y, 'magenta', linewidth=0.4, alpha=1)
	return 	tmp_s, tmp_rho,	tmp_thetaRho


if __name__ == '__main__':
	plotGraph()