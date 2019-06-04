"""
Dijkstra grid based planning

author: Atsushi Sakai(@Atsushi_twi)
"""

import matplotlib.pyplot as plt
import math
import TrajectoryPlanning.PathLattice as PathLattice
import TrajectoryPlanning.CostFunction as CostFunction
import Curves.CubicSpiral as spiral

show_animation = False


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


def dijkstra_planning(start, obstacle, longitudinal_step, lateral_step, longitudinal_number, lateral_number):
	"""
	
	:param start:
	:param obstacle:
	:param longitudinal_step:
	:param lateral_step:
	:param longitudinal_number:
	:param lateral_number:
	:return:
	"""
	reference_line = 0
	nstart = Node(start[0], start[1], 0.0, -1)

	child_node = get_children_node(longitudinal_step, lateral_step, lateral_number)
	
	open_set, closed_set = dict(), dict()  # build empty dict
	open_set[0] = nstart
	
	while True:
		if not open_set:
			break
			
		c_id = min(open_set, key=lambda o: open_set[o].cost)
		
		current = open_set[c_id]
		
		print("current", current)
		# 如果当前点超出范围，退出

		# show graph
		if show_animation:
			plt.plot(current.x, current.y, "xc")
		
		# Remove the item from the open set
		del open_set[c_id]
		# Add it to the closed set
		closed_set[c_id] = current
		
		# expand search grid based on motion model
		for i in range(len(child_node)):
			node = Node(current.x + child_node[i][0], child_node[i][1], current.cost, c_id)
			cost = CostFunction.total_cost(current, node, obstacle, reference_line)
			node.cost = node.cost + cost
			n_id = calc_index(node, nstart, longitudinal_step, lateral_step, lateral_number)
			
			if node.x > 60.0:
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


def calc_index(node, nstart, longitudinal_step, lateral_step, lateral_num):
	"""
	
	:param node: 子节点
	:param nstart: 起点
	:param longitudinal_step:
	:param lateral_step:
	:param lateral_num:横向采样点个数
	:return:
	"""
	# if node.x == nstart.x and node.y == nstart.y:
	#   return 0
	# else:
		
	return (node.y + 2) / lateral_step + ((node.x - nstart.x)/longitudinal_step-1)*lateral_num + 1


def get_children_node(longitudinal_step, lateral_step, lateral_number):
	"""
	
	:param longitudinal_step:
	:param lateral_step:
	:param lateral_number:
	:return:
	"""
	
	motion = [
		[longitudinal_step, (0 - 4) * lateral_step],
		[longitudinal_step, (1 - 4) * lateral_step],
		[longitudinal_step, (2 - 4) * lateral_step],
		[longitudinal_step, (3 - 4) * lateral_step],
		[longitudinal_step, (4 - 4) * lateral_step],
		[longitudinal_step, (5 - 4) * lateral_step],
		[longitudinal_step, (6 - 4) * lateral_step],
		[longitudinal_step, (7 - 4) * lateral_step],
		[longitudinal_step, (8 - 4) * lateral_step]]
	
	return motion


def determine_goal(closed_set):
	"""
	
	:param closed_set:
	:return: 找到最后一列采样点里面的目标点
	"""
	tem_dic = dict()
	for i in range(19, 28):
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


def main():
	print(__file__ + " start!!")
	
	# start and goal position
	longitudinal_number = 3
	lateral_number = 9
	lateral_step = 0.5
	longitudinal_step = 20
	
	start, sampling_point = PathLattice.sampling(longitudinal_number, lateral_number, lateral_step, longitudinal_step)
	
	PathLattice.generate_lattice()
	
	#   print(sampling_point)
	obstacle = [[20, -2], [40, 2]]

	closed_set = dijkstra_planning(start, obstacle, longitudinal_step, lateral_step, longitudinal_number, lateral_number)
	
	goal = determine_goal(closed_set)
	
	rx, ry = calc_final_path(goal, closed_set)
	print("rx, ry: %s" % rx, ry)
	
	for i in range(len(rx)-1):
		
		point_s = [rx[-(i+1)], ry[-(i+1)], 0.0 * math.pi / 180.0, 0.0]
		point_e = [rx[-(i+2)], ry[-(i+2)], 0.0 * math.pi / 180.0, 0.0]
		spiral.test_optimize_trajectory(point_e, point_s, 'k', 1)
		
	if True:
		'''
		for i in range(4):
			plt.plot(rx[i], ry[i], "ro", ms=10)
		'''
		
		plt.plot(obstacle[0][0], obstacle[0][1], "ks", ms=15)
		plt.plot(obstacle[1][0], obstacle[1][1], "ks", ms=15)
		plt.plot(start[0], start[1], "go", ms=10)
		#   plt.plot(sampling_point[:, :, 0], sampling_point[:, :, 1], "bo", ms=5)
		plt.axis("equal")
		plt.show()


if __name__ == '__main__':
	main()
