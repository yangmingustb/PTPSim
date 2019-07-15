"""
计算dijkstra边的权重,路径权重
"""

import math
import FrenetMath.kappaCalculate as calkappa
import Curves.Cubic as cubic
import FrenetMath.FrenetToCartesian as fretoc


#   refLineRho = 0.0	# 参考线为道路中心线
r_circle = 1.0
d_circle = 2.0
obstacle_inflation = 1.5

alpha1 = 100
alpha2 = 1
alpha3 = 10
alpha4 = 0.0


def kappa(father_node, son_node, efficients):
	mean_kappa = calkappa.trajectory_kappa(father_node, son_node, efficients)

	return mean_kappa


def reference_line_cost(father_node, son_node, refl):
	dis1 = abs(father_node.y - refl)
	dis2 = abs(son_node.y - refl)

	cost = (dis1 + dis2)/2.0
	return cost


def collision_risk(father_node, son_node, obstacle):
	"""
	:param father_node:
	:param son_node:
	:param obstacle:
	:return:
	"""
	begin = [father_node.x, father_node.y, 0.0 * math.pi / 180.0]
	end = [son_node.x, son_node.y, 0.0 * math.pi / 180.0]
	s, rho, theta = cubic.Polynomial(begin, end)

	dis = math.inf
	for i in range(0, len(s), 5):
		for j in range(len(obstacle)):

			tmp_dis = math.sqrt((s[i] - obstacle[j][0]) ** 2 + (rho[i] - obstacle[j][1]) ** 2)
			if tmp_dis < r_circle + obstacle_inflation:
				if dis > tmp_dis:
					dis = tmp_dis

				break

	cost = 1/(dis + 0.0001)
	
	return cost


def total_cost(father_node, son_node, obstacle, refl, efficients):
	cost1 = kappa(father_node, son_node, efficients)
	cost2 = reference_line_cost(father_node, son_node, refl)
	cost3 = collision_risk(father_node, son_node, obstacle)


	# print("cost1:%s" % cost1)
	# print("cost2:%s" % cost2)
	# print("cost3:", cost3)
	#   alpha1 = 550

	cost = alpha1 * cost1 + alpha2 * cost2 + alpha3*cost3
	return cost


if __name__ == '__main__':
	pass
