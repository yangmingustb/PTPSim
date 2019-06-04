"""
计算dijkstra边的权重
"""
import math


def length(father_node, son_node):
	
	dis = math.sqrt((father_node.x-son_node.x)**2 + (father_node.y - son_node.y)**2)
	return dis
	

def collision_risk(father_node, son_node, obstacle):
	"""

	:param father_node:
	:param son_node:
	:param obstacle:
	:return:
	"""
	dis1 = math.sqrt((father_node.x - obstacle[0][0]) ** 2 + (father_node.y - obstacle[0][1]) ** 2)
	dis2 = math.sqrt((son_node.x - obstacle[0][0]) ** 2 + (son_node.y - obstacle[0][1]) ** 2)
	
	dis = min(dis1, dis2)
	
	dis3 = math.sqrt((father_node.x - obstacle[1][0]) ** 2 + (father_node.y - obstacle[1][1]) ** 2)
	dis4 = math.sqrt((son_node.x - obstacle[1][0]) ** 2 + (son_node.y - obstacle[1][1]) ** 2)
	
	dis5 = min(dis3, dis4)
	
	cost = 1/(dis + dis5 + 0.0001)
	
	return cost


def reference_line_cost(father_node, son_node, reference_line):
	
	dis1 = abs(father_node.y - reference_line)
	dis2 = abs(son_node.y - reference_line)
	
	cost = dis1+dis2
	return cost
	
	
def total_cost(father_node, son_node, obstacle, reference_line):
	cost1 = collision_risk(father_node, son_node, obstacle)
	cost2 = reference_line_cost(father_node, son_node, reference_line)
	cost3 = length(father_node, son_node)
	print("cost1:%s" % cost1)
	print("cost2:%s" % cost2)
	print("cost3:", cost3)
	#   alpha1 = 550
	alpha1 = 50
	alpha2 = 1
	alpha3 = 1
	cost = alpha1 * cost1 + alpha2 * cost2 + alpha3*cost3
	return cost


if __name__ == '__main__':
	pass
