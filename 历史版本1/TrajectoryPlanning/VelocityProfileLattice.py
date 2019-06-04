"""
generate speed profile lattice
2019/4/5
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import TrajectoryPlanning.SpeedProfileGenerator as speed_generator


def speed_sampling(x_row, y_column, speed_step, longitudinal_step):
	
	start = [0, 10, 0, 0]
	
	#   x_row = 3  # x采样个数
	#   y_column = 6    # y采样个数
	end_set = np.empty(shape=[x_row, y_column, 3])

	for i in range(x_row):
		x_i = (i+1)*longitudinal_step
		for j in range(y_column):
			y_i = j*speed_step
			target_point = [x_i, y_i, 0.0]

			end_set[i, j] = np.array(target_point)

	return start, end_set

	
def generate_lattice():
	"""

	:return:
	"""

	start, end_set = speed_sampling(3, 6, 4, 20)
	#   print(end_set)
	end_size = end_set.shape
	#   print("end_size:")
	#   print(end_size)
	
	son_node_tmp = generate_son_node(start)
	print("son_node_set", son_node_tmp)
	print(len(son_node_tmp))
	#   plt.pause(1)
	#   generate_son_node([20, 0.001, 0, 4.0])
	#   plt.pause(1)
	
	for i in range(len(son_node_tmp)):
		
		k = generate_son_node(son_node_tmp[i])
		print("k", k)
		
		for j in range(len(k)):
			m = generate_son_node(k[j])
			
	#   plt.show()
	
	plt.grid(True)


def generate_son_node(father_node):
	"""
	生成子节点信息，[s,v,acc,tf]
	:param father_node:
	:return:
	"""
	
	son_node_set = []
	for i in range(6):

		p0 = father_node[0] + 20
		p1 = i*4
		if p1 < 0.00000001:
			p1 = 0.00000001
		p2 = 0
		son_node = [p0, p1, p2]
		if father_node[1]+son_node[1] < 0.1:
			continue
		
		if father_node[1] < 0.000001:
			father_node[1] = 0.000001
		p = speed_generator.main(father_node, son_node)
		
		son_node = [p0, p1, p2, p[2][0]]
		
		son_node_set.append(son_node)
		
	return son_node_set


if __name__ == '__main__':
	generate_lattice()
	plt.show()
