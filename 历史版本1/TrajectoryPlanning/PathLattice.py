
import numpy as np
import Curves.CubicSpiral as spiral
import matplotlib.pyplot as plt
import math


def sampling(x_row, y_column, lateral_step, longitudinal_step):
	"""
	
	:param x_row: s采样个数
	:param y_column: lateral采样个数
	:param lateral_step: 采样步长
	:param longitudinal_step: 纵向采样步长
	:return:
	"""
	#   x_row = 2  # x采样个数
	#   y_column = 10    # y采样个数
	start = [0.0, 0.0, 0.0 * math.pi / 180.0, 0.0]
	end_set = np.empty(shape=[x_row, y_column, 4])

	for i in range(x_row):
		x_i = (i+1)*longitudinal_step
		for j in range(y_column):
			y_i = (j-4)*lateral_step
			target_point = [x_i, y_i, 0.0 * math.pi / 180.0, 0.0]

			end_set[i, j] = np.array(target_point)

	return start, end_set


def generate_lattice():
	"""

	:return:
	"""

	#   start = [0.0, 0.0, 0.0 * math.pi / 180.0, 0.0]
	start, end_set = sampling(3, 9, 0.5, 20)
	#   print(end_set)
	end_size = end_set.shape
	#   print("end_size:")
	#   print(end_size)
	# 生成车辆起点到第一列采样点的图
	for i in range(end_size[1]):
		spiral.test_optimize_trajectory(end_set[0, i], start, 'r', 0.1)

	# 采样点之间的图

	for i in range(end_size[0]-1):
		for j in range(end_size[1]):
			#   print([i, j])
			for q in range(end_size[1]):
				#   mptg.test_optimize_trajectory(end_set[1, 0], end_set[0, 1])
				spiral.test_optimize_trajectory(end_set[i+1, j], end_set[i, q], 'r', 0.1)

	plt.axis("equal")
	plt.grid(True)


if __name__ == '__main__':
	generate_lattice()
	plt.show()

