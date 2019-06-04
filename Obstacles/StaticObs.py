"""
vehicle model
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad

width = 1.8
length = 4.0
R = math.sqrt((width / 2) ** 2 + (length / 2) ** 2)


def simVehicle(center_point, heading):
	"""
	使用矩形模拟车辆,前后左右四个点
	:param center_point: 矩形中心
	:param heading: 航向角
	:return:
	"""


	theta1 = math.atan(width/length)
	theta2 = -math.atan(width/length)
	theta3 = math.pi - math.atan(width/length)
	theta4 = -math.pi + math.atan(width/length)

	#   center_position
	a, b = center_point[0], center_point[1]

	# ==========================================
	#   参数方程

	x1 = a + R * np.cos(heading + theta1)
	y1 = b + R * np.sin(heading + theta1)

	x2 = a + R * np.cos(heading + theta2)
	y2 = b + R * np.sin(heading + theta2)

	x3 = a + R * np.cos(heading + theta3)
	y3 = b + R * np.sin(heading + theta3)

	x4 = a + R * np.cos(heading + theta4)
	y4 = b + R * np.sin(heading + theta4)

	x = [x3, x4, x2, x1]
	y = [y3, y4, y2, y1]


	#	subplot = plt.fill(x, y, facecolor='b', alpha=0.8)
	return x, y




if __name__ == "__main__":
	x, y = simVehicle([0, 0], math.pi/4)
	plt.fill(x, y, facecolor='b', alpha=0.8)
	plt.axis("equal")
	plt.show()
	
	
	
	

