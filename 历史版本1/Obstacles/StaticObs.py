"""
静态障碍物列表
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad


def rectangle(center_point, heading):
	"""
	使用矩形模拟车辆,前后左右四个点
	:param center_point: 矩形中心
	:param heading: 航向角
	:return:
	"""
	width = 1.8
	length = 4.0
	
	theta1 = math.atan(width/length)
	theta2 = -math.atan(width/length)
	theta3 = math.pi - math.atan(width/length)
	theta4 = -math.pi + math.atan(width/length)
	
	R = math.sqrt((width/2)**2+(length/2)**2)
	
	#   center_position
	a, b = center_point[0], center_point[1]
	
	# ==========================================
	#   参数方程

	x1 = a + R * np.cos(heading+theta1)
	y1 = b + R * np.sin(heading+theta1)
	
	x2 = a + R * np.cos(heading + theta2)
	y2 = b + R * np.sin(heading + theta2)
	
	x3 = a + R * np.cos(heading + theta3)
	y3 = b + R * np.sin(heading + theta3)
	
	x4 = a + R * np.cos(heading + theta4)
	y4 = b + R * np.sin(heading + theta4)
	
	x = [[x1, y1], [x2, y2], [x3, y3], [x4, y4]]
	
	line1_x = [x[0][0], x[1][0]]
	line1_y = [x[0][1], x[1][1]]
	
	line2_x = [x[0][0], x[2][0]]
	line2_y = [x[0][1], x[2][1]]
	
	line3_x = [x[1][0], x[3][0]]
	line3_y = [x[1][1], x[3][1]]
	
	line4_x = [x[3][0], x[2][0]]
	line4_y = [x[3][1], x[2][1]]
	
	plt.plot(line1_x, line1_y, color='b', linewidth='1', linestyle='-')
	plt.plot(line2_x, line2_y, color='b', linewidth='1', linestyle='-')
	plt.plot(line3_x, line3_y, color='b', linewidth='1', linestyle='-')
	plt.plot(line4_x, line4_y, color='b', linewidth='1', linestyle='-')
	

if __name__ == "__main__":
	rectangle([0, 0], math.pi/4)
	plt.axis("equal")
	plt.show()
	
	
	
	

