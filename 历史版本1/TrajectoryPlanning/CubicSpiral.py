"""
Name:Model trajectory generator
--------------------
reference paper lists:
nagy and kelly 2001
kelly 2003
howard 07 IJRR
matthew 2011 ICRA and Phd thesis 2011
xu wenda ICRA 2012
li xiaohui 2017
-------------------
@author:yang ming ustb
-----------------------
说明：这个问题是一个病态问题，控制量参数化的参数给以微小扰动，
导致计算结果的剧烈变化，对初值也十分敏感。
同时，对差分时的微小扰动也很敏感。
----------------------
以后改动方向：使用最优控制的方法对多项式进行求解，放弃螺旋多项式
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad

# optimization parameter
max_iter = 100
h = np.array([0.001, 0.0001, 0.0001, 0.1]).T  # 差分使用的一个小扰动
cost_th = 0.01
# cost_th = 0.1
show_animation = False


def generate_last_state(start, parameter):
	"""
	计算前向模拟轨迹的末端
	:param start:
	:param parameter:
	:return:
	"""
	
	a = start[3]
	theta0 = start[2]
	b = parameter[0]
	c = parameter[1]
	d = parameter[2]
	sf = parameter[3]
	k = a + b * sf + c * sf ** 2 + d * sf ** 3
	
	theta = theta0 + a * sf + b * sf ** 2 / 2 + c * sf ** 3 / 3 + d * sf ** 4 / 4
	
	# theta = theta * 180.0 / math.pi     # 角度制
	dk = b + 2 * c * sf + 3 * d * sf ** 2
	
	x = quad(lambda s: np.cos(theta0 + a * s + b * s ** 2 / 2.0 + c * s ** 3 / 3.0 + d * s ** 4 / 4.0), 0, sf)
	#   返回的是一个元组（value ,error),取value
	x = x[0]
	y = quad(lambda s: np.sin(theta0 + a * s + b * s ** 2 / 2.0 + c * s ** 3 / 3.0 + d * s ** 4 / 4.0), 0, sf)
	y = y[0]
	
	x = x + start[0]
	y = y + start[1]
	end_state = [x, y, theta, k]
	return end_state


def generate_trajectory(start, parameter):
	"""
	计算前向模拟轨迹
	:param start:
	:param parameter:
	:return:
	"""
	
	a = start[3]
	theta0 = start[2]
	b = parameter[0]
	c = parameter[1]
	d = parameter[2]
	sf = parameter[3]
	k = a + b * sf + c * sf ** 2 + d * sf ** 3
	
	theta = theta0 + a * sf + b * sf ** 2 / 2 + c * sf ** 3 / 3 + d * sf ** 4 / 4
	# theta = theta * 180.0 / math.pi
	dk = b + 2 * c * sf + 3 * d * sf ** 2
	s0 = 0
	step = 0.1
	n = (sf - s0) / step
	n = int(n)
	end_state = []
	for s_n in range(n + 1):
		s_i = s0 + s_n * step
		x = quad(lambda s: np.cos(theta0 + a * s + b * s ** 2 / 2 + c * s ** 3 / 3 + d * s ** 4 / 4), 0, s_i)
		y = quad(lambda s: np.sin(theta0 + a * s + b * s ** 2 / 2 + c * s ** 3 / 3 + d * s ** 4 / 4), 0, s_i)
		#   返回的是一个元组（value ,error),取value
		x = x[0]
		y = y[0]
		
		x = x + start[0]
		y = y + start[1]
		
		end_state.append([x, y])
	
	return end_state


def pi_2_pi(angle):
	"""

	:param angle: 弧度制
	:return: 转换成[-pi,pi]
	"""
	return (angle + math.pi) % (2 * math.pi) - math.pi


def calc_diff(goal, end_state):
	"""
	计算差分
	:param goal:真实目标
	:param end_state:模拟终端
	:return:
	"""
	goal = np.array(goal)
	end_state = np.array(end_state)
	d = np.array([
		goal[0] - end_state[0],
		goal[1] - end_state[1],
		pi_2_pi(goal[2] - end_state[2]),
		goal[3] - end_state[3]])
	
	return d


def calc_jacobian(start, goal, parameters):
	"""
	中心差分计算雅可比矩阵
	:param start:
	:param goal:
	:param parameters:
	:return:
	"""
	goal = np.array(goal)
	start = np.array(start)
	parameters = np.array(parameters)
	
	parameters_b1 = [parameters[0] + h[0], parameters[1], parameters[2], parameters[3]]
	end_statep = generate_last_state(start, parameters_b1)
	dp = calc_diff(goal, end_statep)
	
	parameters_b2 = [parameters[0] - h[0], parameters[1], parameters[2], parameters[3]]
	end_staten = generate_last_state(start, parameters_b2)
	dn = calc_diff(goal, end_staten)
	
	d1 = np.array((dp - dn) / (2.0 * h[0])).reshape(4, 1)
	
	# --------------------------------------
	parameters_c1 = [parameters[0], parameters[1] + h[1], parameters[2], parameters[3]]
	end_statep = generate_last_state(start, parameters_c1)
	dp = calc_diff(goal, end_statep)
	
	parameters_c2 = [parameters[0], parameters[1] - h[1], parameters[2], parameters[3]]
	end_staten = generate_last_state(start, parameters_c2)
	dn = calc_diff(goal, end_staten)
	
	d2 = np.array((dp - dn) / (2.0 * h[1])).reshape(4, 1)
	
	#   -------------------------------------------
	parameters_d1 = [parameters[0], parameters[1], parameters[2] + h[2], parameters[3]]
	end_statep = generate_last_state(start, parameters_d1)
	dp = calc_diff(goal, end_statep)
	
	parameters_d2 = [parameters[0], parameters[1], parameters[2] - h[2], parameters[3]]
	end_staten = generate_last_state(start, parameters_d2)
	dn = calc_diff(goal, end_staten)
	
	d3 = np.array((dp - dn) / (2.0 * h[2])).reshape(4, 1)
	
	#   ---------------------------------------
	parameters_s1 = [parameters[0], parameters[1], parameters[2], parameters[3] + h[3]]
	end_statep = generate_last_state(start, parameters_s1)
	dp = calc_diff(goal, end_statep)
	
	parameters_s2 = [parameters[0], parameters[1], parameters[2], parameters[3] - h[3]]
	end_staten = generate_last_state(start, parameters_s2)
	dn = calc_diff(goal, end_staten)
	
	d4 = np.array((dp - dn) / (2.0 * h[3])).reshape(4, 1)
	
	jacobian = np.hstack((d1, d2, d3, d4))
	
	return jacobian


def selection_learning_param(dp, start, parameters, goal):
	"""
	为了数值稳定性，一般加一个系数
	:param dp:
	:param parameters:
	:param start:
	:param goal:
	:return:
	"""
	goal = np.array(goal)
	start = np.array(start)
	parameters = np.array(parameters)
	
	mincost = float("inf")
	mina = 1.0
	maxa = 2.0
	da = 0.5
	
	for a in np.arange(mina, maxa, da):
		p = parameters + a * dp
		end_state = generate_last_state(start, p)
		end_state_error = calc_diff(goal, end_state)
		cost = np.linalg.norm(end_state_error)
		#   print(cost)
		#   print("cost")
		#   print(a)
		#   print("a")
		if cost <= mincost and a != 0:
			mina = a
			mincost = cost
	
	#  print(mincost, mina)
	#  input()
	
	return mina


def show_trajectory(target, xc, yc):
	plt.clf()
	# plot_arrow(target.x, target.y, target.yaw)
	plt.plot(xc, yc, "-r")
	plt.axis("equal")
	plt.grid(True)
	# plt.pause(0.1)


def optimize_trajectory(target, start, parameters):
	"""
	牛顿迭代法，梯度下降法
	:param target:
	:param start:
	:param parameters:
	:return:
	"""
	
	for i in range(max_iter):
		end_state = generate_last_state(start, parameters)
		end_state_error = np.array(calc_diff(target, end_state)).reshape(4, 1)
		
		#   求出范数,表示终止条件
		cost = np.linalg.norm(end_state_error)
		if cost <= cost_th:
			print("path is ok cost is:" + str(cost))
			break
		
		jacobian = calc_jacobian(start, target, parameters)
		#   print("jacobian's shape:")
		#   print(jacobian.shape)
		# numpy里面的@乘法
		try:
			dp = - np.linalg.inv(jacobian) @ end_state_error
			#   print("dp's shape:")
			#   print(dp.shape)
		except np.linalg.linalg.LinAlgError:
			print("cannot calc path LinAlgError")
			parameters = None
			break
		# alpha = selection_learning_param(dp, start, parameters, target)
		alpha = 1.0
		# 将参数矩阵转换为列的形式
		parameters = np.array(parameters).reshape(4, 1)
		#   print("parameters:")
		#   print("parameters's shape:")
		parameters += alpha * np.array(dp)
		#   print("parameters:")
		#   print("parameters's shape:")
		#   print(parameters.shape)
		#   print(parameters.T)
		
		if i == max_iter and cost >= cost_th:
			parameters = None
			print("cannot calc path")
		
		#   print("iterations:")
		#   print(i)
	return parameters


def test_optimize_trajectory(target, start, line_color, line_width):
	#  设置目标状态
	# target = [10.0, 4.0, 0.0, 0.0]
	#   设置起点状态
	# start = [0.0, 0.0, 0.0, 0.0]
	
	#   初始参数猜测
	p0 = 0.0
	p1 = 0.0
	p2 = 0.0
	p3 = target[0] - start[0]
	parameters_init = [p0, p1, p2, p3]
	
	parameters = optimize_trajectory(target, start, parameters_init)
	# print("参数：")
	# print(parameters)
	end_state = generate_trajectory(start, parameters)
	#   print("终点状态：")
	#   print(end_state)
	
	x = [list1[0] for list1 in end_state]
	y = [list2[1] for list2 in end_state]
	# show_trajectory(target, x, y)
	plt.plot(x, y, color=line_color, linewidth=line_width)
	# plot_arrow(target.x, target.y, target.yaw)
	
	plt.grid(True)
	#   plt.axis("equal")
	#   plt.show()


def generate_trajectories():
	"""
	为了数值稳定性，防止出现“大数吃小数”问题，
	start = [x, y, theta, curvature],角度制，因为弧度太小，引起不稳定；曲率变为
	:return:
	"""
	# start = [0.0, 0.0, 20.0*math.pi / 180.0, 0.0]
	start = [0.0, 0.0, 0.0 * math.pi / 180.0, 0.0]
	for i in range(9):
		target = [20.0, ((i - 4) / 2.0), 50.0 * math.pi / 180.0, 0.0]
		test_optimize_trajectory(target, start, 'r', 1)
	'''
	for i in range(9):
		target = [20.0, ((i - 4) / 2.0), 20.0 * math.pi / 180.0, 0.0]
		test_optimize_trajectory(target, start)

	for i in range(9):
		target = [20.0, ((i - 4) / 2.0), -20.0 * math.pi / 180.0, 0.0]
		test_optimize_trajectory(target, start)
	'''
	
	
def main():
	print(__file__ + " start!!")
	# test_trajectory_generate()
	
	# test_optimize_trajectory()
	generate_trajectories()
	plt.axis("equal")
	plt.grid(True)
	plt.show()


if __name__ == '__main__':
	main()
