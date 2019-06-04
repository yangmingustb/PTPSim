import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad


# optimization parameter
max_iter = 100

cost_th = 0.01
# cost_th = 0.1
show_animation = False
speed_station_animation = 0     # 是否画出v-s图
station_time_animation = 1      # 是否画出s-t图


def generate_end_state(start, parameter):
	"""
	计算前向模拟轨迹的末端
	:param start:
	:param parameter:
	:return:
	"""
	
	s0 = start[0]
	a0 = start[1]
	a1 = start[2]
	t0 = start[3]
	
	a2 = parameter[0][0]
	a3 = parameter[1][0]
	tf = parameter[2][0]
	
	acc = a1+2*a2*(tf-t0)+3*a3*(tf-t0)**2
	speed = a0+a1*(tf-t0)+a2*(tf-t0)**2+a3*(tf-t0)**3
	
	station = a0*(tf-t0)+1/2*a1*(tf-t0)**2+1/3*a2*(tf-t0)**3+1/4*a3*(tf-t0)**4 + s0

	end_state = [station, speed, acc]
	return end_state


def calc_jacobian(start, parameter, end_state):
	"""
	
	:param start: 列表
	:param parameter: 维度3x1的数组
	:param end_state:
	:return:
	"""
	
	s0 = start[0]
	a0 = start[1]
	a1 = start[2]
	t0 = start[3]
	
	a2 = parameter[0][0]
	a3 = parameter[1][0]
	tf = parameter[2][0]
	
	vf = end_state[1]
	af = end_state[2]
	
	j00 = 1/3*(tf-t0)**3
	#   print("joo:", j00)
	j01 = 1/4*(tf-t0)**4
	#   print("j01:", j01)
	j02 = vf
	j10 = (tf-t0)**2
	j11 = (tf-t0)**3
	j12 = af
	j20 = 2*(tf-t0)
	j21 = 3*(tf-t0)**2
	j22 = 2*a2*tf+6*a3*(tf-t0)
	
	jacobian = np.array([[j00, j01, j02], [j10, j11, j12], [j20, j21, j22]])
	#   print("jacobian:", jacobian)
	
	return jacobian


def calc_diff(goal, end_state):
	"""
	计算模拟终端和真实终端的差分
	:param goal:
	:param end_state:
	:return:
	"""
	goal = np.array(goal)
	end_state = np.array(end_state)
	d = [goal[0] - end_state[0],
		goal[1] - end_state[1],
		goal[2] - end_state[2]]
	
	return d


def newton_iteration(target, start, parameters):
	"""
	牛顿迭代法，梯度下降法
	:param target:
	:param start:
	:param parameters:
	:return:
	"""
	
	for i in range(max_iter):
		end_state = generate_end_state(start, parameters)
		#   print("endstate:", end_state)
		end_state_error = np.array(calc_diff(target, end_state)).reshape(3, 1)
		
		#   求出范数,表示终止条件
		cost = np.linalg.norm(end_state_error)
		if cost <= cost_th:
			print("path is ok cost is:" + str(cost))
			break
		
		jacobian = calc_jacobian(start, parameters, end_state)
		#   print("jacobian's shape:")
		#   print(jacobian.shape)
		# numpy里面的@乘法
		try:
			dp = np.linalg.inv(jacobian) @ end_state_error
			#   print("dp's shape:")
			#   print(dp.shape)
		except np.linalg.linalg.LinAlgError:
			print("cannot calc path LinAlgError")
			parameters = None
			break
		# alpha = selection_learning_param(dp, start, parameters, target)
		alpha = 1.0
		# 将参数矩阵转换为列的形式
		parameters = np.array(parameters).reshape(3, 1)
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
	
		print("iterations:")
		print(i)
	return parameters


def main(start, end):
	
	#   start = [0, 10, 0, 0]
	#   end = [20, 30, 0]
	
	delta_station = end[0] - start[0]
	v0 = start[1]
	vf = end[1]

	tf = delta_station * 2 / (v0 + vf)+start[3]
	#   tf = 0.2
	
	init_parameter = [0, 0, tf]
	init_parameter = np.array(init_parameter).reshape(3, 1)
	
	parameters = newton_iteration(end, start, init_parameter)
	#   print("parameters:", parameters)
	
	if speed_station_animation:
		draw_speed_station_graph(start, parameters)
		
	if station_time_animation:
		draw_station_time_graph(start, parameters)
		#   plt.show()
	
	return parameters


def draw_speed_station_graph(start, parameter):
	"""
	
	:param start: 列表
	:param parameter: 维度3x1的数组
	:return:
	"""
	s0 = start[0]
	a0 = start[1]
	a1 = start[2]
	t0 = start[3]
	
	a2 = parameter[0][0]
	a3 = parameter[1][0]
	tf = parameter[2][0]
	
	time_size = 0.1
	n = int((tf-t0)/time_size)
	#   print("n:", n)
	temp_s = []
	temp_v = []
	for i in range(n):
		time_end = t0+time_size*i
		temp_parameter = np.array([a2, a3, time_end]).reshape(3, 1)
		end = generate_end_state(start, temp_parameter)
		temp_s.append(end[0])
		temp_v.append(end[1])
		
	end = generate_end_state(start, parameter)
	temp_s.append(end[0])
	temp_v.append(end[1])
	#   print("temp_s:", temp_s)
	#   print("temp_v:", temp_v)
	plt.plot(temp_s, temp_v)
	#   plt.show()


def draw_station_time_graph(start, parameter):
	
	s0 = start[0]
	a0 = start[1]
	a1 = start[2]
	t0 = start[3]
	
	a2 = parameter[0][0]
	a3 = parameter[1][0]
	tf = parameter[2][0]
	
	time_size = 0.01
	n = int((tf - t0) / time_size)
	#   print("n:", n)
	temp_s = []
	temp_t = []
	temp_a = []
	temp_v = []
	for i in range(n):
		time_end = t0 + time_size * i
		temp_parameter = np.array([a2, a3, time_end]).reshape(3, 1)
		tmp_end = generate_end_state(start, temp_parameter)
		temp_s.append(tmp_end[0])
		temp_a.append(tmp_end[2])
		temp_t.append(time_end)
		temp_v.append(tmp_end[1])
	
	tmp_end = generate_end_state(start, parameter)
	temp_s.append(tmp_end[0])
	temp_v.append(tmp_end[1])
	temp_t.append(tf)
	temp_a.append(tmp_end[2])
	#   print("temp_s:", temp_s)
	#   print("temp_v:", temp_v)
	plt.plot(temp_t, temp_s, color='r', linewidth=1)
	#   plt.plot(temp_t, temp_a, color='b', linewidth=1)
	#   plt.plot(temp_t, temp_v, color='g', linewidth=1)
	plt.plot(temp_s, temp_v, color='c', linewidth=1)
	
	
def set_state_and_run():
	"""
	设置状态并运行
	:return:
	"""
	# 设置起点和终点状态，[s,v,acc,time_0],[s,v,acc]
	#   发现一个怪异的现象，速度不能等于0
	start = [0, 1, 0, 3]
	end = [20, 15, 0]
	
	main(start, end)


if __name__ == '__main__':
	set_state_and_run()
	plt.show()

