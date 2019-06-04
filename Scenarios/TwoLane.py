"""
画出直线车道仿真环境
关于subplot,axes的使用,figure画布，subplot子图,axes子区域
https://blog.csdn.net/C_chuxin/article/details/83994457
author:ming.ustb
date:2019/4/3
"""
import matplotlib.pyplot as plt
import model.simModel as vehicle


lane_width = 3.75


def lane():

	left_line = plt.plot([0, 100], [3.75, 3.75], color='black', linewidth=1, linestyle='-')
	#   plt.plot([0, 100], [3.75 / 2, 3.75 / 2], color='r', linewidth='0.5', linestyle='--')
	center_line = plt.plot([0, 100], [0, 0], color='black', linewidth=1, linestyle='--')
	#   plt.plot([0, 100], [-3.75 / 2, -3.75 / 2], color='r', linewidth='0.5', linestyle='--')
	right_line = plt.plot([0, 100], [-3.75, -3.75], color='black', linewidth=1, linestyle='-')
	plt.axis("equal")


def dynamic_vehicle():

	for station in range(100):
		time_stamp = 0.0+station/100.0
		vehicle.simVehicle([station, 3.75/2.0], 0, 'b', time_stamp)
		plt.pause(0.01)


if __name__ == '__main__':

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	lane()
	dynamic_vehicle()
	plt.axis("equal")
	plt.grid(True)
	plt.show()
