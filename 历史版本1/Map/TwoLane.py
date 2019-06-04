"""
画出直线车道仿真环境
author:ming.ustb
date:2019/4/3
"""
import matplotlib.pyplot as plt
import Obstacles.StaticObs as static


for station in range(100):
	plt.cla()
	plt.plot([0, 100], [3.75, 3.75], color='black', linewidth='2', linestyle='-')
	#   plt.plot([0, 100], [3.75 / 2, 3.75 / 2], color='r', linewidth='0.5', linestyle='--')
	plt.plot([0, 100], [0, 0], color='black', linewidth='1.5', linestyle='--')
	#   plt.plot([0, 100], [-3.75 / 2, -3.75 / 2], color='r', linewidth='0.5', linestyle='--')
	plt.plot([0, 100], [-3.75, -3.75], color='black', linewidth='2', linestyle='-')
	plt.axis("equal")
	
	static.rectangle([station, 3.75/2], 0)
	
	plt.pause(0.01)

plt.axis("equal")
plt.grid(True)
plt.show()


if __name__ == '__main__':
	pass
