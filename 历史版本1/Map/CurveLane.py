"""
画出弯道仿真环境
author:ming.ustb
date:2019/4/3
"""


import numpy as np
import matplotlib.pyplot as plt
import Obstacles.StaticObs as static
import math


# ==========================================
# 画出车道线图
# 圆的基本信息
r = 100.0
# 圆心坐标
a, b = (0., 0.)

# ==========================================
#   参数方程
theta = np.arange(np.pi/2, np.pi, 0.01)
x = a + r * np.cos(theta)
y = b + r * np.sin(theta)

x2_center = a + (r-3.75/2) * np.cos(theta)
y2_center = b + (r-3.75/2) * np.sin(theta)

x2 = a + (r-3.75) * np.cos(theta)
y2 = b + (r-3.75) * np.sin(theta)

x3_center = a + (r+3.75/2) * np.cos(theta)
y3_center = b + (r+3.75/2) * np.sin(theta)
x3 = a + (r+3.75) * np.cos(theta)
y3 = b + (r+3.75) * np.sin(theta)

fig = plt.figure()
axes = fig.add_subplot(111)


#   =====================================
# 画出移动障碍物图
heading = np.arange(-np.pi, -np.pi*3/2, -0.01)
center_x = a + (r+3.75/2) * np.cos(heading)
center_y = b + (r+3.75/2) * np.sin(heading)
real_heading = heading + math.pi*3/2

axes.plot(x, y, color='black', linewidth='2', linestyle='--')
#   axes.plot(x2_center, y2_center, color='r', linewidth='0.5', linestyle='--')
axes.plot(x2, y2, color='black', linewidth='1.5', linestyle='-')
#   axes.plot(x3_center, y3_center, color='r', linewidth='0.5', linestyle='--')
axes.plot(x3, y3, color='black', linewidth='2', linestyle='-')
axes.axis('equal')

for i in range(len(real_heading)):

	static.simVehicle([center_x[i], center_y[i]], real_heading[i])

	plt.pause(0.01)

plt.axis("equal")
plt.grid(True)
plt.show()
