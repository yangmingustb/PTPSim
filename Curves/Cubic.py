"""
cubic polynomial
线性方程组的解法
Aa=b
在笛卡尔空间采样的，多种求解模式都不满意
改为在frenet坐标系下采样.

"""
from scipy.linalg import solve
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import FrenetMath.calHeading as heading

# state_init = [-1.0, -1.0, -180.0 * math.pi / 180.0]
# state_goal = [-10.0, -2.0, -180.0 * math.pi / 180.0]
state_init = [0.0, 0.0, -80.0 * math.pi / 180.0]
state_goal = [20.0, -40.0, -40.0 * math.pi / 180.0]


def Polynomial(x_init, x_goal):

	s0 = x_init[0]
	rho0 = x_init[1]
	theta0 = x_init[2]

	sg = x_goal[0]
	rhog = x_goal[1]
	thetag = x_goal[2]

	A = np.array([[1, s0, s0 ** 2, s0 ** 3],
	              [1, sg, sg ** 2, sg ** 3],
	              [0, 1, 2 * s0, 3 * s0 ** 2],
	              [0, 1, 2 * sg, 3 * sg ** 2]])

	b = np.array([rho0, rhog, np.tan(theta0), np.tan(thetag)])

	a = solve(A, b)
	#   print(a)

	s = np.arange(s0, sg, 0.1)
	rho = a[0] + a[1] * s + a[2] * s ** 2 + a[3] * s ** 3
	frenet_theta = []
	for i in range(len(s)-1):
		state0 = [s[i], rho[i]]
		state1 = [s[i + 1], rho[i + 1]]
		tmp_theta = heading.calHeadingFrenet(state0, state1)
		frenet_theta.append(tmp_theta)

	frenet_theta.extend([frenet_theta[-1]])

	return s, rho, frenet_theta


if __name__ == '__main__':
	start_time = time.time()
	s, rho, frenet_theta = Polynomial(state_init, state_goal)
	print('len(s):', len(s))
	end_time = time.time()
	spiral_time = end_time - start_time
	print("polynomial_time:", spiral_time)
	plt.plot(s, rho)
	plt.axis("equal")
	plt.grid(True)
	plt.show()
