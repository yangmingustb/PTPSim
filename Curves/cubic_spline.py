"""
Cubic spline planner

Author: Atsushi Sakai(@Atsushi_twi)

"""
import math
import numpy as np
import bisect
import matplotlib.pyplot as plt
import FrenetMath.referenceLine as refLine
from scipy.linalg import solve
from scipy import integrate


LaneWidth = 3.75  # [m]
left = 1
right = -1

p0 = [0.0, 0.0, 70.0 * math.pi / 180.0]
p1 = [50.0, 100.0, 40.0 * math.pi / 180.0]
p2 = [100.0, 120.0, 0.0 * math.pi / 180.0]
p3 = [150.0, 100.0, -30.0 * math.pi / 180.0]
p4 = [220.0, 150.0, 70.0 * math.pi / 180.0]
p5 = [300.0, 180.0, -50.0 * math.pi / 180.0]
p6 = [350.0, 150.0, 0.0 * math.pi / 180.0]
p7 = [400.0, 110.0, -40.0 * math.pi / 180.0]
p8 = [430.0, 20.0, -80.0 * math.pi / 180.0]
p9 = [370.0, -80.0, 20.0 * math.pi / 180.0]
p10 = [300.0, -80.0, 0.0 * math.pi / 180.0]
p11 = [200.0, -80.0, 0.0 * math.pi / 180.0]
p12 = [50.0, -100.0, -30.0 * math.pi / 180.0]
p13 = [0, -20, 0.0 * math.pi / 180.0]

p14 = [20, 70, 0.0 * math.pi / 180.0]

x = [p0[0], p14[0], p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0], p8[0],
     p9[0], p10[0], p11[0], p12[0], p13[0], p0[0]]
y = [p0[1], p14[1], p1[1], p2[1], p3[1], p4[1], p5[1], p6[1], p7[1], p8[1],
     p9[1], p10[1], p11[1], p12[1], p13[1], p0[1]]
ds = 0.1  # [m] distance of each intepolated points
s0 = 0


class Spline:
	"""
	Cubic Spline class
	"""

	def __init__(self, x, y):
		self.b, self.c, self.d, self.w = [], [], [], []

		self.x = x
		self.y = y

		self.nx = len(x)  # dimension of x
		h = np.diff(x)

		# calc coefficient c
		self.a = [iy for iy in y]

		# calc coefficient c
		A = self.__calc_A(h)
		B = self.__calc_B(h)
		self.c = np.linalg.solve(A, B)
		#  print(self.c1)

		# calc spline coefficient b and d
		for i in range(self.nx - 1):
			self.d.append((self.c[i + 1] - self.c[i]) / (3.0 * h[i]))
			tb = (self.a[i + 1] - self.a[i]) / h[i] - h[i] * \
			     (self.c[i + 1] + 2.0 * self.c[i]) / 3.0
			self.b.append(tb)

	def calc(self, t):
		"""
		Calc position

		if t is outside of the input x, return None

		"""

		if t < self.x[0]:
			return None
		elif t > self.x[-1]:
			return None

		i = self.__search_index(t)
		dx = t - self.x[i]
		result = self.a[i] + self.b[i] * dx + \
		         self.c[i] * dx ** 2.0 + self.d[i] * dx ** 3.0

		return result

	def calcd(self, t):
		"""
		Calc first derivative

		if t is outside of the input x, return None
		"""

		if t < self.x[0]:
			return None
		elif t > self.x[-1]:
			return None

		i = self.__search_index(t)
		dx = t - self.x[i]
		result = self.b[i] + 2.0 * self.c[i] * dx + 3.0 * self.d[i] * dx ** 2.0
		return result

	def calcdd(self, t):
		"""
		Calc second derivative
		"""

		if t < self.x[0]:
			return None
		elif t > self.x[-1]:
			return None

		i = self.__search_index(t)
		dx = t - self.x[i]
		result = 2.0 * self.c[i] + 6.0 * self.d[i] * dx
		return result

	def __search_index(self, x):
		"""
		search data segment index
		"""
		return bisect.bisect(self.x, x) - 1

	def __calc_A(self, h):
		"""
		calc matrix A for spline coefficient c
		"""
		A = np.zeros((self.nx, self.nx))
		A[0, 0] = 1.0
		for i in range(self.nx - 1):
			if i != (self.nx - 2):
				A[i + 1, i + 1] = 2.0 * (h[i] + h[i + 1])
			A[i + 1, i] = h[i]
			A[i, i + 1] = h[i]

		A[0, 1] = 0.0
		A[self.nx - 1, self.nx - 2] = 0.0
		A[self.nx - 1, self.nx - 1] = 1.0
		#  print(A)
		return A

	def __calc_B(self, h):
		"""
		calc matrix B for spline coefficient c
		"""
		B = np.zeros(self.nx)
		for i in range(self.nx - 2):
			B[i + 1] = 3.0 * (self.a[i + 2] - self.a[i + 1]) / \
			           h[i + 1] - 3.0 * (self.a[i + 1] - self.a[i]) / h[i]
		return B


class Spline2D:
	"""
	2D Cubic Spline class

	"""

	def __init__(self, x, y):
		self.s = self.__calc_s(x, y)
		self.sx = Spline(self.s, x)
		self.sy = Spline(self.s, y)

	def __calc_s(self, x, y):
		dx = np.diff(x)
		dy = np.diff(y)
		self.ds = [math.sqrt(idx ** 2 + idy ** 2)
		           for (idx, idy) in zip(dx, dy)]
		s = [s0]
		s.extend(np.cumsum(self.ds))
		return s

	def calc_position(self, s):
		"""
		calc position
		"""
		x = self.sx.calc(s)
		y = self.sy.calc(s)

		return x, y

	def calc_curvature(self, s):
		"""
		calc curvature
		"""
		dx = self.sx.calcd(s)
		ddx = self.sx.calcdd(s)
		dy = self.sy.calcd(s)
		ddy = self.sy.calcdd(s)
		k = (ddy * dx - ddx * dy) / ((dx ** 2 + dy ** 2) ** (3 / 2))
		return k

	def calc_yaw(self, s):
		"""
		calc yaw
		"""
		dx = self.sx.calcd(s)
		dy = self.sy.calcd(s)
		yaw = math.atan2(dy, dx)
		return yaw


def calc_spline_course(x, y, ds=0.1):
	sp = Spline2D(x, y)
	s = list(np.arange(s0, sp.s[-1], ds))

	rx, ry, ryaw, rk = [], [], [], []
	for i_s in s:
		ix, iy = sp.calc_position(i_s)
		rx.append(ix)
		ry.append(iy)
		ryaw.append(sp.calc_yaw(i_s))
		rk.append(sp.calc_curvature(i_s))

	return rx, ry, ryaw, rk, s


def offsetFromRefLine(x, y, rightLeft, theta):
	"""
	计算与参考线的偏移
	:param x:
	:param y:
	:param rightLeft:
	:param theta:
	:return:
	"""
	x2 = x + rightLeft * LaneWidth * math.cos(theta + math.pi / 2.0)
	y2 = y + rightLeft * LaneWidth * math.sin(theta+ math.pi / 2.0)

	x3 = x + rightLeft *2 * LaneWidth * math.cos(theta + math.pi / 2.0)
	y3 = y + rightLeft *2 * LaneWidth * math.sin(theta+ math.pi / 2.0)

	xref = x + rightLeft * 0.5 * LaneWidth * math.cos(theta + math.pi / 2.0)
	yref = y + rightLeft * 0.5 * LaneWidth * math.sin(theta+ math.pi / 2.0)

	return x2, y2, x3, y3, xref, yref


def main():
	print("Spline 2D test")

	sp = Spline2D(x, y)
	s = np.arange(s0, sp.s[-1], ds)

	rx, ry, ryaw, rk = [], [], [], []
	for i_s in s:
		ix, iy = sp.calc_position(i_s)
		rx.append(ix)
		ry.append(iy)
		ryaw.append(sp.calc_yaw(i_s))
		rk.append(sp.calc_curvature(i_s))
	print('len_s,len_rx:',len(s), len(rx))

	plt.figure(figsize=(3.5, 3.5 * 0.9))  # 单位英寸， 3.5
	plt.axes([0.2, 0.15, 0.75, 0.75])
	plt.axis("equal")
	plt.grid(linestyle="--", linewidth=0.2, alpha=1)

	refLine.wayPointDistribution(rx, ry, ryaw, s)

	'''
	# 如果不使用参考线生成车道线，可以直接生成
	x2_list = []
	y2_list = []
	x3_list = []
	y3_list = []
	xref_list = []
	yref_list = []

	for i in range(len(rx) - 1):
		x2, y2, x3, y3, xref, yref = offsetFromRefLine(rx[i], ry[i], left, ryaw[i])
		x2_list.append(x2)
		y2_list.append(y2)
		x3_list.append(x3)
		y3_list.append(y3)
		xref_list.append(xref)
		yref_list.append(yref)
	'''
	# plt.xlim(-100, 500)
	# plt.ylim(-100, 500)
	font1 = {'family': 'Times New Roman',
	         'weight': 'normal',
	         'size': 10,
	         }
	plt.xlabel('x (m)', font1)
	plt.ylabel('y (m)', font1)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)
	# plt.plot(x, y, "xb", label="input")
	# plt.plot(rx, ry, c="black", ls='-', linewidth= 0.5, alpha=0.8)
	# plt.plot(x2_list, y2_list, c="black", linewidth= 0.3, alpha=0.8, ls='-')
	# plt.plot(x3_list, y3_list, c="black", linewidth= 0.5, alpha=0.8, ls='-')
	# plt.plot(xref_list, yref_list, c="g", linewidth= 0.4, alpha=1, ls='--')
	# plt.legend()
	plt.savefig('/home/ming/桌面/PTPSim/SimGraph/circleRoad2.svg')

	plt.figure(figsize=(3.5, 3.5 * 0.4))  # 单位英寸， 3.5
	plt.axes([0.2, 0.35, 0.75, 0.55])
	# plt.axis("equal")
	plt.grid(linestyle="--", linewidth=0.2, alpha=1)
	# plt.xlim(-100, 500)
	# plt.ylim(-100, 500)
	font1 = {'family': 'Times New Roman',
	         'weight': 'normal',
	         'size': 10,
	         }
	plt.xlabel('s (m)', font1)
	plt.ylabel(r'$\theta$ (deg)', font1)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)
	plt.plot(s, [np.rad2deg(iyaw) for iyaw in ryaw], "-r", label="yaw", linewidth= 0.5, alpha=0.8)
	# plt.legend()
	plt.savefig('/home/ming/桌面/PTPSim/SimGraph/circleRoad2_heading.svg')

	plt.figure(figsize=(3.5, 3.5 * 0.4))  # 单位英寸， 3.5
	plt.axes([0.2, 0.35, 0.75, 0.55])
	plt.grid(linestyle="--", linewidth=0.2, alpha=1)
	# plt.xlim(-100, 500)
	# plt.ylim(-100, 500)
	font1 = {'family': 'Times New Roman',
	         'weight': 'normal',
	         'size': 10,
	         }
	plt.xlabel('s (m)', font1)
	plt.ylabel('kappa (1/m)', font1)
	plt.xticks(fontproperties='Times New Roman', fontsize=10)
	plt.yticks(fontproperties='Times New Roman', fontsize=10)
	plt.plot(s, rk, "-r", label="curvature", linewidth= 0.5, alpha=0.8)
	# plt.legend()
	plt.savefig('/home/ming/桌面/PTPSim/SimGraph/circleRoad2_kappa.svg')


def saveEfficients():
	sp = Spline2D(x, y)
	s = np.arange(s0, sp.s[-1], ds)

	rx, ry, ryaw, rk = [], [], [], []
	for i_s in s:
		ix, iy = sp.calc_position(i_s)
		rx.append(ix)
		ry.append(iy)
		ryaw.append(sp.calc_yaw(i_s))
		rk.append(sp.calc_curvature(i_s))
	# print('len_s,len_rx:', len(s), len(rx))

	efficients = refLine.wayPointDistribution(rx, ry, ryaw, s)

	return efficients


if __name__ == '__main__':

	main()

	plt.show()
