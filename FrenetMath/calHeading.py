import math


def calHeadingXY(state0, state1):
	x0 = state0[0]
	y0 = state0[1]
	x1 = state1[0]
	y1 = state1[1]

	theta = math.atan((y1 - y0) / (x1 - x0))
	return theta


def calHeadingFrenet(state0, state1):
	s0 = state0[0]
	rho0 = state0[1]
	s1 = state1[0]
	rho1 = state1[1]

	frenet_theta = math.atan((rho1-rho0)/(s1-s0))

	return frenet_theta