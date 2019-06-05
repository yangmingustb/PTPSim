import math


def frenetToXY(s, rho, thetaRho, efficients):
	xr, yr, thetar = findEfficients(s, efficients)
	x = xr + rho * math.cos(thetar + math.pi / 2.0)
	y = yr + rho * math.sin(thetar + math.pi / 2.0)
	theta = thetar + thetaRho
	return x, y, theta


def findEfficients(s, efficients):
	"""

	:param s:
	:param efficients:
	:return: pose of reference line
	"""

	s_id = 0
	for i in range(len(efficients)):
		if (s > efficients[i][0]):
			s_id = i
			continue
		else:
			break

	s_start, b_vector, c_vector = efficients[s_id]

	x = b_vector[0] + b_vector[1] * s + b_vector[2] * s ** 2 + b_vector[3] * s ** 3
	d_x = b_vector[1] + 2 * b_vector[2] * s + 3 * b_vector[3] * s ** 2

	y = c_vector[0] + c_vector[1] * s + c_vector[2] * s ** 2 + c_vector[3] * s ** 3
	d_y = c_vector[1] + 2 * c_vector[2] * s + 3 * c_vector[3] * s ** 2

	theta = math.atan2(d_y,d_x)
	return x,y,theta


if __name__ == '__main__':
	pass
