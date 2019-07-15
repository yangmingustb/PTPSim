import math


def frenetToXY(s, rho, thetaRho, efficients):
	xr, yr, thetar = search_in_coefficients(s, efficients)
	x = xr + rho * math.cos(thetar + math.pi / 2.0)
	y = yr + rho * math.sin(thetar + math.pi / 2.0)
	theta = thetar + thetaRho
	return x, y, theta


def search_in_coefficients(s, coefficients):
	"""

	:param s:
	:param coefficients:
	:return: pose of reference line
	"""

	# s_id = 0
	# for i in range(len(coefficients)):
	# 	if (s > coefficients[i][0]):
	# 		s_id = i
	# 		continue
	# 	else:
	# 		break

	s_id = binary_search(coefficients, s)

	s_start, b_vector, c_vector = coefficients[s_id]

	x = b_vector[0] + b_vector[1] * s + b_vector[2] * s ** 2 + b_vector[3] * s ** 3
	d_x = b_vector[1] + 2 * b_vector[2] * s + 3 * b_vector[3] * s ** 2

	y = c_vector[0] + c_vector[1] * s + c_vector[2] * s ** 2 + c_vector[3] * s ** 3
	d_y = c_vector[1] + 2 * c_vector[2] * s + 3 * c_vector[3] * s ** 2

	theta = math.atan2(d_y,d_x)
	return x,y,theta


def binary_search(coefficients, s):

	head = 0
	tail = len(coefficients)-1
	while (head<tail):
		mid = (tail + head) // 2
		if(coefficients[mid][0] < s):
			head = mid+1
		elif(s<coefficients[mid][0]):
			tail = mid-1
		else: return mid
	if(s< coefficients[head][0]): return head-1
	else:return head


if __name__ == '__main__':
	pass
