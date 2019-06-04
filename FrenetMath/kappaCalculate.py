"""
calculate one trajectory kappa
input:two point
output:mean kappa
"""

import math
import Curves.Cubic as cubic
import FrenetMath.FrenetToCartesian as fretoc
import Scenarios.multiLane as mlane


def trajectory_kappa(father_node, son_node, efficients):

	begin = [father_node.x, father_node.y, 0.0 * math.pi / 180.0]
	end = [son_node.x, son_node.y, 0.0 * math.pi / 180.0]
	s, rho, thetaRho = cubic.Polynomial(begin, end)

	tmp_kappa = []
	for i in range(0, len(s), 5):
		x0, y0, theta0 = mlane.frenetToXY(s[i], rho[i],thetaRho[i], efficients)
		x1, y1, theta1 = mlane.frenetToXY(s[i + 1], rho[i + 1], thetaRho[i], efficients)
		x2, y2, theta2 = mlane.frenetToXY(s[i + 2], rho[i + 2], thetaRho[i], efficients)
		k1 = (x1 - x0)*(y2-2*y1+y0)
		k2 = (y1-y0)*(x2-2*x1+x0)
		k3 = ((x1-x0)**2+(y1-y0)**2)**(3.0/2.0)

		if (k3 == 0.0):
			tmp_kappa.append(0)

		else:
			tmp_kappa.append((k1 - k2) / k3)

	sum_kappa = 0
	for i in range(len(tmp_kappa)):
		sum_kappa = sum_kappa + tmp_kappa[i]**2

	mean_kappa = sum_kappa/len(tmp_kappa)

	return mean_kappa



if __name__ == '__main__':
	pass
