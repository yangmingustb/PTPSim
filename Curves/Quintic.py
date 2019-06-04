"""
quintic polynomial
线性方程组的解法
Aa=b

特殊角度,[0,0,pi/2,0.0]
特殊顺序[0,0,-pi,0.0],end = [-10,-2, -pi,0]
一些特殊情况，曲线无法连接

"""
from scipy.linalg import solve
import numpy as np
import matplotlib.pyplot as plt
import math
import time

x_init = [-1.0, -1.0,  90.0* math.pi / 180.0, 0.0]
x_goal = [10.0, 10.0, 90.0 * math.pi / 180.0, 0.0]



def rotation_pose(pose, beta):
    """
    绕原点逆时针旋转beta
    :param xg:
    :param yg:
    :param thetag:
    :return:
    """
    x = pose[0]
    y = pose[1]
    theta = pose[2]
    x_new = x*math.cos(beta)-y*math.sin(beta)
    y_new = y*math.cos(beta)+ x*math.sin(beta)
    theta_new = theta + beta

    return [x_new, y_new, theta_new, pose[3]]


def rotation_position(position, beta):
    """
    绕原点恢复到原来坐标
    :param position:
    :param beta:
    :return:
    """

    x = position[0]
    y = position[1]
    print('y', y)
    beta = -beta
    x_new = x * math.cos(beta) - y * math.sin(beta)
    y_new = y * math.cos(beta) + x * math.sin(beta)

    return [x_new, y_new]


def QuinticPolynomial(x_init, x_goal):

    if (x_init[2]-90* math.pi / 180.0)<0.001 or (x_init[2]+90* math.pi / 180.0) < 0.001 \
        or (x_goal[2]-90* math.pi / 180.0)<0.001 or (x_goal[2]+90* math.pi / 180.0) < 0.001 :

        x_init = rotation_pose(x_init, -45* math.pi / 180.0)
        x_goal = rotation_pose(x_goal, -45 * math.pi / 180.0)

    x0 = x_init[0]
    y0 = x_init[1]
    theta0 = x_init[2]
    kappa0 = x_init[3]

    xg = x_goal[0]
    yg = x_goal[1]
    thetag = x_goal[2]
    kappag = x_goal[3]

    A = np.array([[1, x0, x0**2, x0**3, x0**4, x0**5],
                  [1, xg, xg**2, xg**3, xg**4, xg**5],
                  [0, 1, 2*x0, 3*x0**2, 4*x0**3, 5*x0**4],
                  [0, 1, 2*xg, 3*xg**2, 4*xg**3, 5*xg**4],
                  [0, 0, 2, 6*x0, 12*x0**2, 20*x0**3],
                  [0, 0, 2, 6 * xg, 12 * xg ** 2, 20 * xg ** 3]])

    b = np.array([y0, yg, np.tan(theta0), np.tan(thetag), kappa0 * (1 + np.tan(theta0)) ** (3.0 / 2.0), kappag * (1 + np.tan(thetag)) ** (3.0 / 2.0)])

    a = solve(A, b)
    #   print(a)

    x = np.arange(x0, xg, 0.1)
    y = a[0]+a[1]*x+a[2]*x**2+a[3]*x**3+a[4]*x**4+a[5]*x**5

    tmp_x = x
    tmp_y = y
    if (x_init[2]-90* math.pi / 180.0)<0.001 or (x_init[2]+90* math.pi / 180.0) < 0.001 \
        or (x_goal[2]-90* math.pi / 180.0)<0.001 or (x_goal[2]+90* math.pi / 180.0) < 0.001 :

            tmp_x = []
            tmp_y = []
            for i in range(x.size):
                position = [x[i], y[i]]
                r_position = rotation_position(position, -45* math.pi / 180.0)
                tmp_x.append(r_position[0])
                tmp_y.append(r_position[1])

    return tmp_x, tmp_y


if __name__ == '__main__':
    start_time = time.time()
    tmp_x, tmp_y = QuinticPolynomial(x_init, x_goal)
    end_time = time.time()
    spiral_time = end_time - start_time
    print("polynomial_time:", spiral_time)
    plt.plot(tmp_x, tmp_y)
    plt.axis("equal")
    plt.grid(True)
    plt.show()

