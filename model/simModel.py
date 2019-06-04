"""
MIT License

Copyright (c) 2019 ming

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

===============================

# @File  : simModel.py
# @Author: ming.ustb@outlook.com
# @Date  : 19-4-29
# @GitHub: https://github.com/yangmingustb/PTPSim

"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import quad


width = 1.8
length = 4.0
R = math.sqrt((width / 2) ** 2 + (length / 2) ** 2)


def simVehicle(center_point, heading, color='b', timestamp=1.0):
    """
    使用矩形模拟车辆,前后左右四个点
    3 1
    4 2
    :param center_point: 矩形中心
    :param heading: 航向角
    :return:
    """

    theta1 = math.atan2(width,length)
    theta2 = -math.atan2(width,length)
    theta3 = math.pi + theta2
    theta4 = -math.pi + theta1

    #   center_position
    a, b = center_point[0], center_point[1]

    # ==========================================
    #   参数方程

    x1 = a + R * np.cos(heading + theta1)
    y1 = b + R * np.sin(heading + theta1)
    x2 = a + R * np.cos(heading + theta2)
    y2 = b + R * np.sin(heading + theta2)
    x3 = a + R * np.cos(heading + theta3)
    y3 = b + R * np.sin(heading + theta3)
    x4 = a + R * np.cos(heading + theta4)
    y4 = b + R * np.sin(heading + theta4)
    x = [x3, x4, x2, x1]
    y = [y3, y4, y2, y1]

    plt.fill(x, y, facecolor=color, alpha=timestamp)

    plt.plot([x1, x2], [y1, y2], color='black', linewidth=0.2, linestyle='-',alpha=timestamp)
    plt.plot([x2, x4], [y2, y4], color='black', linewidth=0.2, linestyle='-', alpha=timestamp)
    plt.plot([x1, x3], [y1, y3], color='black', linewidth=0.2, linestyle='-',alpha=timestamp)
    plt.plot([x3, x4], [y3, y4], color='black', linewidth=0.2, linestyle='-',alpha=timestamp)
    plt.axis("equal")


if __name__ == "__main__":
    simVehicle([0, 0], math.pi / 4, 'b', 0.8)
    plt.show()





