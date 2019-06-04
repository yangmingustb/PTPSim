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
@File  : bicycleModel.py
@Author: ming.ustb@outlook.com
@Date  : 19-5-25
@GitHub: https://github.com/yangmingustb/PTPSim
"""

import math
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt


# motion parameter
wheelBase = 2.7  # wheel base [m]
ds = 0.1  # course distanse
v = 10.0 / 3.6  # velocity [m/s]


class State:

    def __init__(self, x=0.0, y=0.0, yaw=0.0, v=0.0):
        """

        :param x:
        :param y:
        :param yaw: 弧度制
        :param v:
        """
        self.x = x
        self.y = y
        self.yaw = yaw
        self.v = v


def pi_2_pi(angle):
    """

    :param angle: 弧度制
    :return: 转换成[-pi,pi]
    """
    return (angle + math.pi) % (2*math.pi) - math.pi


def update(state, v, delta, dt, L):
    """
    :param state:
    :param v:
    :param delta: 车轮转角
    :param dt:
    :param L:
    :return:
    """

    state.v = v
    state.x = state.x + state.v * math.cos(state.yaw) * dt
    state.y = state.y + state.v * math.sin(state.yaw) * dt
    state.yaw = state.yaw + state.v / L * math.tan(delta) * dt
    state.yaw = pi_2_pi(state.yaw)

    return state


def generate_trajectory(path, init_state):

    #   设置初始状态
    state = State(x=2.0, y=0.0, yaw=np.deg2rad(45))
    x, y, yaw = [state.x], [state.y], [state.yaw]
    dt = 0.1
    kp = [0, 0.1, 0.1, 0, 0]
    for ikp in kp:
        state = update(state, v, ikp, dt, wheelBase)
        plt.plot(state.x, state.y)
        x.append(state.x)
        y.append(state.y)
        yaw.append(state.yaw)

    return x, y, yaw


if __name__ == '__main__':
    generate_trajectory(100,0)
    plt.show()