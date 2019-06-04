import math
import numpy as np
import scipy.interpolate

# motion parameter
L = 1.0  # wheel base
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


def generate_trajectory(s, km, kf, k0):

    n = s / ds
    time = s / v  # [s]
    tk = np.array([0.0, time / 2.0, time])
    kk = np.array([k0, km, kf])
    t = np.arange(0.0, time, time / n)
    fkp = scipy.interpolate.interp1d(tk, kk, kind="quadratic")
    kp = [fkp(ti) for ti in t]
    dt = float(time / n)

    #  plt.plot(t, kp)
    #  plt.show()

    #   设置初始状态
    state = State(x=2.0, y=0.0, yaw=np.deg2rad(45))
    x, y, yaw = [state.x], [state.y], [state.yaw]

    for ikp in kp:
        state = update(state, v, ikp, dt, L)
        x.append(state.x)
        y.append(state.y)
        yaw.append(state.yaw)

    return x, y, yaw


def generate_last_state(s, km, kf, k0):

    n = s / ds
    time = s / v  # [s]
    dt = time / n
    tk = np.array([0.0, time / 2.0, time])
    kk = np.array([k0, km, kf])
    t = np.arange(0.0, time, dt)
    #  插值
    fkp = scipy.interpolate.interp1d(tk, kk, kind="quadratic")
    #  得到差值后的转角
    kp = [fkp(ti) for ti in t]
    #  dt = time / n

    #  plt.plot(t, kp)
    #  plt.show()

    #   设置初始状态
    state = State(x=2.0, y=0.0, yaw=np.deg2rad(45))

    [update(state, v, ikp, dt, L) for ikp in kp]
    return state.x, state.y, state.yaw
