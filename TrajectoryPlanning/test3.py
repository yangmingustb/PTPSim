from gekko import GEKKO
import numpy as np
import math
import FrenetMath.referenceLine as refline
import FrenetMath.FrenetToCartesian as ftc
import matplotlib.pyplot as plt


show_animation = 1
rho_max = 3.0  # 行驶通道的左右边界
rho_min = -3.0
s_max = 100.0
reso_s = 2.0
n_s = int(s_max / reso_s)
r_circle = 1.0  # 车体包络圆
d_circle = 2.0  # 圆心间距
obs_inflation = 1.0
safe_distance = obs_inflation + r_circle
kappa_max = 0.187
w_d = 1.0
w_dd = 10.0
w_ddd = 50.0
w_ref = 1.5
rho_r = [0.0 for i in range(n_s + 3)]
# 初值
rho_init_guess = np.zeros((1, n_s + 3))

# 起点状态
x0 = 0
y0 = 0
s0 = 0
rho0 = 0
theta0 = 0
kappa0 = 0
x_init = [x0, y0, s0, rho0, theta0, kappa0]

# 里程s的序列
s = [i for i in np.arange(reso_s, (s_max + 4 * reso_s), reso_s)]

# 障碍物表示
static_obs = np.array([[15, 40], [1, -1]])  # 障碍物的frenet坐标,
obs_shape = static_obs.shape
num_static = obs_shape[1]  # 静态障碍物个数

# 等式约束
xr0, yr0, thetar0 = refline.RefLine(s0)
rho_1 = reso_s * (math.tan(theta0 - thetar0)) + rho0

x1, y1 = ftc.frenetToXY(s0 + reso_s, rho_1)
s2 = s0 + 2 * reso_s
xr2, yr2, thetar2 = refline.RefLine(s2)
m = kappa0 * ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** (3 / 2)
n = (x1 - x0) * (yr2 - 2 * y1 + y0) + (y1 - y0) * (xr2 - 2 * x1 + x0)
q = (x1 - x0) * math.sin(thetar2 + math.pi / 2.0) - (y1 - y0) * math.cos(thetar2 + math.pi / 2.0)
rho_2 = (m - n) / q
# print(rho_1)
# print(rho_2)

# Initialize Model
m = GEKKO()

# help(m)

# define parameter
eq1 = m.Param(value=rho_1)
eq2 = m.Param(value=rho_2)

# initialize variables
rho_vector = [m.Var() for i in range(n_s + 3)]

# initial values
for i in rho_vector:
    i.value = 0

# lower and upper bounds
for i in rho_vector:
    i.lower = rho_min
    i.upper = rho_max

# print('rho_vector:', rho_vector)


def inequality_cons(rho):

    # 曲率约束
    for i in range(n_s):
        x0, y0 = ftc.frenetToXY(s[i], rho[i])
        x1, y1 = ftc.frenetToXY(s[i + 1], rho[i + 1])
        x2, y2 = ftc.frenetToXY(s[i + 2], rho[i + 2])
        k1 = (x1 - x0) * (y2 - 2 * y1 + y0)
        k2 = (y1 - y0) * (x2 - 2 * x1 + x0)
        k3 = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** (3.0 / 2.0)
        kappa1 = -k1 + k2 + k3 * kappa_max
        kappa2 = k1 - k2 + k3 * kappa_max
        m.Equation(kappa1 >= 0)
        m.Equation(kappa2 >= 0)

    # 障碍物约束
    for i in range(n_s + 3):
        if ((s[i] - static_obs[0][0])**2 < 9):

            c1_obs = (s[i] - static_obs[0, 0]) ** 2 + (rho[i] - static_obs[1, 0]) ** 2
            c2_obs = (s[i] + d_circle - static_obs[0, 0]) ** 2 + (rho[i] - static_obs[1, 0]) ** 2
            c3_obs = (s[i] + 2 * d_circle - static_obs[0, 0]) ** 2 + (rho[i] - static_obs[1, 0]) ** 2
            m.Equation(c1_obs >= safe_distance**2)
            m.Equation(c2_obs >= safe_distance**2)
            m.Equation(c3_obs >= safe_distance**2)
        else:
            continue

    # 障碍物2约束
    for i in range(n_s + 3):
        if ((s[i] - static_obs[1][0]) ** 2 < 9):
            c1_obs2 = (s[i] - static_obs[1, 0]) ** 2 + (rho[i] - static_obs[1, 1]) ** 2
            c2_obs2 = (s[i] + d_circle - static_obs[1, 0]) ** 2 + (rho[i] - static_obs[1, 1]) ** 2
            c3_obs2 = (s[i] + 2 * d_circle - static_obs[1, 0]) ** 2 + (rho[i] - static_obs[1, 1]) ** 2

            m.Equation(c1_obs2 >= safe_distance ** 2)
            m.Equation(c2_obs2 >= safe_distance ** 2)
            m.Equation(c3_obs2 >= safe_distance ** 2)
        else:continue


# inEquations
inequality_cons(rho_vector)
m.Equation(rho_vector[0] == eq1)
m.Equation(rho_vector[1] == eq2)


# Objective
def objective(rho):
    f_d = 0
    f_dd = 0
    f_ddd = 0
    f_ref = 0
    for i in range(n_s):
        f_d = f_d + (rho[i + 1] - rho[i]) ** 2 / (reso_s ** 2)
        f_dd = f_dd + (rho[i + 2] - 2 * rho[i + 1] + rho[i]) ** 2 / (reso_s ** 4)
        f_ddd = f_ddd + (rho[i + 3] - 3 * rho[i + 2] + 3 * rho[i + 1] - rho[i]) ** 2 / reso_s ** 6
        f_ref = f_ref + (rho[i] - rho_r[i]) ** 2

    # print('=========================================')
    # print(f_d)
    # print(f_dd)
    # print(f_ddd)
    # print(f_ref)

    f_d = reso_s * w_d * f_d
    f_dd = reso_s * w_dd * f_dd
    f_ddd = reso_s * w_ddd * f_ddd
    f_ref = reso_s * w_ref * f_ref

    f = f_d + f_dd + f_ddd + f_ref

    return f

f = objective(rho_vector)
m.Obj(f)

# Set global options
m.options.IMODE = 3  # steady state optimization

# Solve simulation
m.solve()

# Results
# print('')
# print('Results')
# print('x1: ' + str(x1.value))
# print('x2: ' + str(x2.value))
# print('x3: ' + str(x3.value))
# print('x4: ' + str(x4.value))

plt.figure(figsize=(3.5, 3.5 * 0.618))  # 单位英寸， 3.5
plt.axes([0.15, 0.2, 0.8, 0.6])
plt.grid(linestyle="--", linewidth=0.5, alpha=1)

rho_list = []
for i in rho_vector:
    rho_list.append(i.value)
plt.plot(s, rho_list)

import TrajectoryPlanning.VehicleModel as car

car.simVehicle([static_obs[0][0], static_obs[1][0]], 0 * math.pi / 180, 'r', 1)
car.simVehicle([static_obs[0][1], static_obs[1][1]], 0 * math.pi / 180, 'r', 1)

plt.xlim(-1, 110)
plt.ylim(-4, 4)
plt.show()

