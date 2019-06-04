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

# @File  : PathLattice.py
# @Author: ming.ustb@outlook.com
# @Date  : 19-4-29
# @GitHub: https://github.com/yangmingustb/PTPSim

"""
import numpy as np
import Curves.Cubic as cubicp
import matplotlib.pyplot as plt
import math

longi_num = 5
lateral_num = 9  # 横向采样个数
longi_step = 20.0
lateral_step = 0.5

start = [0.0, 0.0, 0.0 * math.pi / 180.0]
lane_width = 4.0


def sampling(x_row, y_column, lateral_step, longitudinal_step):
    """

    :param x_row: s采样个数
    :param y_column: lateral采样个数
    :param lateral_step: 采样步长
    :param longitudinal_step: 纵向采样步长
    :return:
    """

    end_set = np.empty(shape=[x_row, y_column, 3])

    for i in range(x_row):
        x_i = (i + 1) * longitudinal_step
        for j in range(y_column):
            y_i = (j - lane_width) * lateral_step
            target_point = [x_i, y_i, 0.0 * math.pi / 180.0]

            end_set[i, j] = np.array(target_point)

    return end_set


def generate_lattice():
    end_set = sampling(longi_num, lateral_num, lateral_step, longi_step)
    #   print(end_set)
    end_size = end_set.shape
    # print("end_size:")
    # print(end_size)
    # 生成车辆起点到第一列采样点的图
    for i in range(end_size[1]):
        s, rho, theta = cubicp.Polynomial(start, end_set[0, i])
        plt.scatter(end_set[0, i][0], end_set[0, i][1], color='b', s=2, alpha=0.8, label='sampling points')
        plt.plot(s, rho, 'r', linewidth=0.4, alpha=0.8, label='cubic polynomials')

    # 采样点之间的图
    for i in range(end_size[0] - 1):
        for j in range(end_size[1]):
            #   print([i, j])
            for q in range(end_size[1]):
                #   mptg.test_optimize_trajectory(end_set[1, 0], end_set[0, 1])
                s, rho, theta = cubicp.Polynomial(end_set[i, q], end_set[i + 1, j])
                plt.scatter(end_set[i + 1, j][0], end_set[i + 1, j][1], color='b', s=2, alpha=0.8)
                plt.plot(s, rho, 'r', linewidth=0.4, alpha=0.8)


def plot_arrow(x, y, yaw, length=2, width=0.1):
    """
    arrow函数绘制箭头,表示搜索过程中选择的航向角
    :param x:
    :param y:
    :param yaw:航向角
    :param length:
    :param width:参数值为浮点数，代表箭头尾部的宽度，默认值为0.001
    :return:
    length_includes_head：代表箭头整体长度是否包含箭头头部的长度，默认值为False
    head_width：代表箭头头部的宽度，默认值为3*width，即尾部宽度的3倍
    head_length：代表箭头头部的长度度，默认值为1.5*head_width，即头部宽度的1.5倍
    shape：参数值为'full'、'left'	、'right'，表示箭头的形状，默认值为'full'
    overhang：代表箭头头部三角形底边与箭头尾部直接的夹角关系，通过该参数可改变箭头的形状。
    默认值为0，即头部为三角形，当该值小于0时，头部为菱形，当值大于0时，头部为鱼尾状
    """
    plt.arrow(x, y, length * math.cos(yaw), length * math.sin(yaw), head_length=1.5 * length, head_width=2 * width,
              fc='lime', ec='lime')


if __name__ == '__main__':
    # plt.style.use('ggplot')
    plt.figure(figsize=(3.5, 3.5 * 0.618))  # 单位英寸， 3.5
    plt.axes([0.15, 0.2, 0.8, 0.6])
    plt.grid(linestyle="--", linewidth=0.5, alpha=1)
    # plt.rcParams['font.sans-serif'] = ['Times New Roman']  # 如果要显示中文字体，则在此处设为：SimHei
    # plt.rcParams['axes.unicode_minus'] = False  # 显示负号
    generate_lattice()
    # plt.axis("equal")

    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 10,
             }
    plt.xlabel("s (m)", font1)
    plt.ylabel(r'$\rho$ (m)', font1)

    # 设置坐标刻度值的大小以及刻度值的字体
    plt.tick_params(labelsize=10)
    # x = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0])
    # x = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
    # y = np.array([-2.0, -1, 0.0, 1.0, 2.0])
    # y = np.array([-2.0, -1, 0.0, 1.0, 2.0])

    # xgroup_labels = ['0.0', '20.0', '40.0', '60.0', '80.0', '100.0']  # x轴刻度的标识
    # ygroup_labels = ['-2.0', '-1.0', '0.0', '1.0', '2.0']  # y轴刻度的标识

    # plt.xticks(x, xgroup_labels, fontproperties='Times New Roman', fontsize=10)  # 默认字体大小为10
    # plt.yticks(y, ygroup_labels, fontproperties='Times New Roman', fontsize=10)
    plt.xticks(fontproperties='Times New Roman', fontsize=10)
    plt.yticks(fontproperties='Times New Roman', fontsize=10)
    plt.xlim(-1, 110)
    plt.ylim(-4, 4)

    '''
    plt.legend(loc=0, numpoints=1)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=10, fontweight='bold')  # 设置图例字体的大小和粗细
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)  # 去掉上边框
    ax.spines['right'].set_visible(False)  # 去掉右边框
    '''
    center_line = plt.plot([0, 105], [0, 0], color='lime', linewidth=0.6, linestyle='-', label='reference line')
    plot_arrow(105, 0, 0)

    # plt.legend(loc='upper right')
    # plt.legend()
    # 将文件保存至文件中并且画出图
    plt.savefig('/home/ming/桌面/PTPSim/SimGraph/pathLattice01.svg')
    plt.show()
