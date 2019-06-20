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
"""

import matplotlib.pyplot as plt
import TrajectoryPlanning.SpeedProfileGenerator as speed
import TrajectoryPlanning.SpeedProfileLattice as speedlattice
import TrajectoryPlanning.SpeedProfileCost as speedcost
import numpy

show_node = True
show_speed_lattice = 0
show_speed_result = 1
show_obs = False

max_t = 8  # 8s
max_s = 100  # 采样高度，100m
reso_t = 0.5
reso_s = 2
max_v = 20
inc_s = max_v * reso_t  # s方向最大增量
n_s = int(inc_s / reso_s + 1)  # s方向子节点数量
n_t = int(max_t / reso_t)
max_ns = int(max_s / reso_s) + 1  # 每一列最多的采样节点的个数

start = [0, 0]
init_v = 0  # init_speed
init_acc = 0  # init_acceleration
init_jerk = 0  # init_jerk
init_s = 0
init_t = 0

# max_acc = 2.5  # liu chang liu. 2017 IV. speed profile planning
# ref_speed = 4
static_obs = speedcost.static_obs


class Node:
    def __init__(self, x, y, cost, pind):
        """
        :param x:s里程
        :param y:t时间
        :param cost: 累计最小代价
        :param pind: 指针，指向父节点
        """
        self.x = x
        self.y = y
        self.cost = cost
        self.pind = pind

    def __str__(self):
        return str(self.x) + "," + str(self.y) + "," + str(self.cost) + "," + str(self.pind)


def dijkstra_planning():
    nstart = Node(start[0], start[1], 0.0, -1)
    open_set, closed_set = dict(), dict()  # build empty dict
    open_set[0] = nstart
    while True:
        if not open_set:
            break
        c_id = min(open_set, key=lambda o: open_set[o].cost)
        current = open_set[c_id]
        child_node = generate_son_node(current)
        #   print("c_id:", c_id)
        print("current", current, current.cost)
        #   print("child_node:", child_node)
        # 如果当前点超出范围，退出

        # show graph
        if show_node:
            plt.scatter(current.x, current.y, color='b', s=0.1)  # 把 corlor 设置为空，通过edgecolors来控制颜色

        # Remove the item from the open set
        del open_set[c_id]
        # Add it to the closed set
        closed_set[c_id] = current
        #   print(closed_set[c_id])

        # expand search grid based on motion model
        for i in range(len(child_node)):
            node = Node(child_node[i][0], child_node[i][1], current.cost, c_id)
            cost = speedcost.total_cost(node, closed_set)
            node.cost = node.cost + cost
            n_id = calc_index(node)

            if node.x > max_t:
                break
            if n_id in closed_set:
                continue
            # Otherwise if it is already in the open set
            elif n_id in open_set:
                if open_set[n_id].cost > node.cost:
                    open_set[n_id].cost = node.cost
                    open_set[n_id].pind = c_id
            else:
                open_set[n_id] = node
    #   print("closedset:", len(closed_set))
    #   print(closed_set[552])
    return closed_set


'''
def calc_index(node):
	"""
	分配索引值
	:param node:
	:return:
	"""
	tmp_n = [0]
	for i in range(n_t):
		increment_s = inc_s * (i + 1)
		if increment_s > max_s:
			increment_s = max_s

		n = increment_s / reso_s + 1
		tmp_n.append(n)

	print("----------------------")
	print("tmp_n", tmp_n)
	tmp = 1
	node_id = []
	for i in tmp_n:
		tmp = tmp + i
		node_id.append(tmp)

	node_id = node.y / reso_s + node_id[int(node.x / reso_t) - 1]
	return node_id
'''


def calc_index(node):
    """
    分配索引值
    :param node:
    :return:
    """
    total, node_list, flag = num_node()
    tmp_id = 0
    node_id = []  # 记录每一列最下面第一个节点的索引值
    for i in node_list:
        tmp_id = tmp_id + i
        node_id.append(tmp_id)

    column_id = int(node.x / reso_t) - 1
    node_id = node.y / reso_s + node_id[column_id]
    return node_id


def generate_son_node(current):
    """
    :return: 生成子节点
    """
    motion = []
    for i in range(n_s):
        son_node = [current.x + reso_t, current.y + i * reso_s]
        if son_node[0] > max_t:
            continue
        if son_node[1] > max_s:
            continue

        motion.append(son_node)
    return motion


def determine_goal(closed_set):
    """
    :param closed_set:
    :return: 找到最后边界采样点里面的目标点
    """
    #   将两个边界的点加入tmp_dic判断代价值，取最小值
    tem_dic = dict()
    total, node_list, flag = num_node()  # flag说明第一个到达max_s边界的节点的顺序，在n_t中的顺序
    first = total - node_list[-1]  # 最后一列的最下面节点的序号
    for i in range(first, total):
        tem_dic[i] = closed_set[i]

    num_flag = n_t - flag + 1  # 到达里程边界的节点的列数
    for i in range(num_flag - 1):
        tmp_id = total - 1 - (i + 1) * max_ns
        tem_dic[tmp_id] = closed_set[tmp_id]

    c_id = min(tem_dic, key=lambda o: tem_dic[o].cost)
    goal = tem_dic[c_id]
    return goal


def calc_final_path(ngoal, closedset):
    # generate final course
    rx, ry = [ngoal.x], [ngoal.y]
    pind = ngoal.pind
    while pind != -1:
        n = closedset[pind]
        rx.append(n.x)
        ry.append(n.y)
        pind = n.pind
    return rx, ry


def num_node():
    """
    :return: lattice nodes total number,包括第一个节点，即起始节点也算进去了
    """
    total = 0
    num_list = []
    tmp_flag = []
    for i in range(n_t + 1):
        tmp_num = i * (n_s - 1) + 1
        if (tmp_num - 1) * reso_s < max_s:
            total = total + tmp_num
            num_list.append(tmp_num)
        else:
            total = total + max_ns
            num_list.append(max_ns)
            tmp_flag.append(i)
    flag = tmp_flag[0]
    # print('---------')
    # print(num_list)
    return total, num_list, flag


def SpeedMain():
    print(__file__ + " start!!")
    # plt.style.use('ggplot')
    plt.figure(figsize=(3.5, 3.5 * 0.618))  # 单位英寸， 3.5
    plt.axes([0.15, 0.2, 0.8, 0.75])
    plt.grid(linestyle="--", linewidth=0.5, alpha=0.8)

    closed_set = dijkstra_planning()
    # print(len(closed_set))
    goal = determine_goal(closed_set)
    rx, ry = calc_final_path(goal, closed_set)
    # print("rx, ry: %s" % rx, ry)

    if show_speed_result:
        for i in range(len(rx) - 1):
            point_s = [rx[-(i + 1)], ry[-(i + 1)]]
            point_e = [rx[-(i + 2)], ry[-(i + 2)]]
            speed.speedprofile(point_e, point_s)
    if show_obs:
        plt.plot(static_obs[0], static_obs[1], 'c')
    if show_speed_lattice:
        speedlattice.lattice()

    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 10,
             }
    plt.xlabel("Time (s)", font1)
    plt.ylabel('s (m)', font1)


if __name__ == '__main__':
    SpeedMain()
    # plt.grid(True)
    plt.xlim(-1, 9)
    plt.ylim(-10, 110)
    plt.xticks(fontproperties='Times New Roman', fontsize=10)
    plt.yticks(fontproperties='Times New Roman', fontsize=10)
    plt.savefig('../SimGraph/speedPlanner002.tiff', dpi=600)
    plt.show()
