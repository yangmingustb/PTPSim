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
# @File  :
# @Author: ming.ustb@outlook.com
# @Date  : 19-4-29
# @GitHub: https://github.com/yangmingustb/PTPSim
"""


import matplotlib.pyplot as plt


# optimization parameter
max_iter = 100
cost_th = 0.01
# cost_th = 0.1
show_animation = False

station_time_animation = 1      # 是否画出s-t图


def speedprofile(start, end):

    if station_time_animation:
        x = [start[0], end[0]]
        y = [start[1], end[1]]
        plt.plot(x, y, color='r', linewidth=0.3, alpha=0.8)
        #   plt.show()


if __name__ == '__main__':
    start = [0, 0, 0, 0]
    end = [2, 2, 0]
    speedprofile(start, end)
    plt.show()

