import numpy

x3 = [[1, 2, 3], [2, 3, 4]]
x1 = numpy.array([[1], [2], [3], [4]])
x2 = numpy.array(x1)
x4 = numpy.array(x3)
print("x4:", x4)
x8 = numpy.array(x4)
print("x8:", x8)


x5 = numpy.array([1, 2, 3])
x6 = numpy.array(x1)
print(x2)
print(x4+x5)
print(x6)
print(x1[1, 0])

x = {1: 18, 2: 99}
haha3 = max(x, key=lambda a: a)
del x[1]
print(haha3)
print(x)

dic1 = dict()
dic1[0] = 0
print(dic1)

x5 = numpy.array([[1], [2], [3]]).reshape(1, 3)
print('x5:', x5[0])