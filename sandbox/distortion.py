import scipy.optimize as optimize
import matplotlib.pyplot as pyplot
import numpy


X = numpy.array([8.63, 8.74, 8.85, 8.96, 9.06, 9.17, 9.28, 9.40])
Y = numpy.array([11.91, 12.02, 12.13, 12.23, 12.34, 12.45, 12.55, 12.67])

def line(m, b, n):
    retval = []
    for i in range(n):
        retval.append(m*n+b)
    return numpy.array(retval)

def eval(m, b, y):
    return y-line(m, b, len(y))
        

ans = optimize.leastsq(eval, X, args=(m, b))
