#!/usr/bin/env python

from time import time

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from bsplines import Spline1D


x = np.linspace(0., 1000., 1000)
xtest = np.linspace(0., 1000., 100000)
y = np.sin(x)

t0 = time()
s = Spline1D(x, y)
t = time() - t0
print("create: {} s".format(t))

t0 = time()
s_scipy = InterpolatedUnivariateSpline(x, y, ext=3, k=3)
t = time() - t0
print("create scipy: {} s".format(t))

t0 = time()
s(xtest)
t = time() - t0
print("eval: {} s".format(t))

t0 = time()
s_scipy(xtest)
t = time() - t0
print("eval scipy: {} s".format(t))
