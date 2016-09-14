#!/usr/bin/env python

import os
from time import time

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt

from bsplines import Spline1D

os.makedirs("benchmarks", exist_ok=True)

# spline creation

sizes = np.array([10, 100, 1000, 10000, 100000])
t = np.empty_like(sizes, dtype=np.float64)
t_scipy = np.empty_like(sizes, dtype=np.float64)
for i, n in enumerate(sizes):
    x = np.linspace(0., float(n), n)
    y = np.sin(x)

    nloops = max(10000 // n, 1)

    t0 = time()
    for _ in range(nloops):
        Spline1D(x, y)
    t[i] = (time() - t0) / nloops

    t0 = time()
    for _ in range(nloops):
        InterpolatedUnivariateSpline(x, y, ext=3, k=3)
    t_scipy[i] = (time() - t0) / nloops

fig, ax = plt.subplots()
ax.semilogx(sizes, t / sizes, ls='-', label='Spline1D')
ax.semilogx(sizes, t_scipy / sizes, ls='-', label='UnivariateSpline')
ax.set_xlabel("knots")
ax.set_ylabel("time (s) / knot")
ax.set_title("Spline creation")
ax.legend()
fig.savefig(os.path.join("benchmarks", "1d_create.png"))
plt.close(fig)

    
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
