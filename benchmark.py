#!/usr/bin/env python

import os
from time import time

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt

from bsplines import Spline1D

plt.style.use('bmh')

os.makedirs("benchmarks", exist_ok=True)

def plot_results(results, title, unit, fname):
    fig, ax = plt.subplots()

    for key, result in results.items():
        throughputs = result['sizes'] / (1e6 * result['times'])
        ax.semilogx(result['sizes'], throughputs, ls='-', label=key)

    ax.set_xlabel(unit)
    ax.set_ylabel("throughput ({} / us)".format(unit))
    ax.set_title(title)
    ax.legend(loc='upper left')
    fig.savefig(os.path.join("benchmarks", fname))
    plt.close(fig)


# --------------------------------------------------------------------------
# spline creation
# --------------------------------------------------------------------------

def benchmark_creation_1d(cls, kwargs):
    sizes = np.array([10, 30, 100, 1000, 10000, 100000])
    nloops = np.array([max(100000 // sz, 1) for sz in sizes])
    times = np.empty_like(sizes, dtype=np.float64)
    for i, n in enumerate(sizes):
        x = np.linspace(0., float(n), n)
        y = np.sin(x)

        t0 = time()
        for _ in range(nloops[i]):
            cls(x, y, **kwargs)
        times[i] = (time() - t0) / nloops[i]

    return {'sizes': sizes,
            'nloops': nloops,
            'times': times}

results = {
    'bsplines.Spline1D': benchmark_creation_1d(Spline1D, {}),
    'SciPy UnivariateSpline':
    benchmark_creation_1d(InterpolatedUnivariateSpline, {'ext': 3, 'k': 3})
}
plot_results(results, "1-d spline creation", "knots", "1d_create.png")


# ----------------------------------------------------------------------------
# spline evaluation
# ----------------------------------------------------------------------------

def benchmark_eval_1d(cls, kwargs):
    sizes = np.array([10, 30, 100, 1000, 10000, 100000])
    nloops = np.array([max(100000 // sz, 1) for sz in sizes])
    times = np.empty_like(sizes, dtype=np.float64)
    for i, n in enumerate(sizes):
        x = np.linspace(0., 1000., 1000)
        y = np.sin(x)
        s = cls(x, y, **kwargs)

        xtest = np.linspace(0., 1000., n)
        t0 = time()
        for _ in range(nloops[i]):
            s(xtest)
        times[i] = (time() - t0) / nloops[i]

    return {'sizes': sizes,
            'nloops': nloops,
            'times': times}

results = {
    'bsplines.Spline1D': benchmark_eval_1d(Spline1D, {}),
    'SciPy UnivariateSpline':
    benchmark_eval_1d(InterpolatedUnivariateSpline, {'ext': 3, 'k': 3})
}
plot_results(results, "1-d spline evaluation", "points", "1d_eval.png")

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
