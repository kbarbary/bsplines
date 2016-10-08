#!/usr/bin/env python

import os
from time import time
from collections import OrderedDict
import json

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline, CubicSpline as SciPyCubicSpline

from bsplines import Spline1D, USpline1D, Spline2D


# -----------------------------------------------------------------------------
# utilities
# -----------------------------------------------------------------------------

def timeit(f, args, kwargs, target=0.05):
    """Time execution of f(*args, **kwargs).

    Returns best of 3 run times as (time_per_call, nloops)."""

    # determine number of loops to run.
    t = 0.0
    nloops = 1
    while True:
        t0 = time()
        for _ in range(nloops):
            f(*args, **kwargs)
        t = time() - t0
        if t > target or nloops >= 10**9:
            break
        nloops *= 10

    # production runs
    times = [0., 0., 0.]
    for i in range(3):
        t0 = time()
        for _ in range(nloops):
            f(*args, **kwargs)
        times[i] = time() - t0
        
    return min(times) / nloops, nloops
    

def print_results(results, title, unit):
    """
    Parameters
    ----------
    results : dict
        Dictionary where key is spline name and value is a dictionary of
        timing results. The dictionary contains keys ``sizes``, ``times``,
        which are both arrays.
    title : str
        Axes title.
    unit : str
        Unit of ``sizes`` (e.g., knots or points).
    """

    # check that all `sizes` arrays are the same.
    allsizes = list(result['sizes'] for result in results.values())
    for i in range(1, len(allsizes)):
        if not allsizes[0] == allsizes[i]:
            raise ValueError("Results must have same sizes for printing")

    # header
    print("\n" + title + " (ms)")
    print("{:10s}".format(unit), end='')
    for key in results.keys():
        print("{:25s}".format(key), end='')
    print("\n" + "-" * 60)
    
    sizes = allsizes[0]
    for i in range(len(sizes)):
        print("{:8d}  ".format(sizes[i]), end='')
        for key in results.keys():
            print("{:10.6f}".format(1000 * results[key]['times'][i]) + " "*15,
                  end='')
        print()

    print("-"*60)


def save_results(results, title, unit, fname):
    with open(fname, 'w') as f:
        json.dump({'title': title, 'unit': unit, 'results': results}, f)

# -----------------------------------------------------------------------------
# benchmarks
# -----------------------------------------------------------------------------
        
def benchmark_creation_1d(cls, kwargs):
    sizes = np.array([10, 30, 100, 1000, 10000, 100000])
    nloops = np.empty_like(sizes)
    times = np.empty_like(sizes, dtype=np.float64)
    for i, n in enumerate(sizes):
        x = np.linspace(0., float(n), n)
        y = np.sin(x)

        times[i], nloops[i] = timeit(cls, (x, y), kwargs)

    return {'sizes': sizes.tolist(),
            'nloops': nloops.tolist(),
            'times': times.tolist()}


def benchmark_eval_1d(cls, kwargs):
    sizes = np.array([10, 30, 100, 1000, 10000, 100000])
    nloops = np.empty_like(sizes)
    times = np.empty_like(sizes, dtype=np.float64)
    for i, n in enumerate(sizes):
        x = np.linspace(0., 1000., 1000)
        y = np.sin(x)
        s = cls(x, y, **kwargs)

        xp = np.linspace(0., 1000., n)
        times[i], nloops[i] = timeit(s, (xp,), {})

    return {'sizes': sizes.tolist(),
            'nloops': nloops.tolist(),
            'times': times.tolist()}


def benchmark_create_2d(cls, kwargs):
    sizes = np.array([5, 10, 30, 100, 300, 1000])
    nloops = np.empty_like(sizes)
    times = np.empty_like(sizes, dtype=np.float64)
    for i, n in enumerate(sizes):
        x = np.linspace(0., float(n), n)
        y = np.linspace(0., float(n), n)
        z = np.sin(x) + np.cos(y).reshape((n, 1))

        times[i], nloops[i] = timeit(cls, (x, y, z), kwargs)

    return {'sizes': sizes.tolist(),
            'nloops': nloops.tolist(),
            'times': times.tolist()}

def benchmark_eval_2d(cls, kwargs):
    nknots = 100
    x = np.linspace(0., float(nknots), nknots)
    y = np.linspace(0., float(nknots), nknots)
    z = np.sin(x) + np.cos(y).reshape((nknots, 1))
    s = cls(x, y, z, **kwargs)
    
    sizes = np.array([3, 10, 30, 100, 300, 1000])
    nloops = np.empty_like(sizes)
    times = np.empty_like(sizes, dtype=np.float64)
    for i, n in enumerate(sizes):
        xp = np.linspace(0., float(n), n)
        yp = np.linspace(0., float(n), n)

        times[i], nloops[i] = timeit(s, (xp, yp), {})

    return {'sizes': sizes.tolist(),
            'nloops': nloops.tolist(),
            'times': times.tolist()}


if __name__ == "__main__":

    # results stored here in pickle files for plotting in docs
    os.makedirs("benchmarks", exist_ok=True)

    # Spline 1d creation
    results = OrderedDict([
        ('bsplines.USpline1D', benchmark_creation_1d(USpline1D, {})),
        ('bsplines.Spline1D', benchmark_creation_1d(Spline1D, {})),
         ('SciPy CubicSpline', benchmark_creation_1d(SciPyCubicSpline, {})),
        ('SciPy UnivariateSpline',
         benchmark_creation_1d(InterpolatedUnivariateSpline, {'ext': 3, 'k': 3}))
    ])
    print_results(results, "1-d spline creation", "knots")
    save_results(results, "1-d spline creation", "knots",
                 os.path.join("benchmarks", "1d_create.json"))

    # spline 1d evaluation
    results = OrderedDict([
        ('bsplines.USpline1D', benchmark_eval_1d(USpline1D, {})),
        ('bsplines.Spline1D', benchmark_eval_1d(Spline1D, {})),
        ('SciPy CubicSpline', benchmark_eval_1d(SciPyCubicSpline, {})),
        ('SciPy UnivariateSpline',
         benchmark_eval_1d(InterpolatedUnivariateSpline, {'ext': 3, 'k': 3}))
    ])

    print_results(results, "1-d spline evaluation", "points")
    save_results(results, "1-d spline evaluation", "points",
                 os.path.join("benchmarks", "1d_eval.json"))


    # 2-d creation
    results = OrderedDict([
        ('bsplines.Spline2D', benchmark_create_2d(Spline2D, {})),
        ('SciPy RectBivariateSpline',
         benchmark_create_2d(RectBivariateSpline, {'kx': 3, 'ky': 3}))
    ])

    print_results(results, "2-d spline creation", "knots")


    # 2-d eval
    results = OrderedDict([
        ('bsplines.Spline2D', benchmark_eval_2d(Spline2D, {})),
        ('SciPy RectBivariateSpline',
         benchmark_eval_2d(RectBivariateSpline, {'kx': 3, 'ky': 3}))
    ])

    print_results(results, "2-d spline evaluation", "points")
