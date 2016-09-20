#!/usr/bin/env python

import os
from time import time
from collections import OrderedDict
import pickle

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

from bsplines import Spline1D


# -----------------------------------------------------------------------------
# utilities
# -----------------------------------------------------------------------------

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
        if not np.all(allsizes[0] == allsizes[i]):
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


def save_results(results, fname):
    with open(fname, 'wb') as f:
        pickle.dump(results, f)


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

if __name__ == "__main__":

    # results stored here in pickle files for plotting in docs
    os.makedirs("benchmarks", exist_ok=True)


    results = OrderedDict([
        ('bsplines.Spline1D', benchmark_creation_1d(Spline1D, {})),
        ('SciPy UnivariateSpline',
         benchmark_creation_1d(InterpolatedUnivariateSpline, {'ext': 3, 'k': 3}))
    ])
    print_results(results, "1-d spline creation", "knots")
    save_results(results, os.path.join("benchmarks", "1d_create.pik"))


    results = OrderedDict([
        ('bsplines.Spline1D', benchmark_eval_1d(Spline1D, {})),
        ('SciPy UnivariateSpline',
         benchmark_eval_1d(InterpolatedUnivariateSpline, {'ext': 3, 'k': 3}))
    ])

    print_results(results, "1-d spline evaluation", "points")
    save_results(results, os.path.join("benchmarks", "1d_eval.pik"))
