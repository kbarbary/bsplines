#!/usr/bin/env py.test
import numpy as np
from numpy.testing import assert_allclose

import bsplines


# slow recursive implementation of basis functions to compare to:
def b(k, x, i, t):
    """B-spline of order k at the i-th index."""
    if k == 0:
        return 1.0 if (t[i] <= x < t[i+1]) else 0.0

    return ((x - t[i]) / (t[i+k] - t[i]) * b(k-1, x, i, t) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * b(k-1, x, i+1, t))


def db(k, x, i, t):
    """First derivative of B-spline of order k at the i-th index."""

    if k == 0:
        return 0.0

    return ((x - t[i]) / (t[i+k] - t[i]) * db(k-1, x, i, t) +
            1.0 / (t[i+k] - t[i]) * b(k-1, x, i, t) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * db(k-1, x, i+1, t) +
            -1.0 / (t[i+k+1] - t[i+1]) * b(k-1, x, i+1, t))

def ddb(k, x, i, t):
    """Second derivative of B-spline of order k at the i-th index."""

    if k == 0:
        return 0.0

    return ((x - t[i]) / (t[i+k] - t[i]) * ddb(k-1, x, i, t) +
            2.0 / (t[i+k] - t[i]) * db(k-1, x, i, t) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * ddb(k-1, x, i+1, t) +
            -2.0 / (t[i+k+1] - t[i+1]) * db(k-1, x, i+1, t))


def test_b3():
    """Test that basis function matches naive implementation"""
    t = np.array([0., 1., 2., 3., 4.])  # knots
    x = np.linspace(-1., 5., 100)
    y = np.array([bsplines.b3(xi, 0, t) for xi in x])
    assert_allclose(np.array([bsplines.b3(xi, 0, t) for xi in x]),
                    np.array([b(3, xi, 0, t) for xi in x]))


def test_db3():
    """Test that basis function matches naive implementation"""
    t = np.array([0., 1., 2., 3., 4.])  # knots
    x = np.linspace(-1., 5., 100)
    y = np.array([bsplines.db3(xi, 0, t) for xi in x])
    y_true = np.array([db(3, xi, 0, t) for xi in x])
    assert_allclose(y, y_true)

def test_ddb3():
    """Test that basis function matches naive implementation"""
    t = np.array([0., 1., 2., 3., 4.])  # knots
    x = np.linspace(-1., 5., 100)
    y = np.array([bsplines.ddb3(xi, 0, t) for xi in x])
    y_true = np.array([ddb(3, xi, 0, t) for xi in x])
    assert_allclose(y, y_true)

def test_spline1d_cubic():
    x = np.array([-1., -0.5, 0.5, 2.2, 3.8, 4.])
    y = x**3

    # 2nd deriv boundary conditions from true function
    bcs = (("deriv2", 6. * x[0]), ("deriv2", 6. * x[-1]))
    spline = bsplines.Spline1D(x, y, bc=bcs)

    xtest = np.linspace(-1., 4., 101)
    ytest_true = xtest**3
    ytest = spline(xtest)

    #from matplotlib import pyplot as plt
    #plt.plot(x, y, ls='None', marker='o')
    #plt.plot(xtest, ytest, marker='None', ls='-')
    #plt.show()
    
    for i in range(len(ytest)):
        print(xtest[i], ytest_true[i], ytest[i])
    assert_allclose(ytest, ytest_true)

#test_spline1d_cubic()
