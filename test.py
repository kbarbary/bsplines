#!/usr/bin/env py.test
import numpy as np
from numpy.testing import assert_allclose
import pytest

from bsplines import Spline1D, USpline1D, DomainError

def test_spline1d_notaknot():
    x = np.array([-1., -0.5, 0.5, 2.2, 3.8, 4.])
    y = x**3

    spline = Spline1D(x, y, bcs='notaknot')

    xtest = np.linspace(-1., 4., 101)
    ytest_true = xtest**3
    ytest = spline(xtest)

    assert_allclose(ytest, ytest_true, atol=1e-14)

def test_uspline1d_notaknot():
    x = np.array([-1., 0., 1., 2., 3., 4.])+0.1
    y = x**3

    spline = USpline1D((x[0], x[-1]), y,
                       bcs='notaknot')
    spline2 = Spline1D(x, y, bcs='notaknot')
    print(spline2.coefficients)

    xp = np.linspace(-1., 4., 101)
    yp = xp**3

    assert_allclose(spline(xp), yp, atol=1e-14)


def test_spline1d_cubic():
    x = np.array([-1., -0.5, 0.5, 2.2, 3.8, 4.])
    y = x**3

    # 2nd deriv boundary conditions from true function
    bcs = (("deriv2", 6. * x[0]), ("deriv2", 6. * x[-1]))
    spline = Spline1D(x, y, bcs=bcs)

    xtest = np.linspace(-1., 4., 101)
    ytest_true = xtest**3
    ytest = spline(xtest)

    #from matplotlib import pyplot as plt
    #plt.plot(x, y, ls='None', marker='o')
    #plt.plot(xtest, ytest, marker='None', ls='-')
    #plt.show()
    
    assert_allclose(ytest, ytest_true, atol=1e-14)


def test_1d_input_errors():
    """Test some errors in spline construction"""

    # raise an error about not monotonically increasing x array
    x = [1., 0.999, 3., 4., 5.]
    y = [2., 2., 2., 2., 2.]
    with pytest.raises(ValueError) as excinfo:
        Spline1D(x, y)
    assert 'monotonic' in excinfo.value.args[0]


def test_1d_extrapolation():
    """Test that raising a domain error works."""

    x = [1., 2., 3., 4., 5.]
    y = [2., 2., 2., 2., 2.]
    s = Spline1D(x, y, extend='raise')
    with pytest.raises(DomainError):
        s(np.array([0.9, 1.1]))

