#!/usr/bin/env py.test
import numpy as np
from numpy.testing import assert_allclose
import pytest

from bsplines import Spline1D, USpline1D, Spline2D, DomainError

def test_spline1d_notaknot():
    x = np.array([-1., -0.5, 0.5, 2.2, 3.8, 4.])
    y = x**3

    spline = Spline1D(x, y, bcs='notaknot')

    xtest = np.linspace(-1., 4., 101)
    ytest_true = xtest**3
    ytest = spline(xtest)

    assert_allclose(ytest, ytest_true, atol=1e-14)

def test_spline1d_minimal():
    """Test minimal number of knots."""
    x = np.array([-0.5, 0.5, 2.2])
    y = x**3

    spline = Spline1D(x, y, bcs=('notaknot', ('deriv2', 6.*x[-1])))

    xp = np.linspace(-0.5, 2.2, 101)
    yp = xp**3
    assert_allclose(spline(xp), yp, atol=1e-14)


def test_uspline1d_notaknot():
    x = np.array([-1., 0., 1., 2., 3., 4.])
    y = x**3

    s = USpline1D((x[0], x[-1]), y, bcs='notaknot')

    xp = np.linspace(-1., 4., 101)
    assert_allclose(s(xp), xp**3, atol=1e-14)


def test_spline1d_cubic():
    x = np.array([-1., -0.5, 0.5, 2.2, 3.8, 4.])
    y = x**3

    # 2nd deriv boundary conditions from true function
    bcs = (("deriv2", 6. * x[0]), ("deriv2", 6. * x[-1]))
    s = Spline1D(x, y, bcs=bcs)

    xp = np.linspace(-1., 4., 101)
    assert_allclose(s(xp), xp**3, atol=1e-14)


def test_1d_input_errors():
    """Test some errors in spline construction"""

    # raise an error about not monotonically increasing x array
    x = [1., 0.999, 3., 4., 5.]
    y = [2., 2., 2., 2., 2.]
    with pytest.raises(ValueError) as excinfo:
        Spline1D(x, y)
    assert 'monotonic' in excinfo.value.args[0]


def test_1d_extend_raise():
    """Test that raising a domain error works."""

    x = [1., 2., 3., 4., 5.]
    y = [2., 2., 2., 2., 2.]
    s = Spline1D(x, y, extend='raise')
    with pytest.raises(DomainError):
        s(np.array([0.9, 1.1]))


def test_spline2d():
    data = np.ones((5, 4))
    s = Spline2D([1., 2., 3., 4., 5.], [1., 2., 3., 5.], data)
    xp = np.array([1.5, 2.5])
    yp = np.array([1.5, 2.5])
    assert_allclose(s(xp, yp), np.ones((len(xp), len(yp))))
