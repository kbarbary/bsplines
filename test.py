#!/usr/bin/env py.test
import numpy as np
from numpy.testing import assert_allclose

import bsplines

def test_spline1d_cubic():
    x = np.array([-1., -0.5, 0.5, 2.2, 3.8, 4.])
    y = x**3

    # 2nd deriv boundary conditions from true function
    bcs = (("deriv2", 6. * x[0]), ("deriv2", 6. * x[-1]))
    spline = bsplines.Spline1D(x, y, bcs=bcs)

    xtest = np.linspace(-1., 4., 101)
    ytest_true = xtest**3
    ytest = spline(xtest)

    #from matplotlib import pyplot as plt
    #plt.plot(x, y, ls='None', marker='o')
    #plt.plot(xtest, ytest, marker='None', ls='-')
    #plt.show()
    
    assert_allclose(ytest, ytest_true, atol=1e-14)

#test_spline1d_cubic()
