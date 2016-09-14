"""Fast cubic basis splines in Python"""
from libc.stdio cimport sprintf
import numpy as np
cimport numpy as np

__version__ = "0.1.0"
__all__ = ["b3", "Spline1D"]


cdef extern from "bspl.h":
    double bspl_b3(double x, int i, double *t)
    double bspl_db3(double x, int i, double *t)
    double bspl_ddb3(double x, int i, double *t)

    ctypedef struct bspl_array:
        double *data
        int length
        int stride

    ctypedef struct bspl_spline1d:
        double *knots
        double *coeffs
        int n

    ctypedef struct bspl_bc:
        int deriv
        double value

    ctypedef struct bspl_bcs:
        bspl_bc left
        bspl_bc right

    bspl_spline1d* bspl_create_spline1d(bspl_array x, bspl_array y, bspl_bcs bcs)
    double bspl_eval_spline1d(bspl_spline1d *spline, double x)
    void bspl_free_spline1d(bspl_spline1d *spline)


def b3(double x, int i, double[:] t):    
    return bspl_b3(x, i, &t[0]);

def db3(double x, int i, double[:] t):    
    return bspl_db3(x, i, &t[0]);

def ddb3(double x, int i, double[:] t):    
    return bspl_ddb3(x, i, &t[0]);


cdef class Spline1D:
    """
    Spline1D(x, y, bc='natural')

    One dimensional cubic basis spline.
    
    Parameters
    ----------
    x : 1-d `~numpy.ndarray`
    y : 1-d `~numpy.ndarray`
    """

    cdef bspl_spline1d *ptr   # pointer to c struct
    
    def __cinit__(self, np.ndarray x, np.ndarray y, bcs=((2, 0.0), (2, 0.0))):

        # require 1-d arrays
        if (x.ndim != 1 or y.ndim != 1):
            raise ValueError("x and y must be 1-d arrays")

        # convert to double arrays if needed
        cdef double[:] x_ = np.require(x, dtype=np.dtype(np.double))
        cdef double[:] y_ = np.require(y, dtype=np.dtype(np.double))

        cdef bspl_array x_arr = bspl_array(&x_[0], x_.shape[0],
                                           x_.strides[0]//sizeof(double))
        cdef bspl_array y_arr = bspl_array(&y_[0], y_.shape[0],
                                           y_.strides[0]//sizeof(double))

        # parse boundary conditions
        cdef bspl_bcs bcs_ = bspl_bcs(bspl_bc(bcs[0][0], bcs[0][1]),
                                      bspl_bc(bcs[1][0], bcs[1][1]))
        
        self.ptr = bspl_create_spline1d(x_arr, y_arr, bcs_)


    #def __init__(self, np.ndarray x, np.ndarray y):
    #    pass


    def __dealloc__(self):
        if self.ptr is not NULL:
            bspl_free_spline1d(self.ptr)

    def __call__(self, double x):
        return bspl_eval_spline1d(self.ptr, x)

    def coefficients(self):
        cdef double[:] view = <double[:(self.ptr.n+2)]>(self.ptr.coeffs)
        return np.array(view)  # copy

