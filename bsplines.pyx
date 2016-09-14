"""Fast cubic basis splines in Python"""
from libc.stdio cimport sprintf
import numpy as np
cimport numpy as np

__version__ = "0.1.0"
__all__ = ["b3", "Spline1D"]


cdef extern from "bspl.h":
    double bs_b3(double x, int i, double *t)
    double bs_db3(double x, int i, double *t)
    double bs_ddb3(double x, int i, double *t)

    ctypedef struct bs_array:
        double *data
        int length
        int stride

    ctypedef struct bs_spline1d:
        double *knots
        double *coeffs
        int n

    ctypedef struct bs_bc:
        int deriv
        double value

    ctypedef struct bs_bcs:
        bs_bc left
        bs_bc right

    bs_spline1d* bs_create_spline1d(bs_array x, bs_array y, bs_bcs bcs)
    double bs_eval_spline1d(bs_spline1d *spline, double x)
    void bs_free_spline1d(bs_spline1d *spline)


def b3(double x, int i, double[:] t):    
    return bs_b3(x, i, &t[0]);

def db3(double x, int i, double[:] t):    
    return bs_db3(x, i, &t[0]);

def ddb3(double x, int i, double[:] t):    
    return bs_ddb3(x, i, &t[0]);


cdef class Spline1D:
    """
    Spline1D(x, y, bc='natural')

    One dimensional cubic basis spline.
    
    Parameters
    ----------
    x : 1-d `~numpy.ndarray`
    y : 1-d `~numpy.ndarray`
    """

    cdef bs_spline1d *ptr   # pointer to c struct
    
    def __cinit__(self, np.ndarray x, np.ndarray y, bcs=((2, 0.0), (2, 0.0))):

        # require 1-d arrays
        if (x.ndim != 1 or y.ndim != 1):
            raise ValueError("x and y must be 1-d arrays")

        # convert to double arrays if needed
        cdef double[:] x_ = np.require(x, dtype=np.dtype(np.double))
        cdef double[:] y_ = np.require(y, dtype=np.dtype(np.double))

        cdef bs_array x_arr = bs_array(&x_[0], x_.shape[0],
                                           x_.strides[0]//sizeof(double))
        cdef bs_array y_arr = bs_array(&y_[0], y_.shape[0],
                                           y_.strides[0]//sizeof(double))

        # parse boundary conditions
        cdef bs_bcs bcs_ = bs_bcs(bs_bc(bcs[0][0], bcs[0][1]),
                                      bs_bc(bcs[1][0], bcs[1][1]))
        
        self.ptr = bs_create_spline1d(x_arr, y_arr, bcs_)


    #def __init__(self, np.ndarray x, np.ndarray y):
    #    pass


    def __dealloc__(self):
        if self.ptr is not NULL:
            bs_free_spline1d(self.ptr)

    def __call__(self, double x):
        return bs_eval_spline1d(self.ptr, x)

    def coefficients(self):
        cdef double[:] view = <double[:(self.ptr.n+2)]>(self.ptr.coeffs)
        return np.array(view)  # copy

