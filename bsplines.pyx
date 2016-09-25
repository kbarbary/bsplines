"""Fast cubic basis splines in Python"""
from libc.stdio cimport sprintf
import numpy as np
cimport numpy as np

__version__ = "0.1.0"
__all__ = ["Spline1D"]


cdef extern from "bs.h":
    ctypedef enum bs_errorcode:
        BS_OK           = 0
        BS_OUTOFMEMORY  = 1
        BS_DOMAINERROR  = 2
        BS_NOTMONOTONIC = 3
        BS_LENGTHMISMATCH = 4

    ctypedef struct bs_array:
        double *data
        int length
        int stride

    ctypedef enum bs_bctype:
        BS_DERIV1
        BS_DERIV2
        
    ctypedef struct bs_bc:
        bs_bctype type
        double value

    ctypedef struct bs_bcs:
        bs_bc left
        bs_bc right

    ctypedef enum bs_exttype:
        BS_EXTRAPOLATE
        BS_CONSTANT
        BS_RAISE

    ctypedef struct bs_ext:
        bs_exttype type
        double value

    ctypedef struct bs_exts:
        bs_ext left
        bs_ext right
        
    ctypedef struct bs_spline1d:
        double *knots
        double *coeffs
        int n

    bs_errorcode bs_spline1d_create(bs_array x, bs_array y, bs_bcs bcs, bs_exts exts, bs_spline1d **out)
    bs_errorcode bs_spline1d_eval(bs_spline1d *spline, bs_array x, bs_array out)
    void bs_spline1d_free(bs_spline1d *spline)


class DomainError(Exception):
    pass

#------------------------------------------------------------------------------
# helpers

cdef int assert_ok(bs_errorcode code) except -1:
    """raise an exception of the error code is not OK"""
    if code == BS_OK:
        return 0
    elif code == BS_OUTOFMEMORY:
        raise MemoryError()
    elif code == BS_DOMAINERROR:
        raise DomainError()
    elif code == BS_NOTMONOTONIC:
        raise ValueError("input array(s) not monotonically increasing")
    elif code == BS_LENGTHMISMATCH:
        raise ValueError("input arrays have different lengths")
    else:
        raise Exception("unrecognized error")


cdef int parse_boundary_conditions(pybc, bs_bcs *out) except -1:
    """Parse boundary conditions"""
    cdef bs_bctype lefttype, righttype
    
    if pybc == "natural":
        out.left = out.right = bs_bc(BS_DERIV2, 0.0)
    elif len(pybc) == 2:
        left, right = pybc

        if left[0] == "deriv1":
            lefttype = BS_DERIV1
        else:
            lefttype = BS_DERIV2
        out.left = bs_bc(lefttype, left[1])

        if right[0] == "deriv1":
            righttype = BS_DERIV1
        else:
            righttype = BS_DERIV2
        out.right = bs_bc(righttype, right[1])

    else:
        raise ValueError("unrecognized boundary condition: " + repr(pybc))


# -----------------------------------------------------------------------------
# 1-d spline

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
    
    def __cinit__(self, np.ndarray x, np.ndarray y, bc="natural"):

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
        cdef bs_bcs bcs
        parse_boundary_conditions(bc, &bcs)

        # extensions
        cdef bs_exts exts = bs_exts(bs_ext(BS_CONSTANT, y_[0]),
                                    bs_ext(BS_CONSTANT, y_[-1]))

        cdef bs_errorcode code
        code = bs_spline1d_create(x_arr, y_arr, bcs, exts, &self.ptr)
        assert_ok(code)


    #def __init__(self, np.ndarray x, np.ndarray y):
    #    pass


    def __dealloc__(self):
        bs_spline1d_free(self.ptr)

    def __call__(self, double[:] x):
        cdef double[:] outview
        
        out = np.empty(len(x))
        outview = out

        cdef bs_array x_arr = bs_array(&x[0], x.shape[0],
                                       x.strides[0]//sizeof(double))
        cdef bs_array out_arr = bs_array(&outview[0], outview.shape[0],
                                         outview.strides[0]//sizeof(double))

        cdef bs_errorcode code = bs_spline1d_eval(self.ptr, x_arr, out_arr);
        assert_ok(code)

        return out

    def coefficients(self):
        cdef double[:] view = <double[:(self.ptr.n+2)]>(self.ptr.coeffs)
        return np.array(view)  # copy

