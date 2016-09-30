"""Fast cubic basis splines in Python"""
from libc.math cimport fabs
import numpy as np
cimport numpy as np

__version__ = "0.1.0"
__all__ = ["Spline1D", "USpline1D", "DomainError"]


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

    ctypedef struct bs_range:
        double min
        double max

    ctypedef enum bs_bctype:
        BS_DERIV1
        BS_DERIV2
        BS_NOTAKNOT
        
    ctypedef struct bs_bc:
        bs_bctype type
        double value

    ctypedef struct bs_bcs:
        bs_bc left
        bs_bc right

    ctypedef enum bs_exttype:
        BS_EXTRAPOLATE
        BS_CONSTANT
        BS_VALUE
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

    ctypedef struct bs_uspline1d:
        bs_range x
        double didx
        double *coeffs
        int n
        bs_exts exts

    bs_errorcode bs_uspline1d_create(bs_range x, bs_array y,
                                     bs_bcs bcs, bs_exts exts,
                                     bs_uspline1d **out)
    bs_errorcode bs_uspline1d_eval(bs_uspline1d *spline, bs_array x,
                                   bs_array out)
    void bs_uspline1d_free(bs_uspline1d *spline)

    
class DomainError(Exception):
    """Raised when spline input(s) are outside spline boundaries."""
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


cdef inline bs_array to_bs_array(double[:] x):
    return bs_array(&x[0], x.shape[0], x.strides[0]//sizeof(double))


cdef int try_parse_boundary_condition(pybc, bs_bc *parsed_bc):
    """Try to parse a single Python value as a single boundary condition

    Returns 1 if parsing is sucessful, else 0.
    """

    if pybc == "notaknot":
        parsed_bc.type = BS_NOTAKNOT
        parsed_bc.value = 0.0  # not used
        return 1
    if pybc == "natural":
        parsed_bc.type = BS_DERIV2
        parsed_bc.value = 0.0
        return 1
    elif pybc == "flat":
        parsed_bc.type = BS_DERIV1
        parsed_bc.value = 0.0
        return 1
    elif len(pybc) == 2:
        if pybc[0] == "deriv1":
            parsed_bc.type = BS_DERIV1
            parsed_bc.value = pybc[1]  # need to check if numeric?
            return 1
        elif pybc[0] == "deriv2":
            parsed_bc.type = BS_DERIV2
            parsed_bc.value = pybc[1]  # need to check if numeric?
            return 1

    return 0

    
cdef int try_parse_boundary_conditions(pybcs, bs_bcs *parsed_bcs):
    """Try to parse boundary conditions. Returns 1 if sucessful, else 0."""

    # try to parse the expression as a boundary condition: if it works,
    # apply to both right and left.
    if try_parse_boundary_condition(pybcs, &parsed_bcs.left):
        parsed_bcs.right = parsed_bcs.left
        return 1

    # next, try to parse the expression as a left and right bc separately.
    if (len(pybcs) == 2 and
        try_parse_boundary_condition(pybcs[0], &parsed_bcs.left) and
        try_parse_boundary_condition(pybcs[1], &parsed_bcs.right)):
        return 1

    return 0


cdef int try_parse_extend(pyextend, bs_ext *parsed_ext):

    if pyextend == 'extrapolate':
        parsed_ext.type = BS_EXTRAPOLATE
        return 1
    elif pyextend == 'constant':
        parsed_ext.type = BS_CONSTANT
        return 1
    elif pyextend == 'raise':
        parsed_ext.type = BS_RAISE
        return 1
    else:
        try:
            value = float(pyextend)
        except:
            return 0
        parsed_ext.type = BS_VALUE
        parsed_ext.value = value
        return 1


cdef int try_parse_extends(pyextends, bs_exts *parsed_exts):

     # Try to parse the expression as a bs_ext. If it works, apply to both
     # right and left.
     if try_parse_extend(pyextends, &parsed_exts.left):
         parsed_exts.right = parsed_exts.left
         return 1

     if (len(pyextends) == 2 and
         try_parse_extend(pyextends[0], &parsed_exts.left) and
         try_parse_extend(pyextends[0], &parsed_exts.right)):
         return 1

     return 0


cdef int is_uniform(double[:] x):
    cdef int i
    cdef double dx = x[1] - x[0]
    cdef int ok = 1

    for i in range(2, len(x)):
        ok &= (fabs((x[i] - x[i-1]) - dx) < 1e-9 * dx)

    return ok

# -----------------------------------------------------------------------------
# 1-d splines

cdef class Spline1D:
    """
    Spline1D(x, y, bcs='notaknot', extend='constant')

    One dimensional cubic basis spline.
    
    Parameters
    ----------
    x : `numpy.ndarray` (1-d)
        Abscissa values.

    y : `numpy.ndarray` (1-d)
        Ordinate values to interpolate through.

    bcs : str or tuple
        One of:

        - ``'notaknot'``: The "not a knot" condition: the third derivative is
          enforced to be continuous at the first interior knot point.
        - ``'natural'``: set second derivative to zero at boundary.
        - ``'flat'``: set first derivative to zero at boundary.
        - ``('deriv1', value)``: set first derivative to ``value`` at boundary.
        - ``('deriv2', value)``: set second derivative to ``value`` at
          boundary.

        Can also be a 2-tuple of the above values specifying the left and
        right boundary conditions separately.

    extend : str or float or tuple
        Controls what happens when a value outside the spline domain is
        passed.

        - ``'extrapolate'`` Extrapolate past the last knot.
        - ``'constant'`` Return the value at the outermost knot.
        - ``'raise'`` Raise a ``bspline.DomainError``.
        - numeric value: return this value.

        Can also be a 2-tuple of the above values specifying the left and
        right behavior separately.


    Examples
    --------
    
    Specifying boundary conditions::

        Spline1D(x, y, bcs=('deriv1', 4.0)) # set first derivative to 4 at
                                            # both ends

        Spline1D(x, y, bcs=('natural', ('deriv2', 0)))

    Specifying beavior outside domain::

        Spline1D(x, y, extend='extrapolate')
        Spline1D(x, y, extend=('raise', 0.0))

    """

    cdef bs_spline1d *ptr   # pointer to c struct
    
    def __cinit__(self, x, y, bcs='notaknot', extend='constant'):

        cdef bs_bcs parsed_bcs
        cdef bs_exts parsed_exts
        cdef bs_errorcode code
        
        # convert to double arrays if needed
        cdef double[:] x_ = np.asarray(x, dtype=np.float64)
        cdef double[:] y_ = np.asarray(y, dtype=np.float64)

        # require 1-d arrays
        if (x_.ndim != 1 or y_.ndim != 1):
            raise ValueError("x and y must be 1-d arrays")

        cdef bs_array x_arr = to_bs_array(x_)
        cdef bs_array y_arr = to_bs_array(y_)

        # parse boundary conditions
        if not try_parse_boundary_conditions(bcs, &parsed_bcs):
            raise ValueError("unrecognized boundary condition(s): "+repr(bcs))

        # parse exterior behavior
        if not try_parse_extends(extend, &parsed_exts):
            raise ValueError("unrecognized extend: "+repr(extend))

        code = bs_spline1d_create(x_arr, y_arr, parsed_bcs, parsed_exts,
                                  &self.ptr)
        assert_ok(code)


    #def __init__(self, np.ndarray x, np.ndarray y):
    #    pass


    def __dealloc__(self):
        bs_spline1d_free(self.ptr)

    def __call__(self, double[:] x):
        
        cdef bs_errorcode code
        
        out = np.empty(len(x), dtype=np.float64)
        cdef double[:] outview = out

        code = bs_spline1d_eval(self.ptr, to_bs_array(x),
                                to_bs_array(outview))
        assert_ok(code)

        return out

    def coefficients(self):
        """Return the spline coefficients as a copy"""
        cdef double[:] view = <double[:(self.ptr.n+2)]>(self.ptr.coeffs)
        return np.array(view)  # copy


cdef class USpline1D:
    """
    USpline1D(x, y, bcs='notaknot', extend='constant')

    One dimensional cubic basis spline with uniform grid spacing.
    
    Parameters
    ----------
    x : (float, float)
        2-tuple giving inclusive lower and upper boundaries of abscissa.

    y : `numpy.ndarray` (1-d)
        Ordinate values to interpolate through.

    bcs : str or tuple
        One of:

        - ``'notaknot'``: The "not a knot" condition: the third derivative is
          enforced to be continuous at the first interior knot point.
        - ``'natural'``: set second derivative to zero at boundary.
        - ``'flat'``: set first derivative to zero at boundary.
        - ``('deriv1', value)``: set first derivative to ``value`` at boundary.
        - ``('deriv2', value)``: set second derivative to ``value`` at
          boundary.

        Can also be a 2-tuple of the above values specifying the left and
        right boundary conditions separately.

    extend : str or float or tuple
        Controls what happens when a value outside the spline domain is
        passed.

        - ``'extrapolate'`` Extrapolate past the last knot.
        - ``'constant'`` Return the value at the outermost knot.
        - ``'raise'`` Raise a ``bspline.DomainError``.
        - numeric value: return this value.

        Can also be a 2-tuple of the above values specifying the left and
        right behavior separately.


    Examples
    --------
    
    ::

        USpline1D((xmin, xmax), y)

    See Also
    --------
    Spline1D

    """

    cdef bs_uspline1d *ptr   # pointer to c struct
    
    def __cinit__(self, x, y, bcs='notaknot', extend='constant'):

        cdef bs_bcs parsed_bcs
        cdef bs_exts parsed_exts
        cdef bs_errorcode code
        cdef double[:] x_, y_
        cdef double xmin, xmax

        # convert to double array if needed
        y_ = np.asarray(y, dtype=np.float64)
        if (y_.ndim != 1):
            raise ValueError("y must be 1-d")

        if len(x) <= 2:
            xmin, xmax = x
        else:
            # check if x is an array with uniform spacing
            x_ = np.asarray(x, dtype=np.float64)
            if x_.ndim != 1:
                raise ValueError("x must be 1-d")
            if not is_uniform(x_):
                raise ValueError("x does not have uniform spacing")
            if not (len(x_) == len(y_)):
                raise ValueError("x is array-like but x and y have different "
                                 "lengths")
            xmin, xmax = x_[0], x_[-1]

        # parse boundary conditions
        if not try_parse_boundary_conditions(bcs, &parsed_bcs):
            raise ValueError("unrecognized boundary condition(s): "+repr(bcs))

        # parse exterior behavior
        if not try_parse_extends(extend, &parsed_exts):
            raise ValueError("unrecognized extend: "+repr(extend))

        code = bs_uspline1d_create(bs_range(xmin, xmax), to_bs_array(y_),
                                   parsed_bcs, parsed_exts, &self.ptr)
        assert_ok(code)


    #def __init__(self, np.ndarray x, np.ndarray y):
    #    pass


    def __dealloc__(self):
        bs_uspline1d_free(self.ptr)

    def __call__(self, double[:] x):
        
        cdef bs_errorcode code
        
        out = np.empty(len(x), dtype=np.float64)
        cdef double[:] outview = out

        code = bs_uspline1d_eval(self.ptr, to_bs_array(x),
                                 to_bs_array(outview))
        assert_ok(code)

        return out

    def coefficients(self):
        """Return the spline coefficients as a copy"""
        cdef double[:] view = <double[:(self.ptr.n+2)]>(self.ptr.coeffs)
        return np.array(view)  # copy
