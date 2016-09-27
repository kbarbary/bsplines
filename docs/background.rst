=======================
Mathematical Background
=======================

Solving for coefficients
------------------------

At any given point, there are at most four non-zero basis
functions. For a point :math:`x` between knot points :math:`x_i` and
:math:`x_{i+1}`, these are :math:`b_{i-3}`, :math:`b_{i-2}`,
:math:`b_{i-1}` and :math:`b_i`, so the value of the spline is given
by

.. math::

   y(x) = c_{i-3} b_{i-3}(x) + c_{i-2} b_{i-2}(x) + c_{i-1} b_{i-1}(x) + c_{i} b_{i}(x)

To solve for the coefficients, we demand that the spline passes
through the data points: that :math:`y(x_i) = y_i`, giving us

.. math::

   c_{i-3} b_{i-3}(x_i) + c_{i-2} b_{i-2}(x_i) + c_{i-1} b_{i-1}(x_i) = y_i

Note that the last term drops out because :math:`b_i(x_i)` is zero.
However, this isn't enough to uniquely solve for the coefficients: We
need some additional constraints, which is where boundary conditions
come in.

**Derivative boundary conditions**

One option is to specify that a **derivative** of the spline obtain some given value
at the endpoints. Taking the derivative of the above equation gives

.. math::

   c_{i-3} b'_{i-3}(x_i) + c_{i-2} b'_{i-2}(x_i) + c_{i-1} b'_{i-1}(x_i) = y'_i


Specifying the derivative at the endpoints gives us the system of equations:

.. math::

   \left( \begin{array}{cccccc}
   b'_{-3}(x_0) & b'_{-2}(x_0) & b'_{-1}(x_0) &            \\
   b_{-3}(x_0)  & b_{-2}(x_0)  & b_{-1}(x_0)  &            \\
                & b_{-2}(x_1)  & b_{-1}(x_1)  & b_{0}(x_1) \\
                \\
                &              &              & \ddots     \\
                \\
                &              & b_{N-5}(x_{N-2})  & b_{N-4}(x_{N-2})  & b_{N-3}(x_{N-2})  &   \\
                &              &                   & b_{N-4}(x_{N-1})  & b_{N-3}(x_{N-1})  & b_{N-2}(x_{N-1}) \\
                &              &                   & b'_{N-4}(x_{N-1})  & b'_{N-3}(x_{N-1})  & b'_{N-2}(x_{N-1})
   \end{array} \right)
   \left(\begin{array}{c}
   c_{-3} \\
   c_{-2} \\
   c_{-1} \\
   \\
   \vdots \\
   \\
   c_{N-4} \\
   c_{N-3} \\
   c_{N-2}
   \end{array} \right)
   =
   \left(\begin{array}{c}
   y'_0 \\
   y_0 \\
   y_1 \\
   \\
   \vdots \\
   \\
   y_{N-2} \\
   y_{N-1} \\
   y'_{N-1}
   \end{array} \right)

which we can now solve for the coefficients. For second derivative
conditions, replace :math:`b'(x)` by :math:`b''(x)` and :math:`y'` by
:math:`y''`. Specifically, setting the second derivative to zero at
the endpoints is known as "natural" boundary conditions.

**"Not-a-knot" boundary condition**

Another option for a boundary condition is the "not-a-knot" condition:
that the *third* derivative is continuous at the first interior knot
point:

.. math::

   c_{-3} b'''_{-3}(x_<) + c_{-2} b'''_{-2}(x_<) + c_{-1} b'''_{-1}(x_<) + c_{0} b'''_{0}(x_<) = \\
   c_{-2} b'''_{-2}(x_1) + c_{-1} b'''_{-1}(x_1) + c_{0} b'''_{0}(x_1) + c_{1} b'''_{1}(x_1)

where :math:`x_<` is a value just smaller than :math:`x_1`. The third
derivatives of the basis functions are constant between knots
(:math:`x_i \le x < x_{i+1}`), so :math:`x_<` can be replaced by
:math:`x_0`. Subtract the right hand side from both sides:

.. math::

   c_{-3} b'''_{-3}(x_0) + c_{-2} (b'''_{-2}(x_0) - b'''_{-2}(x_1))  + c_{-1} (b'''_{-1}(x_0) - b'''_{-1}(x_1))  + {} \\
   + c_{0} (b'''_{0}(x_0) - b'''_{0}(x_1)) - c_{1} b'''_{1}(x_1) = 0

The first row in the matrix in the original equation then becomes:

.. math::

   \left( \begin{array}{cccccc}
   b'''_{-3}(x_0) & b'''_{-2}(x_0) - b'''_{-2}(x_1) & b'''_{-1}(x_0) - b'''_{-1}(x_1) & b'''_{0}(x_0) - b'''_{0}(x_1) & -b'''_{1}(x_1) & \ldots
   \end{array} \right)

Similarly, the last row becomes:

.. math::
   
   \left( \begin{array}{ccccc}
   \ldots & b'''_{N-6}(x_{N-3}) & b'''_{N-5}(x_{N-3}) - & b'''_{N-4}(x_{N-3}) -  & b'''_{N-3}(x_{N-3}) - & -b'''_{N-2}(x_{N-2}) \\
          &                     & \qquad  b'''_{N-5}(x_{N-2}) & \qquad   b'''_{N-4}(x_{N-2}) & \qquad  b'''_{N-3}(x_{N-2}) &
   \end{array} \right)


and the first and last elements on the right hand side are replaced with 0.

We can use the third and forth rows to reduce the first row to three
non-zero elements to get this in a form like the original equation,
and similarly for the last row.
