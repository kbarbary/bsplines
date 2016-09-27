#include <stdlib.h>
#include <bs.h>

//-----------------------------------------------------------------------------
// search functions
//
// These all return i such that x>= values[i] and x<values[i+1].
// Return -1 if x < values[0].
// Return n-1 if x >= values[n-1].
//-----------------------------------------------------------------------------

// This version assumes we already know that x >= values[start].
// (Use start=-1 for no knowledge.)
static int find_index_from(double *values, int n, double x, int start)
{
  if (start < -1) start = -1;
  if (start > n-1) start = n-1;
  int i = start + 1;
  while (i < n && x >= values[i]) i++;
  return i-1;
}

/*
static int find_index(double *values, int n, double x)
{
  return find_index_from(values, n, x, -1);
}
*/

// find index using binary search
static int find_index_binary(double *values, int n, double x)
{
  int lo = 0;
  int hi = n;
  int mid = n/2;
  
  if (x < values[0]) return -1;
  if (x >= values[n-1]) return n-1;

  while (hi - lo > 1) {
    if (x >= values[mid]) lo = mid;
    else                  hi = mid;
    mid = lo + (hi - lo) / 2;
  }

  return mid;
}


//-----------------------------------------------------------------------------
// Compute the 4 basis functions that are nonzero, assuming t[i] <= x < t[i+1].
// These are: b_{3, i-3}(x), b_{3, i-2}(x), b_{3, i-1}(x), b_{3, i}(x)
//
// This is faster than computing them separately, as some parts of the
// calculation are shared. These can be derived by "manually inlining"
// recursive function calls in the formula for the basis function.
// (See tests for recursive version).
//
// consts[4*i] through consts[4*i+3] stores four constants used in calculation
// (constants are a function of knot spacings).
//
// t indicies from (i-2) to (i+3) are used.
//-----------------------------------------------------------------------------

static void b3nonzeros(double x, int i, double* restrict t,
                       double* restrict consts, double out[4])
{
    double* restrict c = consts + 4*i;

    double dx1 = x - t[i-2];
    double dx2 = x - t[i-1];
    double dx3 = x - t[i];
    double dx4 = t[i+1] - x;
    double dx5 = t[i+2] - x;
    double dx6 = t[i+3] - x;

    double tmp1 = dx4 * dx4 * c[0];
    double tmp2 = dx3 * dx3 * c[1];
    double tmp3 = dx2 * dx4 * c[2] + dx5 * dx3 * c[3];
    
    out[0] = dx4 * tmp1;
    out[1] = dx1 * tmp1 + dx5 * tmp3;
    out[2] = dx6 * tmp2 + dx2 * tmp3;
    out[3] = dx3 * tmp2;
}


// derivatives of previous function
static void db3nonzeros(double x, int i, double* restrict t,
                        double* restrict consts, double out[4])
{
    double* restrict c = consts + 4*i;

    double dx1 = x - t[i-2];
    double dx2 = x - t[i-1];
    double dx3 = x - t[i];
    double dx4 = t[i+1] - x;
    double dx5 = t[i+2] - x;
    double dx6 = t[i+3] - x;
  
    double tmp1 = dx4 * c[0];
    double tmp2 = dx3 * c[1];
    double tmp3 = dx2 * c[2];
    double tmp4 = dx5 * c[3];
  
    out[0] = -3.0 * dx4 * tmp1;
  
    out[1] = ((        dx4 - 2.0 * dx1) * tmp1 +
              (-       dx4 -       dx5) * tmp3 +
              (- 2.0 * dx3 +       dx5) * tmp4 +
              dx5 * dx4 * c[2]);

    out[2] = ((-     dx3 + 2.0 * dx6) * tmp2 +
              (2.0 * dx4 -       dx2) * tmp3 +
              (      dx3 +       dx2) * tmp4
              - dx2 * dx3 * c[3]);

    out[3] = 3.0 * dx3 * tmp2;
}


// second derivatives
static void d2b3nonzeros(double x, int i, double* restrict t,
                        double* restrict consts, double out[4])
{
    double* restrict c = consts + 4*i;

    double dx1 = x - t[i-2];
    double dx2 = x - t[i-1];
    double dx3 = x - t[i];
    double dx4 = t[i+1] - x;
    double dx5 = t[i+2] - x;
    double dx6 = t[i+3] - x;
  
    out[0] = 6.0 * dx4 * c[0];
  
    out[1] = (- 2.0 * dx4 * c[0]
              - 2.0 * (dx4 - dx1) * c[0]
              -       (dx4 - dx2) * c[2]
              +       (-dx5 - dx4) * c[2]
              -       (dx5 - dx2) * c[2]
              - 2.0 * (dx5 - dx3) * c[3]
              - 2.0 * dx5 * c[3]);

    out[2] = (- 2.0 * dx3 * c[1]
              + 2.0 * (dx6 - dx3) * c[1] 
              + 2.0 * (dx4 - dx2) * c[2]
              - 2.0 * dx2 * c[2]
              +       (dx5 - dx3) * c[3]
              -       (dx2 + dx3) * c[3]
              +       (dx5 - dx2) * c[3]);

    out[3] = 6.0 * dx3 * c[1];
}

// third derivatives
static void d3b3nonzeros(int i, double* restrict consts, double out[4])
{
    double* restrict c = consts + 4*i;
  
    out[0] = -6.0 * c[0];
    out[1] =  6.0 * (c[0] + c[2] + c[3]);
    out[2] = -6.0 * (c[1] + c[2] + c[3]);
    out[3] =  6.0 * c[1];
}


//-----------------------------------------------------------------------------
// knots
//-----------------------------------------------------------------------------

// fill spline knots based on x array (includes padding on either
// end of array).
static double* alloc_knots(bs_array x)
{
  int N = x.length;
  double *knots = malloc((N + 5) * sizeof(double));

  // move pointer past initial two-element padding.
  knots += 2;
  
  // copy x into main part of knots
  for (int i=0; i < N; i++)
    knots[i] = x.data[i * x.stride];

  // fill padded area before beginning
  knots[-2] = knots[0] - 2.0 * (knots[1] - knots[0]);
  knots[-1] = knots[0] - 1.0 * (knots[1] - knots[0]);

  // fill padded area after end.
  knots[N]   = knots[N-1] + 1.0 * (knots[N-1] - knots[N-2]);
  knots[N+1] = knots[N-1] + 2.0 * (knots[N-1] - knots[N-2]);
  knots[N+2] = knots[N-1] + 3.0 * (knots[N-1] - knots[N-2]);

  return knots;
}


static void free_knots(double *knots) {
  free(knots - 2);
}


// constants used when evaluating a spline.
// constaants + 4*i is a pointer to the four constants used
// when evaluating the spline in the range knots[i] <= x < knots[i+1].
static double* alloc_constants(double *knots, int n) {
    double *constants = malloc(4 * n * sizeof(double));

    for (int i=0; i<n; i++) {
        constants[4*i+0] = 1.0 / ((knots[i+1] - knots[i-2]) *
                                  (knots[i+1] - knots[i-1]) *
                                  (knots[i+1] - knots[i  ]));
        
        constants[4*i+1] = 1.0 / ((knots[i+3] - knots[i  ]) *
                                  (knots[i+2] - knots[i  ]) *
                                  (knots[i+1] - knots[i  ]));

        constants[4*i+2] = 1.0 / ((knots[i+2] - knots[i-1]) *
                                  (knots[i+1] - knots[i-1]) *
                                  (knots[i+1] - knots[i  ]));

        constants[4*i+3] = 1.0 / ((knots[i+2] - knots[i-1]) *
                                  (knots[i+2] - knots[i  ]) *
                                  (knots[i+1] - knots[i  ]));
    }

    return constants;
}

//-----------------------------------------------------------------------------
// not-a-knot boundary condition preprocessing
//
// We want to solve a matrix like this:
//
// | x x x x x          |
// | x x x              |
// |   x x x            |
// |     x x x          |
// |        ...         |
// |            x x x   |
// |              x x x |
// |          x x x x x |
//
// So we'll eliminate the trailing two elements in the first row
// and the leading two elements in the last row. Then we can feed it to
// the standard solve().
//
// strategy is basically:
// subtract (row 3) * (A[4,0] / A[4,3]) from row 0
// subtract (row 2) * (A[3,0] / A[3,2]) from row 0


//-----------------------------------------------------------------------------
// solve()
//
// Solve A * x = b for x. The solution is stored in b.
//
// A is an almost tridiagonal n x n matrix with this form:
//
// | x x x              |
// | x x x              |
// |   x x x            |
// |        ...         |
// |            x x x   |
// |              x x x |
// |              x x x |
//
// Only the non-zero elements on each row are stored, so A has 3n elements.
// Rows are contiguous.
//-----------------------------------------------------------------------------

static void solve(double* restrict A, double* restrict b, int n)
{
  // divide 0th row by upper left element
  double t = A[0];
  b[0] /= t;
  A[2] /= t;  
  A[1] /= t;
  // A[0] = 1.0; but not used again.

  // subtract (first element of row 1) x (row 0) from row 1
  // to eliminate first element of row 1.
  t = A[3*1+0];
  b[1]     -= t * b[0];
  A[3*1+2] -= t * A[2];
  A[3*1+1] -= t * A[1];
  // A[3*1+0] = 0.0; but not used again.
  
  // divide row 1 by first nonzero element, to set it to 1.
  t = A[3*1+1];
  b[1]     /= t;
  A[3*1+2] /= t;
  // A[3*1+1] = 1.0; but not used again.

  for (int i=2; i<n-1; i++) {

    // subtract (first element of new row) * (previous row) from new row
    // to eliminate first element.
    t = A[3*i+0];
    b[i]        -= t * b[i-1];
    // A[3*i+2] -= t * 0.0  // no-op b/c previous row is zero.
    A[3*i+1]    -= t * A[3*(i-1)+2];
    // A[3*i+0] -= t * 1.0;  // (previous row is 1.0) but not used again.

    // divide new row by first non-zero element
    t = A[3*i+1];
    b[i]     /= t;
    A[3*i+2] /= t;
    // A[3*i+1] = 1.0;
  }

  // last row is different:
  // subtract first element of last row * 3rd to last row from last row
  b[n-1]          -= A[3*(n-1)+0] * b[n-3];
  // A[3*(n-1)+2] -= A[3*(n-1)+0] * 0.0; // no-op
  A[3*(n-1)+1]    -= A[3*(n-1)+0] * A[3*(n-3)+2];
  A[3*(n-1)+0]    = 0.0;
  
  // subtract first non-zero element * previous row from last row
  b[n-1]       -= A[3*(n-1)+1] * b[n-2];
  A[3*(n-1)+2] -= A[3*(n-1)+1] * A[3*(n-2)+2];
  A[3*(n-1)+1] = 0.0;

  // divide row by 1st non-zero element
  b[n-1]       /= A[3*(n-1)+2];
  A[3*(n-1)+2] =  1.0;

  // back substitute
  for (int i=n-2; i>0; i--) {
    b[i] -= b[i+1] * A[3*i+2];
  }

  // first row is different
  b[0] -= b[1] * A[1] + b[2] * A[2];
}

//-----------------------------------------------------------------------------
// bs_array
//-----------------------------------------------------------------------------

static int is_monotonic(bs_array x)
{
  int ok = 1;
  for (int i=1; i<x.length; i++) {
    ok &= (x.data[i*x.stride] >= x.data[(i-1)*x.stride]);
  }
  return ok;
}


//-----------------------------------------------------------------------------
// spline1d
//-----------------------------------------------------------------------------
bs_errorcode bs_spline1d_create(bs_array x, bs_array y, bs_bcs bcs,
                                bs_exts exts, bs_spline1d **out)
{
  // Initialize output pointer, in case we error out.
  *out = NULL;
  
  // checks
  if (x.length != y.length) return BS_LENGTHMISMATCH;
  if (!is_monotonic(x)) return BS_NOTMONOTONIC;

  bs_spline1d* spline = malloc(sizeof(bs_spline1d));
  
  int N = x.length;
  int M = N + 2;
  
  spline->knots = alloc_knots(x);
  spline->n = N;
  spline->exts = exts;
  spline->consts = alloc_constants(spline->knots, N);

  // process "constant" extends
  if (spline->exts.left.type == BS_CONSTANT) {
      spline->exts.left.type = BS_VALUE;
      spline->exts.left.value = y.data[0];
  }
  if (spline->exts.right.type == BS_CONSTANT) {
      spline->exts.right.type = BS_VALUE;
      spline->exts.right.value = y.data[(N-1)*y.stride];
  }
  
  // sparse matrix (last element not used, but we write to it)
  double *A = malloc((3*M + 1) * sizeof(double));

  // The first row is the left boundary condition
  if (bcs.left.type == BS_DERIV1) {
      db3nonzeros(spline->knots[0], 0, spline->knots, spline->consts, A);
  }
  else {  // type must be BS_DERIV2
      d2b3nonzeros(spline->knots[0], 0, spline->knots, spline->consts, A);
  }

  // At a knot i, the spline has three non-zero components:
  // b_{i-3}, b_{i-2}, b{i-1}. (b_i is zero at knot i).
  for (int i=0; i<N; i++) {
      b3nonzeros(spline->knots[i], i, spline->knots, spline->consts,
                 A + 3*(i+1));
  }

  // derivatives at final point
  if (bcs.right.type == BS_DERIV1) {
      db3nonzeros(spline->knots[N-1], N-1, spline->knots, spline->consts,
                  A + 3*(M-1));
  }
  else {  // type must be BS_DERIV2
      d2b3nonzeros(spline->knots[N-1], N-1, spline->knots, spline->consts,
                   A + 3*(M-1));
  }

  // right hand side:
  double *b = malloc(M * sizeof(double));
  b[0] = bcs.left.value;
  for (int i=0; i<N; i++) {
      b[i+1] = y.data[i * y.stride];
  }
  b[M-1] = bcs.right.value;
  
  // Solve
  solve(A, b, M);
  free(A);
  spline->coeffs = b;

  *out = spline;
  return BS_OK;
}


void bs_spline1d_free(bs_spline1d* spline)
{
  if (spline != NULL) {
    free_knots(spline->knots);
    free(spline->consts);
    free(spline->coeffs);
    free(spline);
  }
}


bs_errorcode bs_spline1d_eval(bs_spline1d *spline, bs_array x, bs_array out)
{
  // ensure that x is increasing
  if (!is_monotonic(x)) return BS_NOTMONOTONIC;
  
  // for first index, it could be anywhere, so use binary search
  int i = find_index_binary(spline->knots, spline->n, x.data[0]);

  double xval;
  double b3vals[4];
  for (int j=0; j<x.length; j++) {
    xval = x.data[j*x.stride];
    i = find_index_from(spline->knots, spline->n, xval, i);

    // index outside left boundary
    if (i == -1) {
      switch (spline->exts.left.type) {
      case BS_EXTRAPOLATE:
        i = 0;
        break;
      case BS_VALUE:
        out.data[j * out.stride] = spline->exts.left.value;
        continue;
      case BS_RAISE:
        return BS_DOMAINERROR;
      }
    }

    // index outside right boundary
    else if (i == spline->n - 1) {
      switch (spline->exts.right.type) {
      case BS_EXTRAPOLATE:
        i = spline->n - 2;
        break;
      case BS_VALUE:
        out.data[j * out.stride] = spline->exts.right.value;
        continue;
      case BS_RAISE:
        return BS_DOMAINERROR;
      }
    }

    // if we get this far, we're either extrapolating or xval is in range.
    b3nonzeros(xval, i, spline->knots, spline->consts, b3vals);
    out.data[j*out.stride] = (spline->coeffs[i]   * b3vals[0] +
                              spline->coeffs[i+1] * b3vals[1] +
                              spline->coeffs[i+2] * b3vals[2] +
                              spline->coeffs[i+3] * b3vals[3]);
  }

  return BS_OK;
}


//-----------------------------------------------------------------------------
// unit basis versions of b3nonzeros and friends
// knot locations in this basis are [0, 1, ..., N-1]
// For 0 <= x < 1 return b_{-3}(x), b_{-2}(x), b_{-1}(x), b_0(x)
// For i <= x < i+1 subtract i from x first. to get b_{i-3}(x), b_{i-2}(x), ...
// (works because all the basis functions are the same with a shift,
//  so b_{-3}(x) = b_{i-3}(i+x).
//-----------------------------------------------------------------------------

static const ONESIXTH = 0.1666666666666666666;

static void b3unonzeros(double x, double out[4])
{

    double dx1 = x + 2.0;
    double dx2 = x + 1.0;
    double dx4 = 1.0 - x;
    double dx5 = 2.0 - x;
    double dx6 = 3.0 - x;

    double tmp1 = ONESIXTH * dx4 * dx4;
    double tmp2 = ONESIXTH * x * x;
    double tmp3 = ONESIXTH * (dx2 * dx4 + dx5 * x);
    
    out[0] = dx4 * tmp1;
    out[1] = dx1 * tmp1 + dx5 * tmp3;
    out[2] = dx6 * tmp2 + dx2 * tmp3;
    out[3] = x   * tmp2;
}


// derivatives of previous function
static void db3unonzeros(double x, double out[4])
{

    double dx1 = x + 2.0;
    double dx2 = x + 1.0;
    double dx4 = 1.0 - x;
    double dx5 = 2.0 - x;
    double dx6 = 3.0 - x;
  
    double tmp1 = ONESIXTH * dx4;
    double tmp2 = ONESIXTH * x;
    double tmp3 = ONESIXTH * dx2;
    double tmp4 = ONESIXTH * dx5;
  
    out[0] = -3.0 * dx4 * tmp1;
  
    out[1] = ((        dx4 - 2.0 * dx1) * tmp1 +
              (-       dx4 -       dx5) * tmp3 +
              (- 2.0 * x +       dx5) * tmp4 +
              ONESIXTH * dx5 * dx4);

    out[2] = ((-     x + 2.0 * dx6) * tmp2 +
              (2.0 * dx4 -       dx2) * tmp3 +
              (      x +       dx2) * tmp4
              - ONESIXTH * dx2 * x);

    out[3] = 3.0 * x * tmp2;
}


// second derivatives
static void d2b3unonzeros(double x, double out[4])
{
    double dx1 = x + 2.0;
    double dx2 = x + 1.0;
    double dx4 = 1.0 - x;
    double dx5 = 2.0 - x;
    double dx6 = 3.0 - x;
  
    out[0] = ONESIXTH * 6.0 * dx4;
  
    out[1] = ONESIXTH * (- 2.0 * dx4
                         - 2.0 * (dx4 - dx1)
                         -       (dx4 - dx2)
                         +       (-dx5 - dx4)
                         -       (dx5 - dx2)
                         - 2.0 * (dx5 - x)
                         - 2.0 * dx5);

    out[2] = ONESIXTH * (- 2.0 * x
                         + 2.0 * (dx6 - x)
                         + 2.0 * (dx4 - dx2)
                         - 2.0 * dx2
                         +       (dx5 - x)
                         -       (dx2 + x)
                         +       (dx5 - dx2));

    out[3] = ONESIXTH * 6.0 * x;
}

// third derivatives
static void d3b3nonzeros(double out[4])
{
  
    out[0] = -1.0;
    out[1] =  3.0;
    out[2] = -3.0;
    out[3] =  1.0;
}
