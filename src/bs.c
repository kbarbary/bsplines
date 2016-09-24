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
// Explicitly written out version of i-th cubic B-spline
// basis function, b_{3, i}. This can be derived by "manually inlining"
// recursive function calls in the formula for the basis function.
// (See tests for recursive version).
//
// This assumes that t[i] through t[i+4] (inclusive) are valid indicies.
//-----------------------------------------------------------------------------

// b_{3,i}(x) for t[i] <= x < t[i+1]
double b3_0(double x, int i, double *t) {
  return
    (x - t[i]) * (x - t[i]) * (x - t[i]) /
    ((t[i+3] - t[i]) * (t[i+2] - t[i]) * (t[i+1] - t[i]));
}

// b_{3, i}(x) for t[i+1] <= x < t[i+2]
double b3_1(double x, int i, double *t) {
  return
    (x - t[i]) / (t[i+3] - t[i]) *
    ((x - t[i]) * (t[i+2] - x) / ((t[i+2] - t[i]) * (t[i+2] - t[i+1]))
     +
     (t[i+3] - x) * (x - t[i+1]) / ((t[i+3] - t[i+1]) * (t[i+2] - t[i+1])))
    +
    (t[i+4] - x) * (x - t[i+1]) * (x - t[i+1]) /
    ((t[i+4] - t[i+1]) * (t[i+3] - t[i+1]) * (t[i+2] - t[i+1]));
}

// b_{3, i}(x) for t[i+2] <= x < t[i+3]
double b3_2(double x, int i, double *t) {
  return
    (x - t[i]) * (t[i+3] - x) * (t[i+3] - x) /
    ((t[i+3] - t[i]) * (t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
    +
    (t[i+4] - x) / (t[i+4] - t[i+1]) *
    ((x - t[i+1]) * (t[i+3] - x) / ((t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
     +
     (t[i+4] - x) * (x - t[i+2]) / ((t[i+4] - t[i+2]) * (t[i+3] - t[i+2])));
}

// b_{3, i}(x) for t[i+3] <= x < t[i+4]
double b3_3(double x, int i, double *t) {
  return
    (t[i+4] - x) * (t[i+4] - x) * (t[i+4] - x) /
    ((t[i+4] - t[i+1]) * (t[i+4] - t[i+2]) * (t[i+4] - t[i+3]));
}

// b_{3,i}(x)
double bs_b3(double x, int i, double *t) {
  if (x < t[i]) return 0.0;
  else if (x < t[i+1]) return b3_0(x, i, t);
  else if (x < t[i+2]) return b3_1(x, i, t);
  else if (x < t[i+3]) return b3_2(x, i, t);
  else if (x < t[i+4]) return b3_3(x, i, t);
  else return 0.0;
}


//-----------------------------------------------------------------------------
// version of above with precomputed differences
//
// dtinv[3i+j] stores 1./(t[i+1+j] - t[i])

// b_{3,i}(x) for t[i] <= x < t[i+1]
double b3_0_pre(double x, int i, double *t, double *dtinv) {
  return
    (x - t[i]) * (x - t[i]) * (x - t[i]) *
    dtinv[3*i+2] * dtinv[3*i+1] * dtinv[3*i];
}

// b_{3, i}(x) for t[i+1] <= x < t[i+2]
double b3_1_pre(double x, int i, double *t, double *dtinv) {
  return
    (x - t[i]) * dtinv[3*i+2] *
    ((x - t[i]) * (t[i+2] - x) * dtinv[3*i+1] * dtinv[3*(i+1)]
     +
     (t[i+3] - x) * (x - t[i+1]) * dtinv[3*(i+1)+1] * dtinv[3*(i+1)])
    +
    (t[i+4] - x) * (x - t[i+1]) * (x - t[i+1]) *
    dtinv[3*(i+1)+2] * dtinv[3*(i+1)+1] * dtinv[3*(i+1)];
}

// b_{3, i}(x) for t[i+2] <= x < t[i+3]
double b3_2_pre(double x, int i, double *t, double *dtinv) {
  return
    (x - t[i]) * (t[i+3] - x) * (t[i+3] - x) *
    dtinv[3*i+2] * dtinv[3*(i+1)+1] * dtinv[3*(i+2)]
    +
    (t[i+4] - x) * dtinv[3*(i+1)+2] *
    ((x - t[i+1]) * (t[i+3] - x) * dtinv[3*(i+1)+1] * dtinv[3*(i+2)]
     +
     (t[i+4] - x) * (x - t[i+2]) * dtinv[3*(i+2)+1] * dtinv[3*(i+2)]);
}

// b_{3, i}(x) for t[i+3] <= x < t[i+4]
double b3_3_pre(double x, int i, double *t, double *dtinv) {
  return
    (t[i+4] - x) * (t[i+4] - x) * (t[i+4] - x) *
    dtinv[3*(i+1)+2] * dtinv[3*(i+2)+1] * dtinv[3*(i+3)];
}

// b_{3, i-3}(x), b_{3, i-2}(x), b_{3, i-1}(x), b_{3, i}(x)
// assuming t[i] <= x < t[i+1]
void b3vec_pre(double x, int i, double *t, double *dtinv, double out[4])
{
  double dx0 = x - t[i-2];
  double dx1 = x - t[i-1];
  double dx2 = x - t[i];
  double dx3 = t[i+1] - x;
  double dx4 = t[i+2] - x;
  double dx5 = t[i+3] - x;
  
  double c1 = dtinv[3*(i-2)+2] * dtinv[3*(i-1)+1] * dtinv[3*i];
  double c2 = dtinv[3* i   +2] * dtinv[3* i   +1] * dtinv[3*i];
  double c3 = dtinv[3*(i-1)+2] * dtinv[3*(i-1)+1] * dtinv[3*i];
  double c4 = dtinv[3*(i-1)+2] * dtinv[3* i   +1] * dtinv[3*i];
  
  double tmp1 = dx3 * dx3 * c1;
  double tmp2 = dx2 * dx2 * c2;
  double tmp3 = dx1 * dx3 * c3 + dx4 * dx2 * c4;
  
  out[0] = dx3 * tmp1;
  out[1] = dx0 * tmp1 + dx4 * tmp3;
  out[2] = dx5 * tmp2 + dx1 * tmp3;
  out[3] = dx2 * tmp2;
}

//-----------------------------------------------------------------------------
// first derivatives of above functions
//-----------------------------------------------------------------------------

// d/dx b_{3,i}(x) for t[i] <= x < t[i+1]
double db3_0(double x, int i, double *t) {
  return
    3.0 * (x - t[i]) * (x - t[i]) /
    ((t[i+3] - t[i]) * (t[i+2] - t[i]) * (t[i+1] - t[i]));
}

// d/dx b_{3, i}(x) for t[i+1] <= x < t[i+2]
double db3_1(double x, int i, double *t) {
  return
    1.0 / (t[i+3] - t[i]) *
    ((x - t[i]) * (t[i+2] - x) / ((t[i+2] - t[i]) * (t[i+2] - t[i+1]))
     +
     (t[i+3] - x) * (x - t[i+1]) / ((t[i+3] - t[i+1]) * (t[i+2] - t[i+1])))
    +
    (x - t[i]) / (t[i+3] - t[i]) *
    ((t[i+2] + t[i] - 2.0 * x) / ((t[i+2] - t[i]) * (t[i+2] - t[i+1]))
     +
     (t[i+3] + t[i+1] - 2.0 * x) / ((t[i+3] - t[i+1]) * (t[i+2] - t[i+1])))    
    +
    (2.0 * (t[i+4] - x) * (x - t[i+1]) - (x - t[i+1]) * (x - t[i+1])) /
    ((t[i+4] - t[i+1]) * (t[i+3] - t[i+1]) * (t[i+2] - t[i+1]));
}

// d/dx b_{3, i}(x) for t[i+2] <= x < t[i+3]
double db3_2(double x, int i, double *t) {
  return
    ((t[i+3] - x) * (t[i+3] - x) - 2.0 * (x - t[i]) * (t[i+3] - x)) /
    ((t[i+3] - t[i]) * (t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
    -
    1.0 / (t[i+4] - t[i+1]) *
    ((x - t[i+1]) * (t[i+3] - x) / ((t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
     +
     (t[i+4] - x) * (x - t[i+2]) / ((t[i+4] - t[i+2]) * (t[i+3] - t[i+2])))
    +
    (t[i+4] - x) / (t[i+4] - t[i+1]) *
    ((t[i+3] + t[i+1] - 2.0 * x) / ((t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
     +
     (t[i+4] + t[i+2] - 2.0 * x) / ((t[i+4] - t[i+2]) * (t[i+3] - t[i+2])));
}

// d/dx b_{3, i}(x) for t[i+3] <= x < t[i+4]
double db3_3(double x, int i, double *t) {
  return
    -3.0 * (t[i+4] - x) * (t[i+4] - x) /
    ((t[i+4] - t[i+1]) * (t[i+4] - t[i+2]) * (t[i+4] - t[i+3]));
}

double bs_db3(double x, int i, double *t) {
  if (x < t[i]) return 0.0;
  else if (x < t[i+1]) return db3_0(x, i, t);
  else if (x < t[i+2]) return db3_1(x, i, t);
  else if (x < t[i+3]) return db3_2(x, i, t);
  else if (x < t[i+4]) return db3_3(x, i, t);
  else return 0.0;
}

//-----------------------------------------------------------------------------
// second derivatives
//-----------------------------------------------------------------------------

// d2/dx2 b_{3,i}(x) for t[i] <= x < t[i+1]
double ddb3_0(double x, int i, double *t) {
  return
    6.0 * (x - t[i]) /
    ((t[i+3] - t[i]) * (t[i+2] - t[i]) * (t[i+1] - t[i]));
}

// d2/dx2 b_{3, i}(x) for t[i+1] <= x < t[i+2]
double ddb3_1(double x, int i, double *t) {
  return
    1.0 / (t[i+3] - t[i]) *
    ((t[i+2] + t[i] - 2.0 * x) / ((t[i+2] - t[i]) * (t[i+2] - t[i+1]))
     +
     (t[i+3] + t[i+1] - 2.0 * x) / ((t[i+3] - t[i+1]) * (t[i+2] - t[i+1])))

    +
    
    1.0 / (t[i+3] - t[i]) *
    ((t[i+2] + t[i] - 2.0 * x) / ((t[i+2] - t[i]) * (t[i+2] - t[i+1]))
     +
     (t[i+3] + t[i+1] - 2.0 * x) / ((t[i+3] - t[i+1]) * (t[i+2] - t[i+1])))    

    +
    
    (x - t[i]) / (t[i+3] - t[i]) *
    (-2.0 / ((t[i+2] - t[i]) * (t[i+2] - t[i+1]))
     -2.0 / ((t[i+3] - t[i+1]) * (t[i+2] - t[i+1])))    

    +

    (2.0 * (t[i+4] + t[i+1] - 2.0 * x) - 2.0 * (x - t[i+1])) /
    ((t[i+4] - t[i+1]) * (t[i+3] - t[i+1]) * (t[i+2] - t[i+1]));
}

// d2/dx2 b_{3, i}(x) for t[i+2] <= x < t[i+3]
double ddb3_2(double x, int i, double *t) {
  return
    (-2.0 * (t[i+3] - x) - 2.0 * (t[i+3] + t[i] - 2.0 * x)) /
    ((t[i+3] - t[i]) * (t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))

    +

    -1.0 / (t[i+4] - t[i+1]) *
    ((t[i+3] + t[i+1] - 2.0 * x) / ((t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
     +
     (t[i+4] + t[i+2] - 2.0 * x) / ((t[i+4] - t[i+2]) * (t[i+3] - t[i+2])))

    +
     
    -1.0 / (t[i+4] - t[i+1]) *
    ((t[i+3] + t[i+1] - 2.0 * x) / ((t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
     +
     (t[i+4] + t[i+2] - 2.0 * x) / ((t[i+4] - t[i+2]) * (t[i+3] - t[i+2])))

    +
     
    (t[i+4] - x) / (t[i+4] - t[i+1]) *
    (-2.0 / ((t[i+3] - t[i+1]) * (t[i+3] - t[i+2]))
     -2.0 / ((t[i+4] - t[i+2]) * (t[i+3] - t[i+2])));
}

// d2/dx2 b_{3, i}(x) for t[i+3] <= x < t[i+4]
double ddb3_3(double x, int i, double *t) {
  return
    6.0 * (t[i+4] - x) /
    ((t[i+4] - t[i+1]) * (t[i+4] - t[i+2]) * (t[i+4] - t[i+3]));
}

double bs_ddb3(double x, int i, double *t) {
  if (x < t[i]) return 0.0;
  else if (x < t[i+1]) return ddb3_0(x, i, t);
  else if (x < t[i+2]) return ddb3_1(x, i, t);
  else if (x < t[i+3]) return ddb3_2(x, i, t);
  else if (x < t[i+4]) return ddb3_3(x, i, t);
  else return 0.0;
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

  // move pointer pointer past initial two-element padding.
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
  knots -= 2;
  free(knots);
}

// dtinv[3i+j] stores 1./(t[i+1+j] - t[i])
static double* alloc_dtinv(double *knots, int n) {
  double *dtinv = malloc(3 * (n+2) * sizeof(double));
  dtinv += 6;
  
  for (int i=-2; i<n; i++) {
    dtinv[3*i+0] = 1.0 / (knots[i+1] - knots[i]);
    dtinv[3*i+1] = 1.0 / (knots[i+2] - knots[i]);
    dtinv[3*i+2] = 1.0 / (knots[i+3] - knots[i]);
  }

  return dtinv;
}

static void free_dtinv(double *dtinv) {
  dtinv -= 6;
  free(dtinv);
}
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
  b[0] -= b[1] * A[1] - b[2] * A[2];
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
  spline->dtinv = alloc_dtinv(spline->knots, N);

  // sparse matrix
  double *A = malloc(3 * M * sizeof(double));

  // The first row is the constraint on a derivative of the spline at x[0]
  if (bcs.left.type == BS_DERIV1) {
    A[0] = db3_3(spline->knots[0], -3, spline->knots);
    A[1] = db3_2(spline->knots[0], -2, spline->knots);
    A[2] = db3_1(spline->knots[0], -1, spline->knots);
  }
  else {  // type must be BS_DERIV2
    A[0] = ddb3_3(spline->knots[0], -3, spline->knots);
    A[1] = ddb3_2(spline->knots[0], -2, spline->knots);
    A[2] = ddb3_1(spline->knots[0], -1, spline->knots);
  }

  // At a knot i, the spline has three non-zero components:
  // b_{i-3}, b_{i-2}, b{i-1}. (b_i is zero at knot i).
  for (int i=0; i<N; i++) {
    A[3*(i+1) + 0] = b3_3(spline->knots[i], i-3, spline->knots);
    A[3*(i+1) + 1] = b3_2(spline->knots[i], i-2, spline->knots);
    A[3*(i+1) + 2] = b3_1(spline->knots[i], i-1, spline->knots);
  }

  // derivatives at final point
  if (bcs.right.type == BS_DERIV1) {
    A[3*(M-1) + 0] = db3_3(spline->knots[N-1], (N-1)-3, spline->knots);
    A[3*(M-1) + 1] = db3_2(spline->knots[N-1], (N-1)-2, spline->knots);
    A[3*(M-1) + 2] = db3_1(spline->knots[N-1], (N-1)-1, spline->knots);
  }
  else {  // type must be BS_DERIV2
    A[3*(M-1) + 0] = ddb3_3(spline->knots[N-1], (N-1)-3, spline->knots);
    A[3*(M-1) + 1] = ddb3_2(spline->knots[N-1], (N-1)-2, spline->knots);
    A[3*(M-1) + 2] = ddb3_1(spline->knots[N-1], (N-1)-1, spline->knots);
  }

  // right hand side:
  double *b = malloc(M * sizeof(double));
  b[0] = bcs.left.value;
  for (int i=0; i<N; i++)
    b[i+1] = y.data[i * y.stride];
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
    free_dtinv(spline->dtinv);
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
      case BS_CONSTANT:
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
      case BS_CONSTANT:
        out.data[j * out.stride] = spline->exts.right.value;
        continue;
      case BS_RAISE:
        return BS_DOMAINERROR;
      }
    }

    // if we get this far, we're either extrapolating or xval is in range.
    b3vec_pre(xval, i, spline->knots, spline->dtinv, b3vals);
    out.data[j*out.stride] = (spline->coeffs[i]   * b3vals[0] +
                              spline->coeffs[i+1] * b3vals[1] +
                              spline->coeffs[i+2] * b3vals[2] +
                              spline->coeffs[i+3] * b3vals[3]);
  }

  return BS_OK;
}
