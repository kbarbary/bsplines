#include <stdlib.h>
#include <math.h>
#include <bs.h>

#include <stdio.h> //debug

//debug
void print_a_and_b(double first[5], double last[5],
                   double *A, double  *b, int M)
{
    printf("\nfirst: [ %f  %f  %f  %f  %f ]\n",
           first[0], first[1], first[2], first[3], first[4]);

    for (int i=0; i<M; i++)
        printf("row %d : | %f  %f  %f |    | %f |\n",
               i, A[3*i+0], A[3*i+1], A[3*i+2], b[i]);

    printf("last: [ %f  %f  %f  %f  %f ]\n",
           last[0], last[1], last[2], last[3], last[4]);

}


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
                       double* restrict consts, double out[restrict 4])
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
  int N = x.size;
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
// solve_simple()
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
// Rows are contiguous in memory: e.g., A[0] through A[2]
// stores the second row (first row with three elements).
//
//-----------------------------------------------------------------------------

static void solve_simple(double* restrict A, double* restrict b, int n)
{
  // divide first row by upper left element
  double t = A[0];
  b[0] /= t;
  A[2] /= t;  
  A[1] /= t;
  A[0] = 1.0; // but not used again.

  // subtract (first element of row 1) x (row 0) from row 1
  // to eliminate first element of row 1.
  t = A[3*1+0];
  b[1]     -= t * b[0];
  A[3*1+2] -= t * A[2];
  A[3*1+1] -= t * A[1];
  A[3*1+0] = 0.0; // but not used again.
  
  // divide row 1 by first nonzero element, to set it to 1.
  t = A[3*1+1];
  b[1]     /= t;
  A[3*1+2] /= t;
  A[3*1+1] = 1.0; // but not used again.

  for (int i=2; i<n-1; i++) {

    // subtract (first element of new row) * (previous row) from new row
    // to eliminate first element.
    t = A[3*i+0];
    b[i]        -= t * b[i-1];
    // A[3*i+2] -= t * 0.0  // no-op b/c previous row is zero.
    A[3*i+1]    -= t * A[3*(i-1)+2];
    A[3*i+0] = 0.0;  // (previous row is 1.0) but not used again.

    // divide new row by first non-zero element
    t = A[3*i+1];
    b[i]     /= t;
    A[3*i+2] /= t;
    A[3*i+1] = 1.0;
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
// solve()
//
// Solve A * x = b for x. The solution is stored in b.
//
// A is a matrix like this:
//
// | x x x x x          |
// | x x x              |
// |   x x x            |
// |        ...         |
// |          x x x     |
// |            x x x   |
// |              x x x |
// |          x x x x x |
//
// A is stored compactly in 3*n elements, with row i corresponding to A[3*i],
// with the exception of the first and last rows which are passed in
// separately because they are too large to be stored this way.
//
// Note that the first 3 and last 3 elements of A are initially empty as
// these row values are stored in `first` and `last`. In fact the last 3
// elements of A are not used at all.
//
//-----------------------------------------------------------------------------

static void solve(double first[restrict 5], double last[restrict 5],
                   double* restrict A, double* restrict b, int n)
{
    // rows 1, 2, 3: divide by first non-zero
    //
    // x x x x x | y       x x x x x | y
    // x x x     | y       1 x x     | y
    //   x x x   | y  -->    1 x x   | y
    //     x x x | y           1 x x | y
    
    for (int i=1; i<4; i++) {
        b[i]     /= A[3*i];
        A[3*i+2] /= A[3*i];
        A[3*i+1] /= A[3*i];
        A[3*i]   = 1.0;
    }

    // eliminate first two elements of first row and divide by first non-zero.
    //
    // x x x x x | y       0 0 1 x x | y
    // 1 x x     | y       1 x x     | y
    //   1 x x   | y  -->    1 x x   | y
    //     1 x x | y           1 x x | y
    b[0]     -= first[0] * b[1];
    first[2] -= first[0] * A[3*1+2];
    first[1] -= first[0] * A[3*1+1];
    first[0] = 0.0;

    b[0]     -= first[1] * b[2];
    first[3] -= first[1] * A[3*2+2];
    first[2] -= first[1] * A[3*2+1];
    first[1] = 0.0;

    b[0]     /= first[2];
    first[4] /= first[2];
    first[3] /= first[2];
    first[2] = 1.0;

    // reduce row 3
    //
    // 0 0 1 x x | y       0 0 1 x x | y
    // 1 x x     | y       1 x x     | y
    //   1 x x   | y  -->    1 x x   | y
    //     1 x x | y           0 1 x | y   
    b[3]     -= A[3*3+0] * b[0];
    A[3*3+2] -= A[3*3+0] * first[4];
    A[3*3+1] -= A[3*3+0] * first[3];
    A[3*3+0] = 0.0;

    b[3]     /= A[3*3+1];
    A[3*3+2] /= A[3*3+1];
    A[3*3+1] = 1.0;

    // permute first three rows:
    // 0 0 1 x x | y       1 x x     | y
    // 1 x x     | y         1 x x   | y
    //   1 x x   | y  -->      1 x x | y
    //     0 1 x | y           0 1 x | y
    double tmp = b[0];
    b[0] = b[1];
    A[3*0+0] = A[3*1+0];
    A[3*0+1] = A[3*1+1];
    A[3*0+2] = A[3*1+2];

    b[1] = b[2];
    A[3*1+0] = A[3*2+0];
    A[3*1+1] = A[3*2+1];
    A[3*1+2] = A[3*2+2];

    b[2] = tmp;
    A[3*2+0] = first[2];
    A[3*2+1] = first[3];
    A[3*2+2] = first[4];
    
    // reduce rest of the middle rows
    for (int i=4; i<n-1; i++) {
        b[i]     -= A[3*i+0] * b[i-1];
        A[3*i+1] -= A[3*i+0] * A[3*(i-1)+2];
        A[3*i+0] = 0.0;

        b[i]     /= A[3*i+1];
        A[3*i+2] /= A[3*i+1];
        A[3*i+1] = 1.0;
    }

    // we now have, e.g.,
    // 1 x x         | y
    //   1 x x       | y
    //     1 x x     | y  (n-5)
    //     0 1 x     | y  (n-4)
    //       0 1 x   | y  (n-3)
    //         0 1 x | y  (n-2)
    //     x x x x x | y  (n-1)

    // eliminate first element of last row using the (n-5)th row.
    b[n-1] -= last[0] * b[n-5];
    if (n-5 < 3) {
        last[2] -= last[0] * A[3*(n-5)+2];
        last[1] -= last[0] * A[3*(n-5)+1];
    }
    else {
        last[1] -= last[0] * A[3*(n-5)+2];
    }
    last[0] = 0.0;

    // eliminate second element of last row using the (n-4)th row.
    b[n-1] -= last[1] * b[n-4];
    if (n-4 < 3) {
        last[3] -= last[1] * A[3*(n-4)+2];
        last[2] -= last[1] * A[3*(n-4)+1];
    }
    else {
        last[2] -= last[1] * A[3*(n-4)+2];
    }
    last[1] = 0.0;

    // eliminate third element of last row using the (n-3)rd row.
    b[n-1] -= last[2] * b[n-3];
    if (n-3 < 3) {
        last[4] -= last[2] * A[3*(n-3)+2];
        last[3] -= last[2] * A[3*(n-3)+1];
    }
    else {
        last[3] -= last[2] * A[3*(n-3)+2];
    }
    last[2] = 0.0;

    // eliminate forth element
    b[n-1] -= last[3] * b[n-2];
    last[4] -= last[3] * A[3*(n-2)+2];
    last[3] = 0.0;

    // normalize last row
    b[n-1] /= last[4];
    last[4] = 1.0;

    // back-substitute
    for (int i=n-2; i>=3; i--) {
        b[i] -= b[i+1] * A[3*i+2];
    }

    // we now have:
    // 1 x x           | y
    //   1 x x         | y
    //     1 x x       | y
    //       1         | y
    //         1       | y
    //          ...
    //
    // eliminate the remaining elements.
    b[2] -= b[3] * A[3*2+1] + b[4] * A[3*2+2];
    b[1] -= b[2] * A[3*1+1] + b[3] * A[3*1+2];
    b[0] -= b[1] * A[3*0+1] + b[2] * A[3*0+2];
}


//-----------------------------------------------------------------------------
// bs_array
//-----------------------------------------------------------------------------

static int is_monotonic(bs_array x)
{
  int ok = 1;
  for (int i=1; i<x.size; i++) {
    ok &= (x.data[i*x.stride] >= x.data[(i-1)*x.stride]);
  }
  return ok;
}


//-----------------------------------------------------------------------------
// spline1d
//-----------------------------------------------------------------------------


static void notaknot_row(double *consts, int i, double row[5])
{
    double buf[4];

    d3b3nonzeros(i-1, consts, row);
    d3b3nonzeros(i, consts, buf);
    row[4] = 0.0;
    for (int i=0; i<4; i++) {
        row[i+1] -= buf[i];
    }
}

// TODO: function to just fill A, so we can do it just once in each
// dimension for spline2d? Problem is:
// - we'd have to also pass around first, last
// - b would be different each time so we'd hae to fill it each time
//   separately, including for boundary conditions.

// Find spline coefficients along one dimension.
// knots and consts are as belong to a spline.
// coeffs is size values.size + 2 with stride cstride.
// A is a buffer with size 3*(values.size + 2)
static void find_1d_coefficients(double* restrict knots,
                                 double* restrict consts,
                                 bs_array values, bs_bcs bcs,
                                 double* restrict A,
                                 double* restrict coeffs)
{
    double first[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double last[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int N = values.size;
    int M = N+2;
    
    // fill rows 1 through M-1 with values of b_{i-3}, b_{i-2}, b{i-1}
    // at knot i.
    for (int i=0; i<N; i++) {
        b3nonzeros(knots[i], i, knots, consts, A + 3*(i+1));
        coeffs[i+1] = values.data[i * values.stride];
    }

    // Left boundary condition
    switch (bcs.left.type) {
    case BS_DERIV1:
        db3nonzeros(knots[0], 0, knots, consts, first);
        coeffs[0] = bcs.left.value;
        break;
    case BS_DERIV2:
        d2b3nonzeros(knots[0], 0, knots, consts, first);
        coeffs[0] = bcs.left.value;
        break;
    case BS_NOTAKNOT:
        notaknot_row(consts, 1, first);
        coeffs[0] = 0.0;
    }

    // Right boundary condition
    switch (bcs.right.type) {
    case BS_DERIV1:
        db3nonzeros(knots[N-1], N-1, knots, consts,
                    last+2);
        coeffs[M-1] = bcs.right.value;
        break;
    case BS_DERIV2:
        d2b3nonzeros(knots[N-1], N-1, knots, consts,
                     last+2);
        coeffs[M-1] = bcs.right.value;
        break;
    case BS_NOTAKNOT:
        notaknot_row(consts, N-2, last);
        coeffs[M-1] = 0.0;
    }

    solve(first, last, A, coeffs, M);
}


bs_errorcode bs_spline1d_create(bs_array x, bs_array y, bs_bcs bcs,
                                bs_exts exts, bs_spline1d **out)
{
    *out = NULL;  // In case of error, ensure that output pointer is NULL.
  
  // checks
  if (x.size != y.size) return BS_SIZEMISMATCH;
  if (!is_monotonic(x)) return BS_NOTMONOTONIC;

  bs_spline1d* spline = malloc(sizeof(bs_spline1d));
  
  int N = x.size;
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

  double *A = malloc(3 * M * sizeof(double));  // sparse row representation
  double *coeffs = malloc(M * sizeof(double)); 

  find_1d_coefficients(spline->knots, spline->consts, y, bcs, A, coeffs);
  free(A);

  spline->coeffs = coeffs;
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
  for (int j=0; j<x.size; j++) {
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
// spline2d
//-----------------------------------------------------------------------------

bs_errorcode bs_spline2d_create(bs_array x, bs_array y, bs_array2d z,
                                bs_bcs xbcs, bs_bcs ybcs, bs_exts xexts,
                                bs_exts yexts, bs_spline2d **out)
{
    *out = NULL;  // In case of error, ensure that output pointer is NULL.

    if ((x.size != z.sizes[0]) || (y.size != z.sizes[1]))
        return BS_SIZEMISMATCH;
    if (!is_monotonic(x) || !is_monotonic(y))
        return BS_NOTMONOTONIC;

    bs_spline2d* spline = malloc(sizeof(bs_spline2d));
  
    int nx = x.size;
    int mx = nx + 2;

    int ny = y.size;
    int my = ny + 2;

    spline->xknots = alloc_knots(x);
    spline->xconsts = alloc_constants(spline->xknots, nx);
    spline->nx = nx;
    spline->xexts = xexts;

    spline->yknots = alloc_knots(y);
    spline->yconsts = alloc_constants(spline->yknots, ny);
    spline->ny = ny;
    spline->yexts = yexts;

    double* coeffs = malloc(mx * my * sizeof(double));
    
    // find coefficients along y (fast axis)
    double *A = malloc(3 * my * sizeof(double));
    for (int i=0; i<nx; i++) {
        bs_array zslice = {z.data + z.strides[1]*i, z.sizes[0], z.strides[0]};
        find_1d_coefficients(spline->yknots, spline->yconsts, zslice,
                             ybcs, A, coeffs+(i*my));
    }
    free(A);

    // find coefficients along x (slow axis);
    A = malloc(3 * mx * sizeof(double));
    double *buf = malloc(mx * sizeof(double));
    for (int i=0; i<my; i++) {
        // for this slice in constant y, the target values are the
        // `nx` coefficients we just found. They are strided in
        // `coeffs` by `my`.
        bs_array coeffs_slice = {coeffs + i, nx, my};
        find_1d_coefficients(spline->yknots, spline->yconsts, coeffs_slice,
                             ybcs, A, buf);
        
        // the results in `buf` are contiguous in x, but we need to
        // copy them back into the coefficients array strided.
        for (int j=0; j<mx; j++) coeffs[i+my*j] = buf[j];
    }
    free(A);
    free(buf);

    spline->coeffs = coeffs;
    *out = spline;

    return BS_OK;
}

void bs_spline2d_free(bs_spline2d* spline)
{
  if (spline != NULL) {
    free_knots(spline->xknots);
    free(spline->xconsts);
    free_knots(spline->yknots);
    free(spline->yconsts);
    free(spline->coeffs);
    free(spline);
  }
}

bs_errorcode bs_spline2d_eval(bs_spline2d *spline, bs_array x, bs_array y,
                              bs_array2d out)
{
    // ensure inputs are increasing.
    if (!is_monotonic(x)) return BS_NOTMONOTONIC;
    if (!is_monotonic(y)) return BS_NOTMONOTONIC;
  
    // for first index, it could be anywhere, so use binary search
    int i = find_index_binary(spline->xknots, spline->nx, x.data[0]);
    int j0 = find_index_binary(spline->yknots, spline->ny, y.data[0]);

    int my = spline->ny + 2; // for indexing coeffs.
    double xb3vals[4];
    double yb3vals[4];
    for (int k=0; k<x.size; k++) {
        double xval = x.data[k*x.stride];
        i = find_index_from(spline->xknots, spline->nx, xval, i);

        // index outside left boundary
        if (i == -1) {
            switch (spline->xexts.left.type) {
            case BS_EXTRAPOLATE:
                i = 0;
                break;
            case BS_VALUE:
                for (int l=0; l<y.size; l++) {
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->xexts.left.value;
                }
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // index outside right boundary
        else if (i == spline->nx - 1) {
            switch (spline->xexts.right.type) {
            case BS_EXTRAPOLATE:
                i = spline->nx - 2;
                break;
            case BS_VALUE:
                for (int l=0; l<y.size; l++) {
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->xexts.right.value;
                }
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // get basis function values for x coordinate
        b3nonzeros(xval, i, spline->xknots, spline->xconsts, xb3vals);

        // x value is in range (or extrapolating); loop over y values:
        int j = j0;
        for (int l=0; l<y.size; l++) {
            double yval = y.data[l*y.stride];
            j = find_index_from(spline->yknots, spline->ny, yval, j);

            // index outside left boundary
            if (j == -1) {
                switch (spline->yexts.left.type) {
                case BS_EXTRAPOLATE:
                    j = 0;
                    break;
                case BS_VALUE:
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->yexts.left.value;
                    continue;
                case BS_RAISE:
                    return BS_DOMAINERROR;
                }
            }

            // index outside right boundary
            else if (j == spline->ny - 1) {
                switch (spline->yexts.right.type) {
                case BS_EXTRAPOLATE:
                    j = spline->ny - 2;
                    break;
                case BS_VALUE:
                    out.data[k * out.strides[0] + l * out.strides[1]] =
                        spline->yexts.right.value;
                    continue;
                case BS_RAISE:
                    return BS_DOMAINERROR;
                }
            }

            // get basis function values for y coordinate
            b3nonzeros(yval, j, spline->yknots, spline->yconsts, yb3vals);

            out.data[k * out.strides[0] + l * out.strides[1]] =
                (spline->coeffs[(i  )*my+j]   * xb3vals[0] * yb3vals[0] +
                 spline->coeffs[(i  )*my+j+1] * xb3vals[0] * yb3vals[1] +
                 spline->coeffs[(i  )*my+j+2] * xb3vals[0] * yb3vals[2] +
                 spline->coeffs[(i  )*my+j+3] * xb3vals[0] * yb3vals[3] +

                 spline->coeffs[(i+1)*my+j]   * xb3vals[1] * yb3vals[0] +
                 spline->coeffs[(i+1)*my+j+1] * xb3vals[1] * yb3vals[1] +
                 spline->coeffs[(i+1)*my+j+2] * xb3vals[1] * yb3vals[2] +
                 spline->coeffs[(i+1)*my+j+3] * xb3vals[1] * yb3vals[3] +
                 
                 spline->coeffs[(i+2)*my+j]   * xb3vals[2] * yb3vals[0] +
                 spline->coeffs[(i+2)*my+j+1] * xb3vals[2] * yb3vals[1] +
                 spline->coeffs[(i+2)*my+j+2] * xb3vals[2] * yb3vals[2] +
                 spline->coeffs[(i+2)*my+j+3] * xb3vals[2] * yb3vals[3] +

                 spline->coeffs[(i+3)*my+j]   * xb3vals[3] * yb3vals[0] +
                 spline->coeffs[(i+3)*my+j+1] * xb3vals[3] * yb3vals[1] +
                 spline->coeffs[(i+3)*my+j+2] * xb3vals[3] * yb3vals[2] +
                 spline->coeffs[(i+3)*my+j+3] * xb3vals[3] * yb3vals[3]);
        }
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

static const double ONESIXTH  = 0.1666666666666666666;

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
static void d3b3unonzeros(double out[4])
{
  
    out[0] = -1.0;
    out[1] =  3.0;
    out[2] = -3.0;
    out[3] =  1.0;
}


bs_errorcode bs_uspline1d_create(bs_range x, bs_array y, bs_bcs bcs,
                                 bs_exts exts, bs_uspline1d **out)
{
    bs_uspline1d* spline = malloc(sizeof(bs_uspline1d));

    int N = y.size;
    int M = N + 2;
    double didx = (N - 1) / (x.max - x.min); // equal to  1 / (step size)

    // values, first and second derivatives of b3_{i-3}, b3_{i-2}, b3_{i-1}
    // at knot i.
    const double b3vals[3] = {1.0/6.0, 2.0/3.0, 1.0/6.0};
    const double db3vals[3] = {-0.5 * didx, 0.0, 0.5 * didx};
    const double d2b3vals[3] = {didx * didx, -2.0 * didx * didx, didx * didx};
    const double notaknot_row[5] = {-1.0, 4.0, -6.0, 4.0, -1.0};
    
    spline->x = x;
    spline->didx = didx;

    spline->n = N;
    spline->exts = exts;

    // process "constant" extends
    if (spline->exts.left.type == BS_CONSTANT) {
        spline->exts.left.type = BS_VALUE;
        spline->exts.left.value = y.data[0];
    }
    if (spline->exts.right.type == BS_CONSTANT) {
        spline->exts.right.type = BS_VALUE;
        spline->exts.right.value = y.data[(N-1)*y.stride];
    }

    double first[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double last[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double *A = malloc(3 * M * sizeof(double));  // sparse row representation
    double *b = malloc(M * sizeof(double));

  // fill rows 1 through M-1 with values of b_{i-3}, b_{i-2}, b{i-1} at knot i.
  for (int i=0; i<N; i++) {
      for (int j=0; j<3; j++) A[3*(i+1)+j] = b3vals[j];
      b[i+1] = y.data[i * y.stride];
  }

  // Left boundary condition
  switch (bcs.left.type) {
  case BS_DERIV1:
      for (int i=0; i<3; i++) first[i] = db3vals[i];
      b[0] = bcs.left.value;
      break;
  case BS_DERIV2:
      for (int i=0; i<3; i++) first[i] = d2b3vals[i];
      b[0] = bcs.left.value;
      break;
  case BS_NOTAKNOT:
      for (int i=0; i<5; i++) first[i] = notaknot_row[i];
      b[0] = 0.0;
  }

  // Right boundary condition
  switch (bcs.right.type) {
  case BS_DERIV1:
      for (int i=0; i<3; i++) last[i+2] = db3vals[i];
      b[M-1] = bcs.right.value;
      break;
  case BS_DERIV2:
      for (int i=0; i<3; i++) last[i+2] = d2b3vals[i];
      b[M-1] = bcs.right.value;
      break;
  case BS_NOTAKNOT:
      for (int i=0; i<5; i++) last[i] = notaknot_row[i];
      b[M-1] = 0.0;
  }


  solve(first, last, A, b, M);
  free(A);
  spline->coeffs = b;

  *out = spline;
  return BS_OK;
}


void bs_uspline1d_free(bs_uspline1d* spline)
{
    if (spline != NULL) {
        free(spline->coeffs);
        free(spline);
    }
}

bs_errorcode bs_uspline1d_eval(bs_uspline1d *spline, bs_array x, bs_array out)
{
    int i;
    double xval;
    double xfloor;
    double b3vals[4];
    for (int j=0; j<x.size; j++) {
        // translate x onto unit basis
        xval = (x.data[j*x.stride] - spline->x.min) * spline->didx;
        xfloor = floor(xval);
        i = (int)xfloor;

        // index outside left boundary
        if (i < 0) {
            switch (spline->exts.left.type) {
            case BS_EXTRAPOLATE:
                i = 0;
                xfloor = 0.0;
                break;
            case BS_VALUE:
                out.data[j * out.stride] = spline->exts.left.value;
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // index outside right boundary
        else if (i >= spline->n-1) {
            switch (spline->exts.right.type) {
            case BS_EXTRAPOLATE:
                i = spline->n - 2;
                xfloor = i;
                break;
            case BS_VALUE:
                out.data[j * out.stride] = spline->exts.right.value;
                continue;
            case BS_RAISE:
                return BS_DOMAINERROR;
            }
        }

        // if we get this far, we're either extrapolating or xval is in range.
        b3unonzeros(xval - xfloor, b3vals);
        out.data[j*out.stride] = (spline->coeffs[i]   * b3vals[0] +
                                  spline->coeffs[i+1] * b3vals[1] +
                                  spline->coeffs[i+2] * b3vals[2] +
                                  spline->coeffs[i+3] * b3vals[3]);
    }

    return BS_OK;
}
