// test for internal C functions
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// we include the source code here rather than linking
// so we can test 'static' functions.
#include "bs.c"

//-----------------------------------------------------------------------------
// test helpers

char assert_message[80] = "\0";
#define ASSERT(condition) {                                     \
        if (!(condition)) {                                     \
            strcpy(assert_message, "failed: " #condition);      \
            return 1;                                           \
        }                                                       \
}

#define ASSERT_CLOSE(ex1, ex2) {                                        \
        if (!(close(ex1, ex2, 1e-14, 1e-14))) {                         \
            strcpy(assert_message, "not close: " #ex1 " and " #ex2);    \
            return 1;                                                   \
        }                                                               \
}

#define ASSERT_ALLCLOSE(ex1, ex2, len) {                                \
        if (!(allclose(ex1, ex2, len, 1e-14, 1e-14))) {                 \
            strcpy(assert_message, "not close: " #ex1 " and " #ex2);    \
            return 1;                                                   \
        }                                                               \
}

// useful for debugging
void print_comparison(double *x, double *y, int n)
{
    for (int i=0; i<n; i++) {
        printf("%.17f   %.17f\n", x[i], y[i]);
    }
}


//-----------------------------------------------------------------------------
// slow, recursive implementation of basis functions to compare to


double b(int k, double x, int i, double *t)
{
    if (k == 0) return (x >= t[i] && x < t[i+1])? 1.0: 0.0;
  
    return ((x - t[i]) / (t[i+k] - t[i]) * b(k-1, x, i, t) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * b(k-1, x, i+1, t));
}


double db(int k, double x, int i, double *t)
{
    if (k == 0) return 0.0;

    return ((x - t[i]) / (t[i+k] - t[i]) * db(k-1, x, i, t) +
            1.0 / (t[i+k] - t[i]) * b(k-1, x, i, t) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * db(k-1, x, i+1, t) +
            -1.0 / (t[i+k+1] - t[i+1]) * b(k-1, x, i+1, t));
}


double d2b(int k, double x, int i, double *t)
{
    if (k == 0) return 0.0;

    return ((x - t[i]) / (t[i+k] - t[i]) * d2b(k-1, x, i, t) +
            2.0 / (t[i+k] - t[i]) * db(k-1, x, i, t) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * d2b(k-1, x, i+1, t) +
            -2.0 / (t[i+k+1] - t[i+1]) * db(k-1, x, i+1, t));
}


double d3b(int k, double x, int i, double *t)
{
    if (k == 0) return 0.0;

    return ((x - t[i])     / (t[i+k] - t[i])     * d3b(k-1, x, i, t) +
            3.0            / (t[i+k] - t[i])     * d2b(k-1, x, i, t) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * d3b(k-1, x, i+1, t) +
            -3.0           / (t[i+k+1] - t[i+1]) * d2b(k-1, x, i+1, t));
}

void slow_b3nonzeros(double x, int i, double *t, double out[4])
{
    out[0] = b(3, x, i-3, t);
    out[1] = b(3, x, i-2, t);
    out[2] = b(3, x, i-1, t);
    out[3] = b(3, x, i  , t);
}

void slow_db3nonzeros(double x, int i, double *t, double out[4])
{
    out[0] = db(3, x, i-3, t);
    out[1] = db(3, x, i-2, t);
    out[2] = db(3, x, i-1, t);
    out[3] = db(3, x, i  , t);
}

void slow_d2b3nonzeros(double x, int i, double *t, double out[4])
{
    out[0] = d2b(3, x, i-3, t);
    out[1] = d2b(3, x, i-2, t);
    out[2] = d2b(3, x, i-1, t);
    out[3] = d2b(3, x, i  , t);
}

void slow_d3b3nonzeros(double x, int i, double *t, double out[4])
{
    out[0] = d3b(3, x, i-3, t);
    out[1] = d3b(3, x, i-2, t);
    out[2] = d3b(3, x, i-1, t);
    out[3] = d3b(3, x, i  , t);
}


//-----------------------------------------------------------------------------
// test helpers

double* linspace(double start, double stop, int n)
{
  double *x = malloc(n * sizeof(double));
  for (int i=0; i<n; i++)
    x[i] = start + (start - stop) * (i / (n-1));

  return x;
}

int close(double x, double y, double rtol, double atol)
{
    return (fabs(x - y) <= rtol * fabs(x) + atol) ? 1 : 0;
}

int allclose(double *x, double *y, int n, double rtol, double atol)
{
    for (int i=0; i<n; i++) {
        if (close(x[i], y[i], rtol, atol) == 0) return 0;
    }
    return 1;
}


//-----------------------------------------------------------------------------
// tests


int test_find_index_binary()
{
    double knots[9] = {0., 0.6, 1.0, 1.1, 1.5, 1.9, 2.3, 3., 4.};

    ASSERT(find_index_binary(knots, 9, 0.) == 0);
    ASSERT(find_index_binary(knots, 9, 4.) == 8);
    ASSERT(find_index_binary(knots, 9, -0.0000001) == -1);
    ASSERT(find_index_binary(knots, 9, 0.9999) == 1);
    ASSERT(find_index_binary(knots, 9, 0.5999) == 0);

    return 0;
}


int test_is_monotonic()
{
    double data[5] = {0., 1., 2., 3., 4.};
    bs_array x = {data, 5, 1};

    ASSERT(is_monotonic(x) == 1);

    x.data[1] = 0.;
    ASSERT(is_monotonic(x) == 1);

    x.data[3] = 2.;
    ASSERT(is_monotonic(x) == 1);

    x.data[1] = -0.00001;
    ASSERT(is_monotonic(x) == 0);

    x.data[3] = 1.99999;
    ASSERT(is_monotonic(x) == 0);

    return 0;
}


int test_b3nonzeros()
{
    double allknots[15] = {-2.1, -2.0, -1.2, 0.0, 0.2, 1.0,
                            1.1,  1.5, 1.9, 2.3, 3.0,
                            4.0,  6.0, 7.0, 8.0, 8.1};

    // knots is length 10 with padding of 3 and 3 before and after.
    double *knots = allknots + 3;

    double *consts = alloc_constants(knots, 10);

    double bvals[4];
    double bvals_ref[4];

    int N = 100;
    double *x = linspace(0.0, 6.0, N);
    for (int i=0; i<N; i++) {

        int j = find_index_binary(knots, 10, x[i]);

        // values
        b3nonzeros(x[i], j, knots, consts, bvals);
        slow_b3nonzeros(x[i], j, knots, bvals_ref);
        ASSERT_ALLCLOSE(bvals, bvals_ref, 4);

        // derivatives
        db3nonzeros(x[i], j, knots, consts, bvals);
        slow_db3nonzeros(x[i], j, knots, bvals_ref);
        ASSERT_ALLCLOSE(bvals, bvals_ref, 4);

        // second derivatives
        d2b3nonzeros(x[i], j, knots, consts, bvals);
        slow_d2b3nonzeros(x[i], j, knots, bvals_ref);
        ASSERT_ALLCLOSE(bvals, bvals_ref, 4);

        // third derivatives
        d3b3nonzeros(j, consts, bvals);
        slow_d3b3nonzeros(x[i], j, knots, bvals_ref);
        ASSERT_ALLCLOSE(bvals, bvals_ref, 4);

    }

    free(consts);  // leaks memory if a test fails, but whatevs.
    free(x);
    return 0;
}
           
           
//-----------------------------------------------------------------------------
// test harness

#define RUNTEST(name) {                                          \
        status = name();                                         \
        printf(#name "() ");                                     \
        if (status) printf("\033[31;1mfailed\033[0m\n");         \
        else        printf("\n");                                \
        total += 1;                                              \
        failed += (status > 0);                                  \
    }


int main()
{
    int status = 0;
    int total = 0;
    int failed = 0;

    RUNTEST(test_find_index_binary);
    RUNTEST(test_is_monotonic);
    RUNTEST(test_b3nonzeros);

    printf("\n");
    if (failed > 0) printf("\033[31;1m");
    else            printf("\033[32;1m");
    printf("%d tests: %d passed, %d failed\n", total, total - failed, failed);
    printf("\033[0m");
    printf("%s\n", assert_message);
    return (failed > 0);
}
