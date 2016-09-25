// test for internal C functions
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// we include the source code here rather than linking
// so we can test 'static' functions.
#include "bs.c"

#define ASSERT(condition) {if (!(condition)) return 1;};


// slow, recursive implementation of basis functions to compare to:
double b(int k, int i, double *t, double x)
{
    if (k <= 0)
        return (x >= t[i] && x < t[i+1])? 1.0: 0.0;
  
    return ((x - t[i]) / (t[i+k] - t[i]) * b(k-1, i, t, x) +
            (t[i+k+1] - x) / (t[i+k+1] - t[i+1]) * b(k-1, i+1, t, x));
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

int close(double x, double y)
{
  double RTOL = 1e-8;
  return (abs(x - y) <= RTOL * x) ? 1 : 0;

}

int allclose(double *x, double *y, int n)
{
  for (int i=0; i<n; i++) {
    if (close(x[i], y[i]) == 0) return 0;
  }
  return 1;
}

//-----------------------------------------------------------------------------
// tests

int test_bs_b3()
{
    int N = 61;
    double knots[5] = {0., 1.3, 2., 3., 4.};
    double *x = linspace(-1., 5., N);
    double *y = malloc(N * sizeof(double));
    double *ytrue = malloc(N * sizeof(double));
  
    for (int j=0; j<N; j++) {
        y[j] = bs_b3(x[j], 0, knots);
        ytrue[j] = b(3, 0, knots, x[j]);
    }

    int status = !allclose(y, ytrue, N);
  
    free(x);
    free(y);
    free(ytrue);

    return status;
}


int test_bfuncs()
{
    int status = 0;
    double knots[9] = {0., 0.6, 1.0, 1.1, 1.5, 1.9, 2.3, 3., 4.};
    double values[3];
    double true_values[3];

    for (int i=3; i<6; i++) {
        values[0] = b3_3(knots[i], i-3, knots);
        values[1] = b3_2(knots[i], i-2, knots);
        values[2] = b3_1(knots[i], i-1, knots);

        true_values[0] = b(3, i-3, knots, knots[i]);
        true_values[1] = b(3, i-2, knots, knots[i]);
        true_values[2] = b(3, i-1, knots, knots[i]);
    
        if (!allclose(true_values, values, 3)) status = 1;
    }

    return status;
}


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

//-----------------------------------------------------------------------------
// test harness

#define RUNTEST(name)                                            \
{                                                                \
    status = name();                                             \
    printf(#name "() ");                                         \
    if (status) printf("\033[31;1mfailed\033[0m\n");             \
    else        printf("\n");                                    \
    total += 1;                                                  \
    failed += (status > 0);                                      \
}


int main()
{
    int status = 0;
    int total = 0;
    int failed = 0;

    RUNTEST(test_bs_b3);
    RUNTEST(test_bfuncs);
    RUNTEST(test_find_index_binary);
    RUNTEST(test_is_monotonic);

    printf("\n");
    if (failed > 0) printf("\033[31;1m");
    else            printf("\033[32;1m");
    printf("%d tests: %d passed, %d failed\n", total, total - failed, failed);
    printf("\033[0m");

    return (failed > 0);
}
