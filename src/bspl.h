typedef struct {
  double *data;
  int length;
  int stride;
} bs_array;

double bs_b3(double x, int i, double *t);
double bs_db3(double x, int i, double *t);
double bs_ddb3(double x, int i, double *t);

//-----------------------------------------------------------------------------
// Boundary conditions
//-----------------------------------------------------------------------------

typedef struct {
  int deriv;
  double value;
} bs_bc;

typedef struct {
  bs_bc left;
  bs_bc right;
} bs_bcs;

//-----------------------------------------------------------------------------
// 1-d splines
//-----------------------------------------------------------------------------

typedef struct {
  double *knots;
  double *coeffs;
  int n;
} bs_spline1d;

bs_spline1d* bs_create_spline1d(bs_array x, bs_array y, bs_bcs bcs);
double bs_eval_spline1d(bs_spline1d *spline, double x);
int bs_evalvec_spline1d(bs_spline1d *spline, bs_array x, bs_array out);
void bs_free_spline1d(bs_spline1d *spline);
