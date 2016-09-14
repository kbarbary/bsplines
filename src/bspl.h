typedef struct {
  double *data;
  int length;
  int stride;
} bspl_array;

double bspl_b3(double x, int i, double *t);
double bspl_db3(double x, int i, double *t);
double bspl_ddb3(double x, int i, double *t);

//-----------------------------------------------------------------------------
// Boundary conditions
//-----------------------------------------------------------------------------

typedef struct {
  int deriv;
  double value;
} bspl_bc;

typedef struct {
  bspl_bc left;
  bspl_bc right;
} bspl_bcs;

//-----------------------------------------------------------------------------
// 1-d splines
//-----------------------------------------------------------------------------

typedef struct {
  double *knots;
  double *coeffs;
  int n;
} bspl_spline1d;

bspl_spline1d* bspl_create_spline1d(bspl_array x, bspl_array y, bspl_bcs bcs);
double bspl_eval_spline1d(bspl_spline1d *spline, double x);
void bspl_free_spline1d(bspl_spline1d *spline);
