#ifndef BS_H
#define BS_H

//-----------------------------------------------------------------------------
// Error codes
//-----------------------------------------------------------------------------

typedef enum {
  BS_OK           = 0,
  BS_OUTOFMEMORY  = 1,
  BS_DOMAINERROR  = 2,
  BS_NOTMONOTONIC = 3,
  BS_LENGTHMISMATCH = 4,
} bs_errorcode;


//-----------------------------------------------------------------------------
// Input data types
//-----------------------------------------------------------------------------

typedef struct {
    double *data;
    int length;
    int stride;
} bs_array;

typedef struct {
    double min; // inclusive
    double max; // inclusive
} bs_range;

//-----------------------------------------------------------------------------
// Boundary conditions
//-----------------------------------------------------------------------------

typedef enum {BS_DERIV1, BS_DERIV2, BS_NOTAKNOT} bs_bctype;

typedef struct {
  bs_bctype type;
  double value;
} bs_bc;

typedef struct {
  bs_bc left;
  bs_bc right;
} bs_bcs;

//-----------------------------------------------------------------------------
// out-of-domain behavior ("extension")
//-----------------------------------------------------------------------------

typedef enum {BS_EXTRAPOLATE, BS_CONSTANT, BS_VALUE, BS_RAISE} bs_exttype;

typedef struct {
  bs_exttype type;
  double value;
} bs_ext;

typedef struct {
  bs_ext left;
  bs_ext right;
} bs_exts;

//-----------------------------------------------------------------------------
// 1-d splines
//-----------------------------------------------------------------------------

typedef struct {
  double *knots;
  double *consts;
  double *coeffs;
  int n;
  bs_exts exts;
} bs_spline1d;

bs_errorcode bs_spline1d_create(bs_array x, bs_array y, bs_bcs bcs,
                                bs_exts exts, bs_spline1d **out);
bs_errorcode bs_spline1d_eval(bs_spline1d *spline, bs_array x, bs_array out);
void         bs_spline1d_free(bs_spline1d *spline);

#endif

typedef struct {
    double xmin;
    double xmax;
    double *coeffs;
    int n;
    bs_exts exts;
} bs_uspline1d;

bs_errorcode bs_uspline1d_create(bs_range x, bs_array y,
                                 bs_bcs bcs, bs_exts exts, bs_spline1d **out);
