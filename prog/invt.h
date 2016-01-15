#ifndef INVT_H__
#define INVT_H__



#include "mtrand.h"
#include "util.h"
#include "invtpar.h"


/* global Metropolis move */
static int mc_metro_g(const double *v, int n, int i)
{
  int j, acc;
  double dv, r;

  j = (int) ( (i + 1 + (n - 1) * rand01()) ) % n;
  dv = v[j] - v[i];
  if ( dv <= 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp(-dv) );
  }

  if ( acc ) {
    i = j;
  }

  return i;
}



/* local Metropolis move */
static int mc_metro_l(const double *v, int n, int i)
{
  int j, acc;
  double dv, r;

  j = (rand01() > 0.5) ? i + 1 : i - 1;
  if ( j < 0 || j >= n ) return i;

  /* periodic boundary condition */
  // j = (j + n) % n;

  dv = v[j] - v[i];
  if ( dv <= 0 ) {
    acc = 1;
  } else {
    r = rand01();
    acc = ( r < exp(-dv) );
  }

  if ( acc ) {
    i = j;
  }

  return i;
}



/* heat-bath move */
static int mc_heatbath(const double *v, double *vac, int n)
{
  int i;
  double vmin = DBL_MAX, r;

  for ( i = 0; i < n; i++ ) {
    if ( v[i] < vmin )
      vmin = v[i];
  }

  vac[0] = 0;
  for ( i = 0; i < n; i++ ) {
    vac[i + 1] = vac[i] + exp(vmin - v[i]);
  }

  r = vac[n] * rand01();

  /* find the bracket that contains r */
  for ( i = n - 1; i >= 0; i-- ) {
    if ( r >= vac[i] )
      break;
  }

  return i;
}



/* multiple-bin update */
static void mbin_update(double *v, int n, int i,
    double a, const double *win, int width)
{
  int j, k;

  v[i] += a * win[0];
  for ( j = 1; j < width; j++ ) {
    k = i - j;
    if ( k < 0 ) {
      k = - k - 1;
    }
    v[k] += a * win[ j ];

    k = i + j;
    if ( k >= n ) {
      k = 2 * n - 1 - k;
    }

    v[k] += a * win[ abs(j) ];
  }
}



/* shift the baseline of the bias potential */
static void shift(double *v, int n)
{
  double s = 0;
  int i;

  for ( i = 0; i < n; i++ ) {
    s += v[i];
  }
  s /= n;
  for ( i = 0; i < n; i++ ) {
    v[i] -= s;
  }
}



/* compute the root-mean-squared error of `v`
 * note: the baseline of `v` will be shifted  */
static double geterror(double *v, int n)
{
  double err = 0;
  int i;

  /* subtract the baseline */
  shift(v, n);

  for ( i = 0; i < n; i++ ) {
    err += v[i] * v[i];
  }
  return sqrt(err / n);
}



/* normalize the error such that the standard deviation is sig */
__inline static double normalize(double *v, int n, double sig)
{
  double err;
  int i;

  err = geterror(v, n);
  for ( i = 0; i < n; i++ ) {
    v[i] *= sig / err;
  }

  return err;
}



/* compute the cosine table */
static double *makecostab(int n)
{
  int i, k;
  double *costab, norm = 1.0, ang;

  xnew(costab, n * n);

  for ( i = 0; i < n; i++ ) {
    costab[i] = norm;
  }

  norm *= sqrt(2);
  ang = M_PI / n;
  for ( k = 1; k < n; k++ ) {
    for ( i = 0; i < n; i++ ) {
      costab[k*n + i] = norm * cos((i + 0.5) * k * ang);
    }
  }

  return costab;
}


/* compute the magnitude of the eigen histogram modes
 * save them to `u` */
static void getcosmodes(double *v, int n, double *u,
    double *costab)
{
  int i, k;
  double s;

  /* this is a simple implementation
   * use FFT for performance */
  for ( k = 0; k < n; k++ ) {
    s = 0;
    for ( i = 0; i < n; i++ ) {
      s += costab[k*n+i] * v[i];
    }
    u[k] = s / n;
  }
}



/* compute the vector from the mode coefficients */
static void fromcosmodes(double *v, int n, double *u,
    double *costab)
{
  int i, k;

  for ( i = 0; i < n; i++ ) {
    v[i] = 0;
  }

  for ( k = 0; k < n; k++ ) {
    for ( i = 0; i < n; i++ ) {
      v[i] += costab[k*n + i] * u[k];
    }
  }
}



#endif /* INVT_H__ */

