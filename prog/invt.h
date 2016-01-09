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


/* normalize the potential */
static void normalize(double *v, int n)
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


/* compute the root-mean-squared error */
static double geterror(double *v, int n)
{
  double err = 0;
  int i;

  /* subtract the baseline */
  normalize(v, n);

  for ( i = 0; i < n; i++ ) {
    err += v[i] * v[i];
  }
  return sqrt(err / n);
}



#endif /* INVT_H__ */

