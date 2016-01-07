#ifndef INVT_H__
#define INVT_H__



#include "mtrand.h"
#include "util.h"
#include "invtpar.h"


/* global Metropolis move */
static int mc_metro_g(double *v, int n, int i)
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
static int mc_metro_l(double *v, int n, int i)
{
  int j, acc;
  double dv, r;

  j = (rand01() > 0.5) ? i + 1 : i - 1;
  /* periodic boundary condition */
  j = (j + n) % n;
  // if ( j < 0 || j >= n ) return i;

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

  for ( i = 0; i < n; i++ ) {
    err += v[i] * v[i];
  }
  return sqrt(err / n);
}



#endif /* INVT_H__ */

