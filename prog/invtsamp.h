#ifndef INVTSAMP_H__
#define INVTSAMP_H__



/* sampling routines for invt.c */



#include "mtrand.h"
#include "util.h"
#include "invtpar.h"



/* global Metropolis move */
static int mc_metro_g(const double *v, int n, int i)
{
  int j, acc;
  double dv, r;

  /*
   * j = (int) ( (i + 1 + (n - 1) * rand01()) ) % n;
   * */
  j = (int) ( n * rand01() );
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
static int mc_metro_l(const double *v, int n, int i, int pbc)
{
  int j, acc;
  double dv, r;

  j = (rand01() > 0.5) ? i + 1 : i - 1;
  if ( pbc ) {
    j = ( j + n ) % n;
  } else {
    if ( j < 0 || j >= n ) return i;
  }

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



typedef struct {
  double x;
  double v;
  double f;
  double ek, ep;

  double dx;
  double dt;
  double tp;
  double thermdt, expndt;
  double dwa, dwb; /* potential parameters sin(x) (dwa - dwb * sin(x)) */
  double *vb; /* point to the bias potential */
} invtmd_t;



static void invtmd_init(invtmd_t *md, int n,
    double dt, double tp, double thermdt,
    double dwa, double dwb, double *vb)
{
  md->dt = dt;
  md->tp = tp;
  md->thermdt = thermdt;
  md->expndt = exp( -0.5 * thermdt );
  md->dwa = dwa;
  md->dwb = dwb;
  md->vb = vb;
  md->x = 2 * M_PI * rand01();
  md->v = randgaus();
  md->v *= sqrt(tp);
  md->f = 0;
  md->dx = 2 * M_PI / n;
  md->ep = 0;
  md->ek = 0;
}



/* compute the normal force */
static void invtmd_force(invtmd_t *md)
{
  double s = sin(md->x), c = cos(md->x);
  double a = md->dwa, b = md->dwb;

  md->ep = s * (a - b * s);
  md->f = c * (-a + 2 * b * s);
}


/* add the bias force */
static void invtmd_forceb(invtmd_t *md)
{
  int i;
  double y, x = md->x, dx = md->dx;
  double tp = md->tp;

  /* assuming x > 0 */
  i = (int) (x / dx);
  y = x - i * dx;
  if ( y >= 0.5 * dx ) {
    md->f += -(md->vb[i+1] - md->vb[i]) * tp / dx;
  } else {
    md->f += -(md->vb[i] - md->vb[i-1]) * tp / dx;
  }
}



/* velocity rescaling thermostat */
static double invtmd_vrescale(invtmd_t *md)
{
  double ek1, ek2, s, r;

  ek1 = 0.5 * md->v * md->v;
  r = randgaus();
  ek2 = ek1 + (1 - md->expndt) * ( r * r * md->tp * 0.5 - ek1 )
      + 2 * r * sqrt( md->expndt * (1 - md->expndt) * ek1 * md->tp * 0.5);
  if ( ek2 < 0 ) {
    ek2 = 0;
  }
  s = sqrt( ek2 / ek1 );
  md->v *= s;
  return ek2;
}



static int invtmd_vv(invtmd_t *md)
{
  double dt = md->dt;

  invtmd_vrescale(md);
  md->v += md->f * dt * 0.5;
  md->x = fmod(md->x + md->v * dt + 10 * 2 * M_PI, 2 * M_PI);
  invtmd_force(md);
  invtmd_forceb(md);
  md->v += md->f * dt * 0.5;
  md->ek = invtmd_vrescale(md);

  return (int) (md->x / md->dx);
}


#endif /* INVTSAMP_H__ */

