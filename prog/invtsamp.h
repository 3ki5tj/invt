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



/* Ornstein-Uhlenbeck process with variance 1/2
 * and the equilibrium distribution is exp(-x^2) / sqrt(pi) */
typedef struct {
  double x;
  double gamdt;
  double expndt;
  double sqrtdt;
  int n;
  double *v;
  double *p; /* probability distribution */
  double *ap; /* accumulative distribution */
} ouproc_t;



static ouproc_t *ouproc_open(double *v, int n, double gamdt)
{
  ouproc_t *ou;

  xnew(ou, 1);
  ou->x = randgaus();
  ou->gamdt = gamdt;
  ou->expndt = exp(-0.5 * gamdt);
  ou->sqrtdt = sqrt(2 * 0.5 * gamdt);
  ou->n = n;
  ou->v = v;
  xnew(ou->p, n);
  xnew(ou->ap, n + 1);

  return ou;
}



static void ouproc_close(ouproc_t *ou)
{
  free(ou->p);
  free(ou->ap);
  free(ou);
}



/* Ornstein-Uhlenbeck process */
static int ouproc_step(ouproc_t *ou)
{
  double r, vmin;
  int i, n = ou->n;

  /* 1. carry out the Ornstein-Uhlenbeck process */
  ou->x *= ou->expndt;
  ou->x += randgaus() * ou->sqrtdt;
  ou->x *= ou->expndt;
  r = 0.5 * (1 + erf(ou->x));

  /* 2. inversely map r to an index of v */
  /* 2a. compute the probability distribution */
  for ( vmin = ou->v[0], i = 1; i < n; i++ ) {
    if ( ou->v[i] < vmin ) {
      vmin = ou->v[i];
    }
  }
  for ( i = 0; i < n; i++ ) {
    ou->p[i] = exp( -ou->v[i] + vmin );
  }

  /* 2b. compute the accumulative distribution */
  ou->ap[0] = 0;
  for ( i = 0; i < n; i++ ) {
    ou->ap[i + 1] = ou->ap[i] + ou->p[i];
  }

  /* 2c. locate the index of the bin containing r
   * theoretically, binary search would be better
   * but since we have only O(n) operations elsewhere,
   * a linear search won't hurt too much, I suppose */
  r *= ou->ap[n];
  for ( i = n - 1; i > 0; i-- ) {
    if ( ou->ap[i] <= r )
      break;
  }

  return i;
}



typedef struct {
  double x;
  double v;
  double f;
  double ek, ep;

  int n;
  double dx;
  double dt;
  double tp;
  double thermdt;
  double expndt, sqrtdt;
#if 0
  double dwa, dwb; /* potential parameters sin(x) (dwa - dwb * sin(x)) */
#endif
  double *vb; /* point to the bias potential */
} invtmd_t;



static void invtmd_init(invtmd_t *md, int n,
    double dt, double tp, double thermdt, double *vb)
{
  md->dt = dt;
  md->tp = tp;
  md->thermdt = thermdt;

  /* parameters for integrating half thermdt */
  md->expndt = exp( -0.25 * thermdt );
  md->sqrtdt = sqrt( 2 * tp * thermdt * 0.5 );

#if 0
  md->dwa = dwa;
  md->dwb = dwb;
#endif
  md->vb = vb;
  md->x = 2 * M_PI * rand01();
  md->v = randgaus();
  md->v *= sqrt(tp);
  md->f = 0;
  md->n = n;
  md->dx = 2 * M_PI / n;
  md->ep = 0;
  md->ek = 0;
}



#if 0
/* compute the normal force */
static void invtmd_force(invtmd_t *md)
{
  double s = sin(md->x), c = cos(md->x);
  double a = md->dwa, b = md->dwb;

  md->ep = s * (a - b * s);
  md->f = c * (-a + 2 * b * s);
}
#endif


/* add the bias force */
static void invtmd_forceb(invtmd_t *md)
{
  int i, j, n = md->n;
  double y, x = md->x, dx = md->dx;
  double tp = md->tp;

  /* assuming x > 0 */
  i = (int) (x / dx);
  if ( i < 0 || i >= n ) {
    fprintf(stderr, "i %d, x %g, dx %g\n", i, x, dx);
    exit(1);
  }
  y = x - i * dx;
  if ( y >= 0.5 * dx ) {
    j = (i + 1) % n;
    md->f += -(md->vb[j] - md->vb[i]) * tp / dx;
  } else {
    j = (i - 1 + n) % n;
    md->f += -(md->vb[i] - md->vb[j]) * tp / dx;
  }
}



/* Langevin thermostat for a half step */
static double invtmd_tstat(invtmd_t *md)
{
  md->v *= md->expndt;
  md->v += randgaus() * md->sqrtdt;
  md->v *= md->expndt;
  return md->v * md->v * 0.5;
}



static int invtmd_vv(invtmd_t *md)
{
  double dt = md->dt;
  int i;

  invtmd_tstat(md);
  md->v += md->f * dt * 0.5;
  //printf("v %g, x %g, f %g, dx %g\n", md->v, md->x, md->f, md->dx);
  md->x = fmod(md->x + md->v * dt + 100 * 2 * M_PI, 2 * M_PI);
  i = (int) (md->x / md->dx);
  if ( i >= md->n || i < 0 ) {
    fprintf(stderr, "md->x %g, v %g, f %g, i %d, md->dx %g\n",
        md->x, md->v, md->f, (int)(md->x/md->dx), md->dx);
  }
  md->f = 0;
  /* invtmd_force(md); */
  invtmd_forceb(md);
  md->v += md->f * dt * 0.5;
  md->ek = invtmd_tstat(md);

  return i;
}


#endif /* INVTSAMP_H__ */

