#ifndef INVT_H__
#define INVT_H__



/* utility functions for invt.c */



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
static double geterror(double *v, int n, const double *p)
{
  double err = 0, x;
  int i;

  /* subtract the baseline */
  shift(v, n);

  for ( i = 0; i < n; i++ ) {
    x = v[i] * v[i];
    /* multiply the normalization factor */
    x *= ( p != NULL ) ? p[i] : 1.0 / n;
    err += x;
  }

  return sqrt(err);
}



/* compute the root-mean-squared error of the modes `u`
 * This routine here is demonstrate the normalization */
__inline static double getuerror(double *u, int n)
{
  double err = 0, x;
  int i;

  for ( i = 0; i < n; i++ ) {
    x = u[i] * u[i];
    err += x;
  }

  return sqrt(err);
}



/* normalize the error such that the standard deviation is sig */
__inline static double normalize(double *v, int n, double sig,
    const double *p)
{
  double err;
  int i;

  err = geterror(v, n, p);
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
      s += costab[k*n + i] * v[i];
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



/* estimate the eigenvalues of the updating scheme
 * for a given window `win` */
static double *esteigvals(int n, int winn, const double *win)
{
  int i, j;
  double *lamarr, x;

  xnew(lamarr, n);

  /* loop over eigenvalues */
  for ( i = 0; i < n; i++ ) {
    lamarr[i] = 1;
    /* loop over windows */
    for ( j = 1; j < winn; j++ ) {
      x = sin(i * j * M_PI * 0.5 / n);
      lamarr[i] -= 4 * win[j] * x * x;
    }
  }

  return lamarr;
}



/* estimate the correlation integrals, Gamma_i,
 * of the eigenmodes of the updating scheme
 * for a given sampling method `sampmethod`
 *
 * Note: this function assumes that the eigenmodes of
 * the updating scheme (eigenvectors of w) are the same
 * as those of the Monte Carlo transition matrix
 * This is true only for perfect sampling and
 * local Monte Carlo sampling */
static double *estgamma(int n, int sampmethod)
{
  int i;
  double *gamma, x;

  xnew(gamma, n);

  gamma[0] = 0;
  for ( i = 1; i < n; i++ ) {
    if ( sampmethod == SAMPMETHOD_METROGLOBAL ) {
      gamma[i] = 1.0; /* n / (n - 1.0); */
    } else if ( sampmethod == SAMPMETHOD_METROLOCAL ) {
      x = tan( i * M_PI * 0.5 / n );
      gamma[i] = 1.0 / (x * x);
    } else if ( sampmethod == SAMPMETHOD_HEATBATH ) {
      gamma[i] = 1.0;
    } else {
      if ( i == 1 ) { /* complain only once */
        fprintf(stderr, "Error: unknown sampling method, %d\n",
            sampmethod);
      }
      gamma[i] = 1.0;
    }
  }

  return gamma;
}




/* print out the values of lambda_i and gamma_i */
static void dumplambdagamma(int n, const double *lamarr,
    const double *gamma, const double *xerr)
{
  int i;

  fprintf(stderr, "   i:   lambda_i      gamma_i\n");
  for ( i = 0; i < n; i++ ) {
    fprintf(stderr, "%4d: %10.6f %10.3f %20.6e\n",
        i + 1, lamarr[i], gamma[i], xerr[i]);
  }
}



/* estimate the error of a constant alpha
 * according to the analytical prediction */
static double esterror0_ez(double alpha,
   int n, int winn, double *win, int sampmethod,
   const char *name, int verbose)
{
  int i;
  double *lamarr, *gamma, *xerr, err;

  xnew(xerr, n);

  /* estimate the eigenvalues of the w matrix,
   * for the updating scheme */
  lamarr = esteigvals(n, winn, win);

  /* estimate the correlation integrals
   * of the eigenmodes of the w matrix,
   * for the updating scheme */
  gamma = estgamma(n, sampmethod);

  err = 0;
  for ( i = 0; i < n; i++ ) {
    xerr[i] = 0.5 * alpha * lamarr[i] * gamma[i];
    err += xerr[i];
  }

  /* print out the values of lambda_i and gamma_i */
  if ( verbose >= 2 ) {
    dumplambdagamma(n, lamarr, gamma, xerr);
  }
  fprintf(stderr, "estimated %s saturated error %g, sqr: %e\n",
      name, sqrt(err), err);

  free(lamarr);
  free(gamma);
  free(xerr);

  return sqrt(err);
}



/* estimate the error of the updating schedule
 * alpha(t) = 1/ [lambda (t + t0)]
 * for a single updating mode
 * according the analytical formula
 *
 * currently, assuming equilibrium values of < x^2 >
 * of alpha(0) = c / t0 at t = 0
 * */
static double esterr1(double lambda, double t, double t0,
    double lambda_i, double gamma_i)
{
  const double tol = 1e-6;
  double r, errsat; /* errsat: saturated error */

  r = lambda_i / lambda;
  errsat = 0.5 * gamma_i * r / (t + t0);
  /* degenerate case */
  if ( r < tol ) {
    return errsat * ( log( (t + t0) / t0 ) + 1 );
  } else {
    return errsat * (1 + 1 / (r * 2 - 1)
           * ( 1 - pow(t0 / (t + t0), 2 * r - 1) ) );
  }
}



/* estimate the error of the updating schedule
 * alpha(t) = 1/ [lambda (t + t0)]
 * according to the analytical prediction
 *
 * currently, assuming equilibrium values of < x^2 >
 * of alpha(0) = c / t0 at t = 0
 * */
static double esterrn(double lambda, double t, double t0,
    int n, const double *lamarr, const double *gamma,
    double *xerr)
{
  int i;
  double x, err = 0;

  for ( i = 1; i < n; i++ ) {
    x = esterr1(lambda, t, t0, lamarr[i], gamma[i]);
    if ( xerr != NULL ) {
      xerr[i] = x;
    }
    err += x;
  }

  return err;
}



/* estimate the final error of the updating schedule
 * alpha(t) = c / (t + t0)
 * according to the analytical prediction
 *
 * currently, assuming equilibrium values of < x^2 >
 * of alpha(0) = c / t0 at t = 0
 * */
static double esterror_ez(double c, double t, double t0,
   int n, int winn, double *win, int sampmethod,
   int verbose)
{
  double *lamarr, *gamma, *xerr, err;

  xnew(xerr, n);

  /* estimate the eigenvalues of the w matrix,
   * for the updating scheme */
  lamarr = esteigvals(n, winn, win);

  /* estimate the correlation integrals
   * of the eigenmodes of the w matrix,
   * for the updating scheme */
  gamma = estgamma(n, sampmethod);

  err = esterrn(1.0 / c, t, t0, n, lamarr, gamma, xerr);

  /* print out the values of lambda_i and gamma_i */
  if ( verbose >= 2 ) {
    dumplambdagamma(n, lamarr, gamma, xerr);
  }
  fprintf(stderr, "estimated final error %g, sqr: %e\n",
      sqrt(err), err);

  free(lamarr);
  free(gamma);
  free(xerr);

  return sqrt(err);
}



#endif /* INVT_H__ */

