#ifndef INVT_H__
#define INVT_H__



/* utility functions for invt.c */



#include "mtrand.h"
#include "util.h"
#include "invtpar.h"



/* multiple-bin update */
static void mbin_update(double *v, int n, int i,
    double a, const double *win, int winn, int pbc)
{
  int j, k;

  v[i] += a * win[0];
  for ( j = 1; j < winn; j++ ) {
    k = i - j;
    if ( k < 0 ) {
      k = pbc ? k + n : - k - 1;
    }
    v[k] += a * win[ j ];

    k = i + j;
    if ( k >= n ) {
      k = pbc ? k - n : 2 * n - 1 - k;
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
  static int once;

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
    } else if ( sampmethod == SAMPMETHOD_MD ) {
      /* borrow the local sampling data */
      x = tan( i * M_PI * 0.5 / n );
      gamma[i] = 1.0 / (x * x);
    } else {
      if ( !once ) { /* complain only once */
        fprintf(stderr, "Error: unknown sampling method, %d\n",
            sampmethod);
        once = 1;
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



/* estimate the square-root error
 * under a constant updating magnitude alpha
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
  if ( verbose >= 1 ) {
    fprintf(stderr, "estimated %s saturated error %g, sqr: %e\n",
        name, sqrt(err), err);
  }

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
  if ( fabs(2 * r - 1) < tol ) {
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
static double esterror_ez(double c, double t, double t0, double a0,
   int n, int winn, double *win, int sampmethod,
   int verbose)
{
  double *lamarr, *gamma, *xerr, err;

  xnew(xerr, n);

  /* estimate the eigenvalues of the w matrix,
   * for the updating scheme */
  lamarr = esteigvals(n, winn, win);

  /* estimate the autocorrelation integrals
   * of the eigenmodes of the w matrix,
   * for the updating scheme */
  gamma = estgamma(n, sampmethod);

  if ( fabs(t0) <= 0 ) {
    t0 = c / a0;
  }
  err = esterrn(1.0 / c, t, t0, n, lamarr, gamma, xerr);

  /* print out the values of lambda_i and gamma_i */
  if ( verbose >= 2 ) {
    dumplambdagamma(n, lamarr, gamma, xerr);
  }
  if ( verbose >= 1 ) {
    fprintf(stderr, "estimated final error %g, sqr: %e\n",
        sqrt(err), err);
  }

  free(lamarr);
  free(gamma);
  free(xerr);

  return sqrt(err);
}



/* estimate the best parameter c for the updating schedule
 * alpha(t) = c / (t + t0)
 * according to the analytical prediction
 * if t0 <= 0, t0 is set to 1 / (c a0)
 * assuming a single local minimum
 * */
static double estbestc(double t, double t0, double a0,
   int n, int winn, double *win, int sampmethod,
   double prec, double *err, int verbose)
{
  double cl, cm, cr, cn;
  double el, em, er, en;
  int it;

  /* set the default precision */
  if ( prec <= 0 ) {
    prec = 1e-8;
  }

  /* specify the initial bracket */
  cl = 0.5;
  cm = 1.0;
  cr = 2.0;

  /* compute the values at the initial bracket */
  el = esterror_ez(cl, t, t0, a0, n, winn, win, sampmethod, 0);
  em = esterror_ez(cm, t, t0, a0, n, winn, win, sampmethod, 0);
  er = esterror_ez(cr, t, t0, a0, n, winn, win, sampmethod, 0);

  for ( it = 1; ; it++ ) {
    if ( verbose ) {
      fprintf(stderr, "%d: %g (%g) - %g (%g) - %g (%g)\n",
          it, cl, el, cm, em, cr, er);
      if ( verbose >= 3 ) {
        fprintf(stderr, "Press enter to continue...\n");
        getchar();
      }
    }

    /* find the minimal value */
    if ( el < em && el < er ) {
      /* el is the minimal of the three, extend to the left */
      cr = cm;
      cm = cl;
      cl = cl * 0.5; /* make sure cl > 0 */
      er = em;
      em = el;
      el = esterror_ez(cl, t, t0, a0, n, winn, win, sampmethod, 0);
    } else if ( er < em && er < el ) {
      /* er is the minimal of the three, extend to the right */
      cl = cm;
      cm = cr;
      cr = cr * 2.0;
      el = em;
      em = er;
      er = esterror_ez(cr, t, t0, a0, n, winn, win, sampmethod, 0);
    } else {
      /* break the loop */
      if ( cr - cl < prec ) {
        break;
      }

      /* em is the minimal of the three */
      if ( cm - cl > cr - cm ) {
        /* refine the left half */
        cn = (cl + cm) * 0.5;
        en = esterror_ez(cn, t, t0, a0, n, winn, win, sampmethod, 0);
        if ( en > em ) {
          /* L - (N - M - R) */
          cl = cn;
          el = en;
        } else {
          /* (L - N - M) - R */
          cr = cm;
          cm = cn;
          er = em;
          em = en;
        }
      } else {
        /* refine the right half */
        cn = (cm + cr) * 0.5;
        en = esterror_ez(cn, t, t0, a0, n, winn, win, sampmethod, 0);
        if ( en > em ) {
          /* (L - M - N) - R */
          cr = cn;
          er = en;
        } else {
          /* L - (M - N - R) */
          cl = cm;
          cm = cn;
          el = em;
          em = en;
        }
      }
    }
  }

  *err = em;

  return cm;
}



/* compute the integrand for `intqt` */
static double intqt_getvel(int n, const double *lamarr, const double *gamma,
    double dq)
{
  int i;
  double lam, y = 0;

  for ( i = 1; i < n; i++ ) {
    lam = lamarr[i];
    y += gamma[i] * lam * lam * exp( 2 * lam * dq );
  }

  return sqrt( y );
}



/* compute the analytically optimal protocol
 * for a period t, and fixed
 * qt = Integral {from 0 to t} alpha(t') dt'
 * the results are saved to tarr[0..m], qarr[0..m], and aarr[0..m]
 * return the square-root error
 * */
static double intqt(double t, double qt,
    int m, double *tarr, double *aarr, double *qarr,
    int n, const double *lamarr, const double *gamma)
{
  double c, dq, err, y, dt;
  int j;

  /* integrate over q over uniform grid */
  dq = qt / m;
  tarr[0] = 0;
  qarr[0] = 0;
  y = intqt_getvel(n, lamarr, gamma, qarr[0] - qt);
  for ( j = 1; j <= m; j++ ) {
    qarr[j] = j * dq;
    tarr[j] = tarr[j - 1] + y * 0.5 * dq;
    y = intqt_getvel(n, lamarr, gamma, qarr[j] - qt);
    tarr[j] += y * 0.5 * dq;
  }

  /* normalize the time array */
  c = t / tarr[m];
  for ( j = 1; j <= m; j++ ) {
    tarr[j] *= c;
  }

  /* differentiate q(t) */
  aarr[0] = y = ( qarr[1] - qarr[0] ) / ( tarr[1] - tarr[0] );
  for ( j = 1; j < m; j++ ) {
    aarr[j] += 0.5 * y;
    y = ( qarr[j + 1] - qarr[j] ) / ( tarr[j + 1] - tarr[j] );
    aarr[j] += 0.5 * y;
  }
  aarr[m] = y;

  /* compute the error */
  err = 0;
  for ( j = 0; j <= m; j++ ) {
    y = intqt_getvel(n, lamarr, gamma, qarr[j] - qt) * aarr[j];
    y = y * y;
    if ( j == 0 ) {
      dt = tarr[1] - tarr[0];
    } else if ( j == m ) {
      dt = tarr[m] - tarr[m - 1];
    } else {
      dt = tarr[j+1] - tarr[j-1];
    }
    err += y * 0.5 * dt;
  }

  return sqrt( err );
}


/* save optimal protocol to file */
static int intqt_save(int m, const double *tarr,
    const double *aarr, const double *qarr,
    double c, double t0,
    const char *fn)
{
  FILE *fp;
  int i;
  double a1, q1;

  if ( fn == NULL || fn[0] == '\0' ) {
    return -1;
  }

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d\n", m);
  for ( i = 0; i < m; i++ ) {
    a1 = c / (tarr[i] + t0);
    q1 = c * log( 1 + tarr[i] / t0 );
    fprintf(fp, "%g\t%20.8e\t%12.6f\t%20.8e\t%12.6f\n",
        tarr[i], aarr[i], qarr[i], a1, q1);
  }

  fclose(fp);

  fprintf(stderr, "saved optimal schedule to %s\n",
     fn);

  return 0;
}



/* compute the square-root residue error */
static double intqt_reserr(double qt, double a0,
    int n, const double *lamarr, const double *gamma)
{
  int i;
  double err = 0;

  for ( i = 1; i < n; i++ ) {
    err += 0.5 * a0 * gamma[i] * lamarr[i]
         * exp(-2.0 * lamarr[i] * qt);
  }
  return sqrt(err);
}



/* return the square-root error from the optimal schedule  */
__inline static double opterror_ez(double c, double t, double a0,
    int m, const char *fnarr,
    int n, int winn, double *win, int sampmethod)
{
  double qt, t0, erra, errr, err;
  double *lamarr, *gamma;
  double *tarr, *qarr, *aarr;

  xnew(tarr, m + 1);
  xnew(aarr, m + 1);
  xnew(qarr, m + 1);

  /* estimate the eigenvalues of the w matrix,
   * for the updating scheme */
  lamarr = esteigvals(n, winn, win);

  /* estimate the autocorrelation integrals
   * of the eigenmodes of the w matrix,
   * for the updating scheme */
  gamma = estgamma(n, sampmethod);

  qt = c * log(1 + a0 * t / c);
  t0 = c / a0;

  /* compute the optimal schedule */
  erra = intqt(t, qt, m, tarr, aarr, qarr, n, lamarr, gamma);

  /* compute the residual error */
  errr = intqt_reserr(qt, a0, n, lamarr, gamma);

  /* save the schedule to file */
  intqt_save(m, tarr, aarr, qarr, c, t0, fnarr);

  free(lamarr);
  free(gamma);
  free(tarr);
  free(qarr);
  free(aarr);

  err = sqrt( erra * erra + errr * errr );
  return err;
}



#endif /* INVT_H__ */

