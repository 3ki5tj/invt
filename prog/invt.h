#ifndef INVT_H__
#define INVT_H__



/* utility functions for invt.c and predict.c */



#include <stdarg.h>
#include "mtrand.h"
#include "util.h"
#include "invtpar.h"



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



/* estimate the eigenvalues of the updating scheme
 * for a given window `win` */
static double *esteigvals(int n, int winn, const double *win)
{
  int i, j;
  double *lambda, x;

  xnew(lambda, n);

  /* loop over eigenvalues */
  for ( i = 0; i < n; i++ ) {
    lambda[i] = 1;
    /* loop over windows */
    for ( j = 1; j < winn; j++ ) {
      x = sin(i * j * M_PI * 0.5 / n);
      lambda[i] -= 4 * win[j] * x * x;
    }
  }

  return lambda;
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




/* print out the values of lambda and gamma and errors */
__inline static void dumperror(int n, const double *lambda,
    const double *gamma, int nerr, ...)
{
  int i, j;
  va_list ap;
  const double *xerr;

  fprintf(stderr, "   i:   lambda        gamma    ");
  for ( j = 0; j < nerr; j++ ) {
    fprintf(stderr, " err_%-10d", j);
  }
  fprintf(stderr, "\n");

  for ( i = 1; i < n; i++ ) {
    fprintf(stderr, "%4d: %10.6f %10.3f",
        i + 1, lambda[i], gamma[i]);
    va_start(ap, nerr);
    for ( j = 0; j < nerr; j++ ) {
      xerr = va_arg(ap, const double *);
      fprintf(stderr, " %14.4e", xerr[i]);
    }
    va_end(ap);
    fprintf(stderr, "\n");
  }
}



/* estimate the square-root error
 * under a constant updating magnitude alpha
 * according to the analytical prediction */
static double esterror_eql(double alpha, int n, double *xerr,
    const double *lambda, const double *gamma)
{
  int i;
  double y, err;

  err = 0;
  for ( i = 0; i < n; i++ ) {
    y = 0.5 * alpha * lambda[i] * gamma[i];
    if ( xerr != NULL ) {
      xerr[i] = y;
    }
    err += y;
  }

  return sqrt(err);
}




/* estimate the error of the updating schedule
 * alpha(t) = c / (t + t0)
 * for a single updating mode
 * according the analytical formula
 *
 * currently, assuming equilibrium values of < x^2 >
 * of alpha(0) = c / t0 at t = 0
 * */
static double esterror_invt1(double t, double c, double a0,
    double lambda, double gamma)
{
  const double tol = 1e-6;
  double t0, r, errsat; /* errsat: saturated error */

  t0 = c / a0;
  r = lambda * c;
  errsat = 0.5 * gamma * r / (t + t0);

  if ( fabs(2 * r - 1) < tol ) {
    /* degenerate case */
    return errsat * ( log( (t + t0) / t0 ) + 1 );
  } else {
    return errsat * (1 + 1 / (r * 2 - 1)
           * ( 1 - pow(t0 / (t + t0), 2 * r - 1) ) );
  }
}



/* estimate the error of the updating schedule
 * alpha(t) = c / (t + t0)
 * according to the analytical prediction
 *
 * currently, assuming equilibrium values of < x^2 >
 * of alpha(0) = c / t0 at t = 0, so t0 = c / a0
 * */
static double esterror_invt(double t, double c, double a0,
    int n, double *xerr,
    const double *lambda, const double *gamma)
{
  int i;
  double x, err = 0;

  for ( i = 1; i < n; i++ ) {
    x = esterror_invt1(t, c, a0, lambda[i], gamma[i]);
    if ( xerr != NULL ) {
      xerr[i] = x;
    }
    err += x;
  }

  return sqrt( err );
}



/* estimate the best parameter c for the updating schedule
 * alpha(t) = c / (t + t0)
 * according to the analytical prediction
 * if t0 <= 0, t0 is set to 1 / (c a0)
 * assuming a single local minimum
 * */
static double estbestc_invt(double t, double a0,
   int n, const double *lambda, const double *gamma,
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
  el = esterror_invt(t, cl, a0, n, NULL, lambda, gamma);
  em = esterror_invt(t, cm, a0, n, NULL, lambda, gamma);
  er = esterror_invt(t, cr, a0, n, NULL, lambda, gamma);

  for ( it = 1; ; it++ ) {
    if ( verbose >= 1 ) {
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
      el = esterror_invt(t, cl, a0, n, NULL, lambda, gamma);
    } else if ( er < em && er < el ) {
      /* er is the minimal of the three, extend to the right */
      cl = cm;
      cm = cr;
      cr = cr * 2.0;
      el = em;
      em = er;
      er = esterror_invt(t, cr, a0, n, NULL, lambda, gamma);
    } else {
      /* break the loop */
      if ( cr - cl < prec ) {
        break;
      }

      /* em is the minimal of the three */
      if ( cm - cl > cr - cm ) {
        /* refine the left half */
        cn = (cl + cm) * 0.5;
        en = esterror_invt(t, cn, a0, n, NULL, lambda, gamma);
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
        en = esterror_invt(t, cn, a0, n, NULL, lambda, gamma);
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




#endif /* INVT_H__ */

