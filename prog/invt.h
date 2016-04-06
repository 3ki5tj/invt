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



#if 0
/* compute the root-mean-squared error of the modes `u`
 * This routine here is to demonstrate the normalization */
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
#endif



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
static double *geteigvals(int n,
    int winn, const double *win, int pbc,
    double tol, int *err, int verbose)
{
  int i, j, nerr = 0;
  double *lambda, x, fac;

  if ( tol <= 0 ) { /* use the default value */
    tol = 100 * DBL_EPSILON;
  }

  xnew(lambda, n);

  /* loop over eigenvalues */
  for ( i = 0; i < n; i++ ) {
    lambda[i] = 1;
    /* loop over windows */
    for ( j = 1; j < winn; j++ ) {
      x = i * j * M_PI / n;
      if ( !pbc ) {
        x *= 0.5;
      }
      x = sin(x);
      fac = 4;
      if ( pbc && j * 2 == n ) fac = 2;
      lambda[i] -= fac * win[j] * x * x;
    }
    if ( lambda[i] < -tol ) {
      nerr += 1;
      if ( verbose ) {
        fprintf(stderr, "Warning: eigenvalue %d: %g is negative\n", i + 1, lambda[i]);
      }
    } else if ( lambda[i] < tol ) {
      lambda[i] = 0;
    }
  }

  if ( err != NULL ) {
    *err = nerr;
  }

  return lambda;
}



/* save the updating window function */
__inline static int savewin(int winn, double *win,
    const char *fn)
{
  FILE *fp;
  int i;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %d\n", winn);
  for ( i = 0; i < winn; i++ ) {
    fprintf(fp, "%d %g\n", i, win[i]);
  }
  fprintf(fp, "%d 0\n", winn);
  fclose(fp);
  return 0;
}


/* save the updating window function
 * in matrix form */
__inline static int savewinmat(int winn, double *win,
    int n, int pbc, const char *fn)
{
  FILE *fp;
  int i, j, ks[3], k, l;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %d %d\n", winn, n);
  for ( i = 1; i <= n; i++ ) {
    for ( j = 1; j <= n; j++ ) {
      double y = 0;

      if ( pbc ) {
        ks[0] = i - j;
        ks[1] = i - j - n;
        ks[2] = i - j + n;
      } else {
        ks[0] = i - j;
        ks[1] = i + j - 1;
        ks[2] = i + j - 2 * n - 1;
      }
      for ( l = 0; l < 3; l++ ) {
        k = ks[l];
        if ( k < winn && k > -winn ) {
          if ( k >= 0 ) {
            y += win[k];
          } else {
            y += win[-k];
          }
        }
      }
      fprintf(fp, "%d %d %g\n", i, j, y);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}



/* modify the window function
 * such that no eigenvalue is negative */
__inline static double *trimwindow(int n,
    int *winn, double *win,
    int pbc, double tol, int verbose)
{
  const int itmax = 10000;
  int i, j, err, it;
  double lam, lamn, lammax, *lambda = NULL;
  double nwin[NBMAX + 1];

  /* A1. copy window function */
  for ( i = 0; i < *winn; i++ ) {
    nwin[i] = win[i];
  }
  for ( i = *winn; i < NBMAX; i++ ) {
    nwin[i] = 0;
  }

  if ( n > NBMAX ) {
    fprintf(stderr, "cannot use this function with %d > %d\n",
      n, NBMAX);
    exit(1);
  }

  for ( it = 0; it < itmax; it++ ) {
    /* B1. compute the eigenvalues from the window */
    lambda = geteigvals(n, *winn, nwin, pbc,
        tol, &err, it == itmax - 1);

    /* B2. change negative eigenvalues to zeros */
    lammax = 0;
    err = 0;
    for ( i = 0; i < n; i++ ) {
      if ( lambda[i] < 0 ) {
        if ( -lambda[i] > lammax ) {
          lammax = -lambda[i];
        }
        lambda[i] = 0;
        err += 1;
      }
    }

    /* B3. inversely Fourier transform to get the window function */

    if ( pbc ) {
      /* periodic boundary condition, lambda_n = 0 */
      lamn = 0;
    } else {
      /* reflective boundary condition */
      /* compute lambda_n such that win[n] == 0 */
      lamn = lambda[0];
      for ( i = 1; i < n; i++ ) {
        lam = lambda[i];
        if ( i % 2 == 0 ) {
          lamn += 2 * lam;
        } else {
          lamn -= 2 * lam;
        }
      }
      if ( n % 2 == 0 ) {
        lamn = -lamn;
      }
    }

    if ( err == 0 || it == itmax - 1 ) {
      fprintf(stderr, "it %d: lambda %g, %g, ..., %g lamn %g\n",
          it, lambda[0], lambda[1], lambda[n - 1], lamn);
    }

    for ( j = 0; j < *winn; j++ ) {
      if ( pbc ) {
        /* periodic boundary condition */
        /* win[j] = Sum { k = 0 to n - 1 } lambda_k cos( 2 k j Pi / n) ); */
        nwin[j] = 0;
        for ( i = 0; i < n; i++ ) {
          lam = lambda[i];
          nwin[j] += lam * cos( j * i * 2 * M_PI / n );
        }
        nwin[j] /= n;
      } else {
        /* reflective boundary condition */
        /* win[j] = ( (lambda_0 + lambda_n) / 2
         *        + Sum { k = 1 to n - 1 } lambda_k cos( k j Pi / n) ); */
        nwin[j] = lambda[0] * 0.5; /* lambda0[0] should be 1.0 */
        if ( j % 2 == 0 ) {
          nwin[j] += lamn * 0.5;
        } else {
          nwin[j] -= lamn * 0.5;
        }

        for ( i = 1; i < n; i++ ) {
          lam = lambda[i];
          nwin[j] += lam * cos( j * i * M_PI / n );
        }
        nwin[j] /= n;
      }

      if ( verbose > 0 && (err == 0 || it == itmax - 1) ) {
        fprintf(stderr, "it %d: win %4d: %22.10e -> %22.10e\n",
            it, j, win[j], nwin[j]);
      }
    }

    //printf("it %d, err %d\n", it, err); getchar();
    if ( err == 0 ) {
      fprintf(stderr, "it %d: %d positive eigenvalues\n",
          it, n);
      break;
    } else {
      fprintf(stderr, "it %d: %d/%d negative eigenvalues, "
          "most negative value %g, modifying the window function\n",
          it, err, n, -lammax);
    }

    free(lambda);
  }

  /* change the window width if the iteration fails */
  if ( it >= itmax ) {
    *winn = n;
  }

  /* copy the window function */
  for ( j = 0; j < *winn; j++ ) {
    win[j] = nwin[j];
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
static double *estgamma(int n, int sampmethod,
    int pbc, double localg)
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
      x = i * M_PI / n;
      if ( !pbc ) {
        x *= 0.5;
      }
      x = tan(x);
      gamma[i] = 1.0 / (2 * localg * x * x) + 1 / (2 * localg) - 1;
    } else if ( sampmethod == SAMPMETHOD_HEATBATH ) {
      gamma[i] = 1.0;
    } else if ( sampmethod == SAMPMETHOD_MD ) {
      /* borrow the local sampling data */
      x = i * M_PI / n;
      if ( !pbc ) {
        x *= 0.5;
      }
      x = tan(x);
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
 * currently, assuming initial values of < x^2 >
 * are given by the equilibirium ones at alpha0,
 * and t0 = 2 c / alpha0.
 * */
static double esterror_invt1(double T, double c, double alpha0,
    double lambda, double gamma)
{
  const double tol = 1e-7;
  double t0, r, errsat; /* errsat: saturated error */

  t0 = 2 * c / alpha0;
  r = lambda * c;
  errsat = gamma * r / (T + t0);

  if ( lambda < 0 ) lambda = 0;

  if ( fabs(2 * r - 1) < tol ) {
    /* degenerate case */
    r += tol;
  }

  //fprintf(stderr, "r %g, 2*r-1 %g, errsat %g\n", r, 2*r - 1, errsat);
  return errsat * (r + (r - 1) * pow(t0 / (T + t0), 2 * r - 1) )
                / (r * 2 - 1);
}



/* estimate the error of the updating schedule
 * alpha(t) = c / (t + t0)
 * according to the analytical prediction
 *
 * currently, assuming initial values of < x^2 >
 * are given by the equilibrium ones at alpha0,
 * and t0 = 2 c / alpha0
 * */
static double esterror_invt(double T, double c, double a0,
    int n, double *xerr,
    const double *lambda, const double *gamma)
{
  int i;
  double x, err = 0;

  for ( i = 1; i < n; i++ ) {
    x = esterror_invt1(T, c, a0, lambda[i], gamma[i]);
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
 * assuming a single local minimum
 * */
static double estbestc_invt(double T, double a0,
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
  el = esterror_invt(T, cl, a0, n, NULL, lambda, gamma);
  em = esterror_invt(T, cm, a0, n, NULL, lambda, gamma);
  er = esterror_invt(T, cr, a0, n, NULL, lambda, gamma);

  for ( it = 1; ; it++ ) {
    if ( verbose >= 1 ) {
      fprintf(stderr, "%d: %g (%g) - %g (%g) - %g (%g)\n",
          it, cl, el, cm, em, cr, er);
      if ( verbose >= 3 ) {
        fprintf(stderr, "Press enter to continue...\n");
        getchar();
      }
    }

    /* update the bracket */
    if ( el < em && el < er ) {
      /* el is the least of the three, extend to the left */
      cr = cm;
      cm = cl;
      cl = cl * 0.5; /* make sure cl > 0 */
      er = em;
      em = el;
      el = esterror_invt(T, cl, a0, n, NULL, lambda, gamma);
    } else if ( er < em && er < el ) {
      /* er is the least of the three, extend to the right */
      cl = cm;
      cm = cr;
      cr = cr * 2.0;
      el = em;
      em = er;
      er = esterror_invt(T, cr, a0, n, NULL, lambda, gamma);
    } else {
      /* break the loop */
      if ( cr - cl < prec ) {
        break;
      }

      /* em is the least of the three */
      if ( cm - cl > cr - cm ) {
        /* refine the left side */
        cn = (cl + cm) * 0.5;
        en = esterror_invt(T, cn, a0, n, NULL, lambda, gamma);
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
        /* refine the right side */
        cn = (cm + cr) * 0.5;
        en = esterror_invt(T, cn, a0, n, NULL, lambda, gamma);
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



/* save the gamma values */
__inline static int savegamma(int n, const double *gamma,
    const char *fn)
{
  int i;
  FILE *fp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write [%s]\n", fn);
    return -1;
  }

  fprintf(fp, "# %d\n", n);
  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "%4d %14.6f\n", i, gamma[i]);
  }

  fclose(fp);

  return 0;
}


__inline static int loadgamma(int n, double *gamma,
    const char *fn)
{
  FILE *fp;
  char buf[128];
  int i, id;
  double x;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fgets(buf, sizeof buf, fp);
  sscanf(buf, "# %d", &i);
  if ( i != n ) {
    fprintf(stderr, "n mismatch, %d != %d (%s)\n",
        n, i, fn);
    fclose(fp);
    return -1;
  }

  for ( i = 0; i < n; i++ ) {
    if ( fgets(buf, sizeof buf, fp) == NULL ) {
      break;
    }
    sscanf(buf, "%d%lf", &id, &x);
    if ( i != id ) {
      fprintf(stderr, "index mismatch %d != %d (%s)\n",
          i, id, fn);
      break;
    }
    fprintf(stderr, "gamma(%d) = %g -> %g\n", i, gamma[i], x);
    gamma[i] = x;
  }

  fclose(fp);

  return 0;
}


#endif /* INVT_H__ */

