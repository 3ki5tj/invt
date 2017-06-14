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

  for ( i = 0; i < n; i++ ) s += v[i];
  s /= n;
  for ( i = 0; i < n; i++ ) v[i] -= s;
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



/* construct a Gaussian window
 * which may not be stable */
static void mkgauswin(double sig, int n, int pbc, double *win, int *winn)
{
  int i, j, winmax;
  double x, *lambda, lamn;

  xnew(lambda, n);
  /* starting from the Fourier transform: the eigenvalues */
  for ( i = 0; i < n; i++ ) {
    x = i * (pbc ? 2 : 1) * M_PI * sig / n;
    lambda[i] = exp(-0.5*x*x);
  }

  /* truncate the Gaussian at 10 * sigma */
  //*winn = (int) (sig * 10 + 0.5);
  //if ( *winn < 1 ) *winn = 1;

  /* limit the maximal bin width */
  winmax = pbc ? (n/2 + 1) : n;
  //if ( *winn > winmax ) *winn = winmax;
  *winn = winmax;

  /* for reflective boundary condition, compute lambda_n such that win[n] == 0 */
  for ( lamn = lambda[0], i = 1; i < n; i++ )
    lamn += lambda[i] * (i % 2 ? -2 : 2);
  if ( n % 2 == 0 ) lamn = -lamn;

  for ( j = 0; j < *winn; j++ ) {
    if ( pbc ) {
      /* periodic boundary condition */
      /* win[j] = Sum { k = 0 to n - 1 } lambda_k cos( 2 k j Pi / n) ); */
      for ( win[j] = 0, i = 0; i < n; i++ )
        win[j] += lambda[i] * cos( j * i * 2 * M_PI / n );
    } else {
      /* reflective boundary condition */
      /* win[j] = ( (lambda_0 + lambda_n) / 2
       *        + Sum { k = 1 to n - 1 } lambda_k cos( k j Pi / n) ); */
      win[j] = (lambda[0] + lamn * (j % 2 ? -1 : 1)) * 0.5;
      for ( i = 1; i < n; i++ )
        win[j] += lambda[i] * cos(j * i * M_PI / n);
    }
    win[j] /= n;
    //x = j / sig; printf("win %4d: %g %g\n", j, win[j], win[j]/(exp(-0.5*x*x)/sig/sqrt(2*M_PI)));
  }
  free(lambda);
}



/* construct the window for an optimal updating scheme */
static void mksincwin(int km, int n, int pbc, double *win, int *winn)
{
  int i;
  double k2m1 = km * 2 - 1, ang, y;

  *winn = pbc ? (n/2 + 1) : n;

  if ( pbc ) {
    /* win[i] = sin((2*km-1)*i*Pi/n)/sin(i*Pi/n)/n; */
    win[0] = 1.0 * k2m1 / n;
    for ( i = 1; i < *winn; i++ ) {
      ang = i * M_PI / n;
      y = sin(k2m1 * ang);
      win[i] = y / sin(ang) / n;
    }
  } else {
    /* win[i] = {(-1)^(n+km+i)+
     *         sin[(2*km-1)*i*Pi/2/n]/sin(i*Pi/2/n)}/(2*n); */
    win[0] = ((((n + km) % 2) ? -1 : 1) + k2m1)/ (2.0 * n);
    for ( i = 1; i < n; i++ ) {
      ang = i * M_PI / (2 * n);
      y = sin(k2m1 * ang);
      win[i] = ( ((n + km + i) % 2 ? -1 : 1)
             + y / sin(ang) ) / (2 * n);
    }
  }
}



/* estimate the eigenvalues of the updating scheme
 * for a given window `win`
 * `tol`: tolerance of negative eigenvalues
 * `err`: number of negative eigenvalues */
static void geteigvals(double *lambda, int n,
    const double *win, int winn, int pbc,
    double tol, int *err, int verbose)
{
  int i, j, nerr = 0;
  double x, fac;

  /* use the default tolerance of negative eigenvalues */
  if ( tol <= 0 ) tol = 100 * DBL_EPSILON;

  /* loop over eigenvalues */
  for ( i = 0; i < n; i++ ) {
    lambda[i] = 1;
    /* loop over kernel indices
     * lambda_i = 1 - 4 Sum_{j=1}^b mu_j sin^2(i*j*PI/(gn))
     * where g == 1 for periodic or g == 2 for non-periodic variable */
    for ( j = 1; j < winn; j++ ) {
      x = i * j * M_PI / n;
      if ( !pbc ) x *= 0.5;
      x = sin(x);
      fac = ( pbc && j * 2 == n ) ? 2 : 4;
      lambda[i] -= fac * win[j] * x * x;
    }
    if ( lambda[i] < -tol ) {
      nerr += 1;
      if ( verbose ) {
        fprintf(stderr, "Warning: eigenvalue %d: %g is negative\n", i + 1, lambda[i]);
      }
    } else if ( lambda[i] < tol ) { // necessary? for sinc only?
      lambda[i] = 0;
    }
  }

  if ( err != NULL ) *err = nerr;
}



/* save the updating window function */
__inline static int savewin(double *win, int winn,
    const char *fn)
{
  FILE *fp;
  int i;
  double ymin = 0, ymax = 0;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %d\n", winn);
  for ( i = 0; i < winn; i++ ) {
    if ( win[i] < ymin ) ymin = win[i];
    if ( win[i] > ymax ) ymax = win[i];
    fprintf(fp, "%d %g\n", i, win[i]);
  }
  fprintf(fp, "%d 0\n", winn);
  fclose(fp);
  fprintf(stderr, "saved window function to %s, min %g, max %g\n",
      fn, ymin, ymax);
  return 0;
}



/* save the updating window function
 * in matrix form */
__inline static int savewinmat(double *win, int winn,
    int n, int pbc, const char *fn)
{
  FILE *fp;
  int i, j, ks[3], k, l;
  double y, ymin = 0, ymax = 0;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fprintf(fp, "# %d %d\n", winn, n);
  for ( i = 1; i <= n; i++ ) {
    for ( j = 1; j <= n; j++ ) {
      y = 0;
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
      if ( y > ymax ) ymax = y;
      if ( y < ymin ) ymin = y;
      fprintf(fp, "%9.6f ", y);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  fprintf(stderr, "saved updating matrix to %s, min %g, max %g\n",
      fn, ymin, ymax);
  return 0;
}



/* modify the window function
 * such that no eigenvalue in `lambda` is negative */
__inline static void stablizewin(double *lambda, int n,
    double *win, int *winn,
    int pbc, double tol, int verbose)
{
  const int itmax = 10000;
  int i, j, err, it;
  double lamn = 0, lammax, *nwin;

  xnew(nwin, n + 1);

  /* A1. copy window function */
  for ( i = 0; i < n; i++ )
    nwin[i] = (i < *winn ? win[i] : 0);

  for ( it = 0; it < itmax; it++ ) {
    /* B1. compute the eigenvalues from the window */
    geteigvals(lambda, n, nwin, *winn, pbc,
        tol, &err, it == itmax - 1);

    /* B2. change negative eigenvalues to zeros */
    lammax = 0; /* magnitude of the most negative eigenvalue */
    err = 0;
    for ( i = 0; i < n; i++ ) {
      if ( lambda[i] < -1e3 * DBL_EPSILON ) {
        if ( -lambda[i] > lammax )
          lammax = -lambda[i];
        lambda[i] = 0;
        err += 1;
      }
    }

    /* B3. inversely Fourier transform to get the window function */

    if ( !pbc ) { /* reflective boundary condition */
      /* compute lambda_n such that win[n] == 0 */
      for ( lamn = lambda[0], i = 1; i < n; i++ )
        lamn += lambda[i] * (i % 2 ? -2 : 2);
      if ( n % 2 == 0 ) lamn = -lamn;
    }

    if ( err == 0 || it == itmax - 1 ) {
      fprintf(stderr, "it %d: lambda %g, %g, ..., %g lamn %g\n",
          it, lambda[0], lambda[1], lambda[n - 1], lamn);
      if ( err == 0 ) break;
    }

    /* construct a new window */
    for ( j = 0; j < *winn; j++ ) {
      if ( pbc ) {
        /* periodic boundary condition */
        /* win[j] = Sum { k = 0 to n - 1 } lambda_k cos( 2 k j Pi / n) ); */
        nwin[j] = 0;
        for ( i = 0; i < n; i++ )
          nwin[j] += lambda[i] * cos( j * i * 2 * M_PI / n );
      } else {
        /* reflective boundary condition */
        /* win[j] = ( (lambda_0 + (-1)^j lambda_n) / 2
         *        + Sum { k = 1 to n - 1 } lambda_k cos( k j Pi / n) ); */
        nwin[j] = (lambda[0] + lamn * (j % 2 ? -1: 1)) * 0.5;
        for ( i = 1; i < n; i++ )
          nwin[j] += lambda[i] * cos(j * i * M_PI / n);
      }
      nwin[j] /= n;

      if ( verbose > 0 && (err == 0 || it == itmax - 1) ) {
        fprintf(stderr, "it %d: win %4d: %22.10e -> %22.10e\n",
            it, j, win[j], nwin[j]);
      }
    }

    if ( err == 0 ) {
      fprintf(stderr, "it %d: %d positive eigenvalues\n",
          it, n);
      break;
    } else {
      fprintf(stderr, "it %d: %d/%d negative eigenvalues, "
          "most negative value %g, modifying the window function\n",
          it, err, n, -lammax);
    }
  }

  /* change the window width if the iteration fails */
  if ( it >= itmax ) *winn = n;

  /* copy the window function */
  for ( j = 0; j < *winn; j++ )
    win[j] = nwin[j];

  free(nwin);
}



/* prepare the window function */
__inline static double *prepwin(double *lambda, int n,
    const double *win0, int winn0, double *win, int *winn,
    int pbc, double gaussig, int kc,
    const char *fnwin, const char *fnwinmat, int verbose)
{
  int i;

  if ( gaussig > 0 ) {
    mkgauswin(gaussig, n, pbc, win, winn);
  } else if ( kc > 0 ) {
    mksincwin(kc, n, pbc, win, winn);
  } else {
    /* copy the user window */
    *winn = winn0;
    for ( i = 0; i < winn0; i++ )
      win[i] = win0[i];
  }

  /* modify the window function such that all eigenvalues
   * lambda[i] are positive-definite */
  stablizewin(lambda, n, win, winn, pbc, 0, verbose);
  if ( fnwin[0] != '\0' ) {
    /* save the window kernel */
    savewin(win, *winn, fnwin);
  }
  if ( fnwinmat[0] != '\0' ) {
    /* save the n x n updating matrix */
    savewinmat(win, *winn, n, pbc, fnwinmat);
  }
  return win;
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
__inline static void estgamma(double *gamma, int n,
    int sampmethod, int pbc, double mvsize)
{
  int i;
  double x;
  static int once;

  if ( sampmethod == SAMPMETHOD_MD ) {
    sampmethod = SAMPMETHOD_METROLOCAL;
  }
  gamma[0] = 0;
  for ( i = 1; i < n; i++ ) {
    if ( sampmethod == SAMPMETHOD_METROGLOBAL ) {
      gamma[i] = 1.0; /* n / (n - 1.0); */
    } else if ( sampmethod == SAMPMETHOD_METROLOCAL ) {
      x = i * M_PI / n;
      if ( !pbc ) x *= 0.5;
      x = tan(x);
      gamma[i] = 1.0 / (mvsize * x * x) + 1 / mvsize - 1;
    } else if ( sampmethod == SAMPMETHOD_GAUSS ) {
      x = 2 * mvsize * i * M_PI / n;
      if ( !pbc ) x *= 0.5;
      x = exp(-x*x/2);
      gamma[i] = (1 + x) / (1 - x);
    } else if ( sampmethod == SAMPMETHOD_HEATBATH ) {
      gamma[i] = 1.0;
    } else {
      if ( !once ) { /* complain only once */
        fprintf(stderr, "Error: unknown sampling method, %d\n",
            sampmethod);
        once = 1;
      }
      gamma[i] = 1.0;
    }
  }
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
__inline static double esterror_eql(double alpha, int n, double *xerr,
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
 * and the default t0 = 2 / alpha0.
 * */
static double esterror_invt1(double T, double c, double a0,
    double t0, double er0, double *ea, double *er,
    double lambda, double gamma)
{
  const double tol = 1e-7;
  double f, y, r, errsat; /* errsat: saturated error */

  r = lambda * c;

  if ( lambda < 0 ) lambda = 0;

  if ( fabs(2 * r - 1) < tol ) {
    /* degenerate case */
    r += tol;
  }

  f = t0/(T + t0);
  y = pow(f, r * 2 - 1);
  errsat = gamma * r / (T + t0);

  *er = (0.5 * a0 * gamma * lambda + er0) * f * y;
  *ea = errsat * r / (r * 2 - 1) * (1 - y);
  //fprintf(stderr, "error %g vs %g\n", er + ea, errsat * (r + (r - 1) * pow(t0 / (T + t0), 2 * r - 1) ) / (r * 2 - 1)); getchar();
  //fprintf(stderr, "r %g, 2*r-1 %g, errsat %g\n", r, 2*r - 1, errsat);
  //return errsat * (r + (r - 1) * pow(t0 / (T + t0), 2 * r - 1) ) / (r * 2 - 1);
  return *er + *ea;
}



/* estimate the error of the updating schedule
 * alpha(t) = c / (t + t0)
 * according to the analytical prediction
 *
 * currently, assuming initial values of < x^2 >
 * are given by the equilibrium ones at a0,
 * */
#define esterror_invt(T, c, a0, t0, n, xerr, lambda, gamma) \
    esterror_invt_x(T, c, a0, t0, n, NULL, NULL, \
        NULL, xerr, NULL, NULL, lambda, gamma)

static double esterror_invt_x(double T, double c, double a0,
    double t0, int n, double *errr, double *erra,
    double *xerr0, double *xerr, double *xerrr, double *xerra,
    const double *lambda, const double *gamma)
{
  int i;
  double x, xr0, xr, xa, E = 0, Er = 0, Ea = 0;

  if ( t0 <= 0 ) t0 = 2 / a0;

  for ( i = 1; i < n; i++ ) {
    xr0 = ( xerr0 != NULL ? xerr0[i] : 0 );
    x = esterror_invt1(T, c, a0, t0, xr0, &xr, &xa, lambda[i], gamma[i]);
    if ( xerr  != NULL ) xerr[i] = x;
    if ( xerrr != NULL ) xerrr[i] = xr;
    if ( xerra != NULL ) xerra[i] = xa;
    E += x;
    Er += xr;
    Ea += xa;
  }

  if ( errr != NULL ) *errr = sqrt( Er );
  if ( erra != NULL ) *erra = sqrt( Ea );
  return sqrt( E );
}



/* estimate the best parameter c for the updating schedule
 * alpha(t) = c / (t + t0)
 * according to the analytical prediction
 * assuming a single local minimum
 * */
__inline static double estbestc_invt(double T, double a0, double t0,
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
  el = esterror_invt(T, cl, a0, t0, n, NULL, lambda, gamma);
  em = esterror_invt(T, cm, a0, t0, n, NULL, lambda, gamma);
  er = esterror_invt(T, cr, a0, t0, n, NULL, lambda, gamma);

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
      el = esterror_invt(T, cl, a0, t0, n, NULL, lambda, gamma);
    } else if ( er < em && er < el ) {
      /* er is the least of the three, extend to the right */
      cl = cm;
      cm = cr;
      cr = cr * 2.0;
      el = em;
      em = er;
      er = esterror_invt(T, cr, a0, t0, n, NULL, lambda, gamma);
    } else {
      /* break the loop */
      if ( cr - cl < prec ) {
        break;
      }

      /* em is the least of the three */
      if ( cm - cl > cr - cm ) {
        /* refine the left side */
        cn = (cl + cm) * 0.5;
        en = esterror_invt(T, cn, a0, t0, n, NULL, lambda, gamma);
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
        en = esterror_invt(T, cn, a0, t0, n, NULL, lambda, gamma);
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
  double gamtot = 0;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write [%s]\n", fn);
    return -1;
  }

  /* compute the total gamma value */
  for ( i = 1; i < n; i++ ) gamtot += gamma[i];

  fprintf(fp, "# %d %g\n", n, gamtot);
  for ( i = 1; i < n; i++ ) {
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
  double x, sum = 0;

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

  gamma[0] = 0;
  for ( i = 1; i < n; i++ ) {
    if ( fgets(buf, sizeof buf, fp) == NULL ) {
      break;
    }
    sscanf(buf, "%d%lf", &id, &x);
    if ( i != id ) {
      fprintf(stderr, "index mismatch %d != %d (%s)\n",
          i, id, fn);
      return -1;
    }
    //fprintf(stderr, "gamma(%d) = %g -> %g\n", i, gamma[i], x);
    gamma[i] = x;
    sum += x;
  }
  fprintf(stderr, "loaded %d gamma values %g,...%g, sum %g, from %s\n",
      i-1, gamma[1], gamma[i-1], sum, fn);
  fclose(fp);

  return 0;
}


#endif /* INVT_H__ */

