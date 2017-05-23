#ifndef METAD_H__
#define METAD_H__


#include "util.h"
#include "cosmodes.h"
#include "eig.h"
#include "intq.h"
#include "invt.h"
#include "cmvar.h" /* for autocorrelation integrals (gamma) */



const double metad_lamcut = 0.1;

/* metadynamics for integer */
typedef struct {
  int imin, imax, idel;
  double xmin, xmax, xdel;
  int n;
  double *v; /* bias potential */
  double *h; /* histogram */
  double *vref; /* reference bias potential */
  double a; /* updating magnitude */
  int pbc; /* periodic boundary condition */
  int winn;
  double *win;
  double *lambda; /* eigenvalues of the updating marix */
  cmvar_t *cm;
  double *gamma; /* correlation integrals */
  double *tmat; /* n x n transition matrix */
  double *tgamma; /* correlation integrals from the transition matrix */
  double *costab; /* cosine transform coefficients */
  double hfl;
  double *vtmp;
  double *hmod;
  intq_t *intq;
  double errref; /* estimated error */
} metad_t;



/* prepare the window function */
static void metad_prepwin(metad_t *metad,
    double gaussig, int okmax, const double *win0, int winn0)
{
  int i, n = metad->n, pbc = metad->pbc;

  xnew(metad->win, n);
  if ( gaussig > 0 ) {
    mkgauswin(gaussig, n, pbc, metad->win, &metad->winn);
  } else if ( okmax >= 0 ) {
    mksincwin(okmax, n, pbc, metad->win, &metad->winn);
  } else {
    /* copy the user window */
    metad->winn = winn0;
    for ( i = 0; i < metad->winn; i++ )
      metad->win[i] = win0[i];
  }

  /* modify the window function such that all eigenvalues
   * lambda[i] are positive-definite */
  stablizewin(metad->lambda, n, metad->win, &metad->winn, pbc, 0.0, 1);
  /* save the window kernel */
  savewin(metad->win, metad->winn, "ker.dat");
  /* save the n x n updating matrix */
  savewinmat(metad->win, metad->winn, n, pbc, "win.dat");
}

static metad_t *metad_open(int imin, int imax, int idel,
    int pbc, double gaussig, int okmax, const double *win, int winn)
{
  metad_t *metad;
  int i, n;

  xnew(metad, 1);
  metad->imin = imin;
  metad->imax = imax;
  metad->idel = idel;
  metad->n = n = (imax - imin) / idel + 1;
  metad->imax = metad->imin + n * metad->idel;
  metad->xmin = 0;
  metad->xmax = 0;
  metad->xdel = 0;
  xnew(metad->v, n);
  xnew(metad->h, n);
  xnew(metad->vref, n);
  for ( i = 0; i < n; i++ ) {
    metad->v[i] = 0;
    metad->h[i] = 0;
    metad->vref[i] = 0;
  }
  metad->a = 1.0 / n;
  metad->pbc = pbc;
  /* prepare the window */
  xnew(metad->lambda, n);
  metad_prepwin(metad, gaussig, okmax, win, winn);
  metad->cm = cmvar_open(metad->n, metad->pbc);
  xnew(metad->gamma, n);
  xnew(metad->tmat, n * n);
  xnew(metad->tgamma, n);
  metad->costab = mkcostab(n, metad->pbc);
  xnew(metad->vtmp, n);
  xnew(metad->hmod, n);
  metad->intq = NULL;
  metad->errref = 0;
  return metad;
}



__inline static metad_t *metad_openf(double xmin, double xmax, double xdel,
    int pbc, double gaussig, int okmax, const double *win, int winn)
{
  metad_t *metad;
  int n = (int)((xmax - xmin) / xdel);
  gaussig /= xdel;
  metad = metad_open(0, n - 1, 1, pbc, gaussig, okmax, win, winn);
  metad->xmin = xmin;
  metad->xdel = xdel;
  metad->xmax = xmin + metad->n * xdel;
  return metad;
}



static void metad_close(metad_t *metad)
{
  free(metad->v);
  free(metad->h);
  free(metad->win);
  free(metad->lambda);
  cmvar_close(metad->cm);
  free(metad->gamma);
  free(metad->tmat);
  free(metad->tgamma);
  free(metad->costab);
  free(metad->vtmp);
  free(metad->hmod);
  if ( metad->intq != NULL ) {
    intq_close( metad->intq );
  }
  free(metad);
}



/* return the bin index of value x */
__inline static int metad_getindex(metad_t *metad, int x)
{
  if ( x < metad->imin || x >= metad->imax )
    return -1;
  return (x - metad->imin) / metad->idel;
}

__inline static int metad_getindexf(metad_t *metad, double x)
{
  if ( x < metad->xmin || x >= metad->xmax )
    return -1;
  return (int) ((x - metad->xmin) / metad->xdel);
}



/* decide if a transition from xold (index iold)
 * to xnew (index iold) is to be accepted */
__inline static int metad_acc(metad_t *metad,
    int iold, int xnew, int *inew)
{
  double dv, r;

  *inew = metad_getindex(metad, xnew);
  /* reject out-of-boundary moves */
  if ( *inew < 0 ) return 0;

  if ( iold == *inew ) return 1;
  dv = metad->v[*inew] - metad->v[iold];
  if ( dv <= 0 ) return 1;
  r = rand01();
  return r < exp(-dv);
}


/* compute the histogram fluctuation (Wang-Landau version) */
__inline static double metad_hfl_wl(metad_t *metad)
{
  int i;
  double hmin, hmax, hi;

  hmin = hmax = metad->h[0];
  for ( i = 1; i < metad->n; i++ ) {
    hi = metad->h[i];
    if ( hi < hmin ) {
      hmin = hi;
    } else if ( hi > hmax ) {
      hmax = hi;
    }
  }

  return (hmax - hmin) / (hmax + hmin + DBL_EPSILON);
}

/* compute the histogram fluctuation */
__inline static double metad_hfl(metad_t *metad)
{
  int i, k, n = metad->n;
  double x, fl = 0, tot = 0;

  for ( i = 0; i < n; i++ ) tot += metad->h[i];
  getcosmodes(metad->h, n, metad->vtmp, metad->costab);
  /* truncate modes with lambda */
  for ( k = 1; k < n; k++ ) {
    if ( metad->lambda[k] < metad_lamcut ) break;
    x = n * metad->vtmp[k] / tot;
    fl += x * x;
  }
  /* compute the filtered histogram */
  for ( ; k < n; k++ ) metad->vtmp[k] = 0;
  fromcosmodes(metad->hmod, n, metad->vtmp, metad->costab);
  return sqrt( fl );
}



/* check if histogram is flat enough to switch to
 * a smaller updating magnitude */
static int metad_wlcheck(metad_t *metad, double fl, double magred)
{
  int i;

  /* compute the histogram fluctuation */
  metad->hfl = metad_hfl(metad);
  //printf("fl %g, %g\n", metad->hfl, metad_hfl_wl(metad));
  //getchar();
  /* return if the histogram not flat enough */
  if ( metad->hfl > fl ) return 0;

  /* reduce the updating magnitude and clear the histogram */
  metad->a *= magred;
  for ( i = 0; i < metad->n; i++ ) metad->h[i] = 0;
  fprintf(stderr, "changing the updating magnitude to %g\n", metad->a);
  return 1;
}



__inline static void metad_updatev_wl(metad_t *metad, int i)
{
  double amp = metad->n * metad->a;
  metad->v[i] += amp;
  metad->h[i] += 1;
}

/* multiple-bin update */
__inline static void metad_updatev(metad_t *metad, int i)
{
  int j, k, n = metad->n;
  double amp = metad->n * metad->a;

  metad->v[i] += amp * metad->win[0];
  /* update the bias potential at the neighbors */
  for ( j = 1; j < metad->winn; j++ ) {
    /* left neighbor */
    k = i - j;
    if ( k < 0 ) k = metad->pbc ? k + n : - k - 1;
    metad->v[k] += amp * metad->win[j];
    if ( j * 2 == n && metad->pbc ) continue;

    /* right neighbor */
    k = i + j;
    if ( k >= n ) k = metad->pbc ? k - n : 2 * n - 1 - k;
    metad->v[k] += amp * metad->win[j];
  }
  metad->h[i] += 1;
}


static void metad_trimv(metad_t *metad, double *v)
{
  int i, n = metad->n;
  double v0 = 0;

  for ( i = 0; i < n; i++ ) v0 += v[i];
  v0 /= n;
  for ( i = 0; i < n; i++ ) v[i] -= v0;
}


__inline static void metad_saveheader(metad_t *metad, FILE *fp)
{
  if ( fabs(metad->xmax - metad->xmin) > 0 ) {
    fprintf(fp, "# %d %g %g %g",
        metad->n, metad->xmin, metad->xmax, metad->xdel);
  } else {
    fprintf(fp, "# %d %d %d %d",
        metad->n, metad->imin, metad->imax, metad->idel);
  }
  fprintf(fp, " %g %d\n", metad->a, metad->pbc);
}

__inline static double metad_getx(metad_t *metad, int i)
{
  if (fabs(metad->xmax - metad->xmin) > 0) { /* float */
    return metad->xmin + (i + 0.5) * metad->xdel;
  } else {
    return metad->imin + i * metad->idel;
  }
}

/* compute the mode-trunction error */
__inline static double metad_errtrunc(metad_t *metad,
    double *v, int *kc, const char *fntrunc)
{
  int k, i, n = metad->n;
  double x, vi, err = 0;
  FILE *fp;

  getcosmodes(v, n, metad->vtmp, metad->costab);
  *kc = -1;
  for ( k = 1; k < n; k++ ) {
    if ( metad->lambda[k] > 0.5 ) continue;
    if ( *kc < 0 ) *kc = k;
    err += metad->vtmp[k] * metad->vtmp[k];
  }

  /* save mode-truncated profile */
  if ( fntrunc != NULL && (fp=fopen(fntrunc, "w")) != NULL ) {
    metad_saveheader(metad, fp);
    for ( i = 0; i < n; i++ ) {
      for ( vi = 0, k = 0; k < *kc; k++ )
        vi += metad->costab[k*n + i] * metad->vtmp[k];
      x = metad_getx(metad, i);
      fprintf(fp, "%g %g %g\n", x, vi, v[i]);
    }
    fclose(fp);
  }

  return err;
}

/* save the bias potential to file */
static int metad_save(metad_t *metad, const char *fn)
{
  FILE *fp;
  int i;
  double x;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  metad_trimv(metad, metad->v);
  metad_saveheader(metad, fp);
  for ( i = 0; i < metad->n; i++ ) {
    x = metad_getx(metad, i);
    fprintf(fp, "%g %g %g %g %g\n", x,
        metad->v[i], metad->h[i], metad->vref[i], metad->hmod[i]);
  }
  fclose(fp);
  return 0;
}



/* load the bias potential from file */
__inline static int metad_load(metad_t *metad, double *v, const char *fn)
{
  FILE *fp;
  int i, isfloat = (fabs(metad->xmax - metad->xmin) > 0), err = -1;
  int n, imin, imax, idel;
  double x, xmin, xmax, xdel;
  static char buf[256];

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fgets(buf, sizeof buf, fp);
  if ( isfloat ) {
    sscanf(buf, "# %d%lf%lf%lf", &n, &xmin, &xmax, &xdel);
    if ( n != metad->n
      || fabs(xmin - metad->xmin) > 1e-6
      || fabs(xmax - metad->xmax) > 1e-6
      || fabs(xdel - metad->xdel) > 1e-6 ) {
      fprintf(stderr, "%s data mismatch\n", fn);
      goto ERR;
    }
  } else {
    sscanf(buf, "# %d%d%d%d", &n, &imin, &imax, &idel);
    if ( n != metad->n
      || imin != metad->imin
      || imax != metad->imax
      || idel != metad->idel ) {
      fprintf(stderr, "%s data mismatch\n", fn);
      goto ERR;
    }
  }
  for ( i = 0; i < metad->n; i++ ) {
    fgets(buf, sizeof buf, fp);
    if ( 2 != sscanf(buf, "%lf%lf", &x, &v[i]) ) {
      fprintf(stderr, "%s corrupted at line %d\n", fn, i);
      goto ERR;
    }
  }
  metad_trimv(metad, v);
  err = 0;
ERR:
  fclose(fp);
  return err;
}



/* update the accumulator for computing the gamma */
__inline static void metad_add_varv(metad_t *metad)
{
  cmvar_add(metad->cm, metad->v);
}

/* compute the autocorrelation integrals (gamma)
 * from the variance of the bias potential (cmvar) */
__inline static void metad_getgamma_varv(metad_t *metad, double alpha0,
    const char *fn)
{
  int i, n = metad->n;

  cmvar_get(metad->cm);
  for ( i = 1; i < n; i++ ) {
    double lam = metad->lambda[i];
    if ( lam < 0.1 ) lam = 0.1;
    metad->gamma[i] = metad->cm->uvar[i]*2/alpha0/lam;
  }
  if ( fn != NULL ) savegamma(n, metad->gamma, fn);
}

/* compute the normalized transition matrix */
__inline static double *metad_normalize_tmat(metad_t *metad)
{
  int i, j, n = metad->n;
  double *mat, *col, x, y, z, dx;

  xnew(mat, n * n);
  xnew(col, n);
  for ( j = 0; j < n; j++ ) {
    for ( col[j] = 0, i = 0; i < n; i++ )
      col[j] += metad->tmat[i*n + j];
    if ( col[j] > 0 ) {
      for ( i = 0; i < n; i++ )
        mat[i*n + j] = metad->tmat[i*n + j] / col[j];
    } else {
      for ( i = 0; i < n; i++ )
        mat[i*n + j] = ( i == j );
    }
  }

  /* make the transition matrix symmetric
   * to satisfy detailed balance */
  for ( i = 0; i < n; i++ ) {
    for ( j = i + 1; j < n; j++ ) {
      x = mat[i*n + j];
      y = mat[j*n + i];
      z = (x + y)/2;
      dx = z - x;
      /* make sure the diagonal elements can afford the adjustment */
      if ( dx > 0 ) {
        if ( mat[j*n+j] < dx )
          dx = mat[j*n+j];
      } else {
        if ( mat[i*n+i] < -dx )
          dx = -mat[i*n+i];
      }
      mat[j*n + i] = mat[i*n + j] = x + dx;
    }
  }

  free(col);
  return mat;
}

/* Gnuplot command:
   plot "tmat.dat" matrix with image
   cf. http://www.gnuplotting.org/tag/matrix/ */
__inline static int metad_save_tmat(const double *mat, int n, const char *fn)
{
  int i, j;
  FILE *fp;
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < n; j++ )
      fprintf(fp, "%7.5f ", mat[i*n+j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

/* compute the autocorrelation integrals (gamma) from the transition matrix */
__inline static void metad_getgamma_tmat(metad_t *metad, double dt,
    const char *fn)
{
  int i, j, k, n = metad->n;
  double *mat, *val, *vec, *g, x, gam;

  mat = metad_normalize_tmat(metad);
  metad_save_tmat(mat, n, "tmat.dat");
  xnew(val, n);
  xnew(vec, n*n);
  xnew(g, n);
  eigsym(mat, val, vec, n);
  for ( i = 0; i < n; i++ ) {
    x = val[i];
    if ( x < 0 ) x = 0;
    x = pow(x, 1.0/dt);
    if ( x > 0.999999 ) x = 0.999999;
    g[i] = (1 + x)/(1 - x);
  }
  for ( k = 1; k < n; k++ ) {
    /* k: mode of updating matrix */
    gam = 0;
    for ( j = 1; j < n; j++ ) {
      /* j: mode of transition matrix */
      /* matrix multiplication of phi and vec */
      for ( x = 0, i = 0; i < n; i++ ) {
        x += metad->costab[k*n + i] * vec[i*n + j];
      }
      gam += g[j] * x * x;
    }
    metad->tgamma[k] = gam / n;
  }

  free(mat);
  free(val);
  free(vec);
  free(g);

  if ( fn != NULL ) savegamma(n, metad->tgamma, fn);
}


/* compute the optimal schedule */
static void metad_getalpha(metad_t *metad, double T, int gammethod,
    const char *fngamma, int sampmethod, double alpha0, double *qT,
    double qprec, int nint, const char *fnalpha)
{
  double *gamma = metad->gamma, t0 = 2/alpha0, y;

  if ( gammethod == GAMMETHOD_LOAD ) {
    if ( loadgamma(metad->n, gamma, fngamma) != 0 )
      exit(1);
  } else if ( gammethod == GAMMETHOD_TMAT ) {
    gamma = metad->tgamma;
  } else if ( gammethod == GAMMETHOD_NONE ) {
    estgamma(gamma, metad->n, sampmethod, metad->pbc, 1.0);
  }
  y = esterror_opt(T, alpha0, 0, qT, qprec,
      nint, &metad->intq, metad->n, -1, metad->pbc,
      metad->lambda, gamma, 0);
  metad->errref = y * y;
  /* save the optimal schedule to file */
  intq_save(metad->intq, 1.0, t0, 0, fnalpha);
}



/* compute the root-mean-squared error of `v`
 * note: the baseline of `v` will be shifted  */
static double metad_geterror(metad_t *metad)
{
  double err = 0, x;
  int i, n = metad->n;

  metad_trimv(metad, metad->v);
  for ( i = 0; i < n; i++ ) {
    x = metad->v[i] - metad->vref[i];
    //printf("i %d, dif %g = %g, %g\n", i, x, metad->v[i], metad->vref[i]);
    err += x * x;
  }
  //getchar();
  err /= n;
  return err;
}



#endif /* METAD_H__ */
