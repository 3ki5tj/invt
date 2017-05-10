#include "util.h"
//#include "corr.h"
#include "cosmodes.h"
#include "eig.h"
#include "intq.h"



/* metadynamics for integer */
typedef struct {
  int xmin, xmax, xdel;
  int n;
  double *v; /* bias potential */
  double *h; /* histogram */
  double *vref; /* reference bias potential */
  double a; /* updating magnitude */
  int pbc; /* periodic boundary condition */
  int winn;
  double *win;
  double *lambda; /* eigenvalues of the updating marix */
  double *gamma; /* correlation integrals */
  double *gamma_tmat; /* correlation integrals from the transition matrix */
  double *costab; /* cosine transform coefficients */
  double hflatness;
  double *tmat; /* n x n transition matrix */
  //corr_t *corr;
  double *vtmp;
  double *hmod;
  intq_t *intq;
  double errref;
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
  //  /* save the window kernel */
  //  savewin(metad->win, metad->winn, m->fnwin);
  //  /* save the n x n updating matrix */
  //  savewinmat(metad->win, metad->winn, n, pbc, m->fnwinmat);
}

static metad_t *metad_open(int xmin, int xmax, int xdel,
    int pbc, double gaussig, int okmax, const double *win, int winn)
{
  metad_t *metad;
  int i, n;

  xnew(metad, 1);
  metad->xmin = xmin;
  metad->xmax = xmax;
  metad->xdel = xdel;
  metad->n = n = (xmax - xmin) / xdel + 1;
  metad->xmax = metad->xmin + n * metad->xdel;
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
  metad_prepwin(metad, gaussig, okmax, win, winn);
  xnew(metad->gamma, n);
  xnew(metad->gamma_tmat, n);
  metad->costab = mkcostab(n, metad->pbc);
  xnew(metad->tmat, n * n);
  xnew(metad->vtmp, n);
  xnew(metad->hmod, n);
  metad->intq = NULL;
  metad->errref = 0;
  return metad;
}



static void metad_close(metad_t *metad)
{
  free(metad->v);
  free(metad->h);
  free(metad->win);
  free(metad->lambda);
  free(metad->gamma);
  free(metad->gamma_tmat);
  free(metad->costab);
  free(metad->tmat);
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
  if ( x < metad->xmin || x >= metad->xmax )
    return -1;
  return (x - metad->xmin) / metad->xdel;
}



/* decide if a transition from xold (index iold)
 * to xnew (index iold) is to be accepted */
static int metad_acc(metad_t *metad, int iold, int xnew,
    int *inew)
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


/* compute the histogram flatness (Wang-Landau version) */
__inline static double metad_hflatness_wl(metad_t *metad)
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



/* compute the histogram flatness */
__inline static double metad_hflatness(metad_t *metad)
{
  int i, k, n = metad->n;
  double x, fl = 0, tot = 0;
  const double lamcut = 0.1;

  for ( i = 0; i < n; i++ ) tot += metad->h[i];
  getcosmodes(metad->h, n, metad->vtmp, metad->costab);
  /* truncate modes with lambda */
  for ( k = 1; k < n; k++ ) {
    if ( metad->lambda[k] < lamcut ) break;
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

  /* compute the histogram flatness */
  metad->hflatness = metad_hflatness(metad);
  //printf("flatness %g, %g\n", metad->hflatness, metad_hflatness_wl(metad));
  //getchar();
  /* return if the histogram not flatness enough */
  if ( metad->hflatness > fl ) return 0;

  /* reduce the updating magnitude and clear the histogram */
  metad->a *= magred;
  for ( i = 0; i < metad->n; i++ ) metad->h[i] = 0;
  fprintf(stderr, "changing the updating magnitude to %g\n", metad->a);
  return 1;
}



static void metad_updatev_wl(metad_t *metad, int i)
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



static int metad_save(metad_t *metad, const char *fn)
{
  FILE *fp;
  int i;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  metad_trimv(metad, metad->v);
  fprintf(fp, "# %d %d %d %d %g\n",
      metad->n, metad->xmin, metad->xmax, metad->xdel, metad->a);
  for ( i = 0; i < metad->n; i++ ) {
    fprintf(fp, "%d %g %g %g %g\n",
        metad->xmin + i * metad->xdel,
        metad->v[i], metad->h[i], metad->vref[i], metad->hmod[i]);
  }
  fclose(fp);
  return 0;
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

  {
    // Gnuplot command:
    //   plot "tmat.dat" matrix with image
    // cf. http://www.gnuplotting.org/tag/matrix/
    FILE *fp = fopen("tmat.dat", "w");
    for ( i = 0; i < n; i++ ) {
      for ( j = 0; j < n; j++ )
        fprintf(fp, "%5.3f ", mat[i*n+j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  free(col);
  return mat;
}

/* compute the gamma values from the transition matrix */
static void metad_getgamma_tmat(metad_t *metad, double *gam, double dt)
{
  int i, j, k, n = metad->n;
  double *mat, *val, *vec, *g, x;

  mat = metad_normalize_tmat(metad);
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
    gam[k] = 0;
    for ( j = 1; j < n; j++ ) {
      /* j: mode of transition matrix */
      /* matrix multiplication of phi and vec */
      for ( x = 0, i = 0; i < n; i++ ) {
        x += metad->costab[k*n + i] * vec[i*n + j];
      }
      gam[k] += g[j] * x * x;
    }
    gam[k] /= n;
  }

  free(mat);
  free(val);
  free(vec);
  free(g);
}


static void metad_getalpha(metad_t *metad, invtpar_t *m, double T)
{
  m->t0 = 2 / m->alpha0;
  metad->errref = esterror_opt(T, m->alpha0, 0, &m->qT, m->qprec,
      m->alpha_nint, &metad->intq, m->n, m->kcutoff, m->pbc,
      metad->lambda, metad->gamma, m->verbose);
  /* save the optimal schedule to file */
  intq_save(metad->intq, 1.0, m->t0,
      m->alpha_resample, m->fnalpha);
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
    err += x * x;
  }
  err /= n;
  return err;
}




