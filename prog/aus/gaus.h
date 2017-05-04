#ifndef GAUS_H__
#define GAUS_H__


/* Generalized ensemble of adaptive umbrella sampling
 * */


#include "mmwl.h"

typedef struct {
  int n;
  double *ave;
  double *sig;
  double *beta1;
  double *beta2;
  double *lnz; /* partition function */
  double alpha; /* updating magnitude */
  double hflatness;
  double t, t0;
  int invt; /* using 1/t schedule */
  mmwl_t *mmwl;
  double *cnt;
  double *acc;
  int xmin, xmax, dx;
  int xn;
  double *hist, *htot;
} gaus_t;



__inline static gaus_t *gaus_open(double xcmin, double xcmax,
    int n, double sig, double beta,
    int xmin, int xmax, int dx)
{
  gaus_t *gaus;
  int i, xn;
  double delx;

  xnew(gaus, 1);
  gaus->n = n;
  xnew(gaus->ave, n);
  xnew(gaus->sig, n);
  xnew(gaus->beta1, n);
  xnew(gaus->beta2, n);
  xnew(gaus->lnz, n);
  xnew(gaus->mmwl, n);
  xnew(gaus->cnt, n);
  xnew(gaus->acc, n);
  delx = (xcmax - xcmin) / (n - 1);
  for ( i = 0; i < n; i++ ) {
    gaus->ave[i] = xcmin + i * delx;
    gaus->sig[i] = sig;
    gaus->beta1[i] = beta;
    gaus->beta2[i] = 0;
    gaus->lnz[i] = 0;
    mmwl_init(gaus->mmwl + i, 1e-3);
    gaus->cnt[i] = 0;
    gaus->acc[i] = 0;
  }
  gaus->xmin = xmin;
  gaus->dx = dx;
  gaus->xn = xn = (xmax - xmin) / gaus->dx + 1;
  gaus->xmax = gaus->xmin + gaus->dx * xn;
  xnew(gaus->hist, n * xn);
  xnew(gaus->htot, xn);
  for ( i = 0; i < n * xn; i++ ) {
    gaus->hist[i] = 0;
  }
  for ( i = 0; i < xn; i++ ) {
    gaus->htot[i] = 0;
  }
  gaus->alpha = 1;
  gaus->t = 0;
  gaus->t0 = 1;
  gaus->invt = 0;
  fprintf(stderr, "%d states Ec %g ... %g\n",
      n, gaus->ave[0], gaus->ave[n-1]);

  return gaus;
}



__inline static void gaus_close(gaus_t *gaus)
{
  free(gaus->ave);
  free(gaus->sig);
  free(gaus->beta1);
  free(gaus->beta2);
  free(gaus->lnz);
  free(gaus->mmwl);
  free(gaus->cnt);
  free(gaus->acc);
  free(gaus->hist);
  free(gaus->htot);
  free(gaus);
}



/* low level histogram */
static void savehistlow(int id, double *h, int n,
    double xmin, double dx, double ave, double sig, FILE *fp)
{
  int i;
  double tot = 0, f = 0, x, x2, f0;

  for ( i = 0; i < n; i++ )
    tot += h[i];
  if ( tot <= 0 ) return;
  f0 = 1./sqrt(2*M_PI)/sig;
  for ( i = 0; i < n; i++ ) {
    if ( h[i] <= 0 ) continue;
    x = xmin + i*dx;
    if ( sig > 0 ) {
      x2 = (x - ave) / sig;
      f = exp(-0.5 * x2 * x2) * f0;
    }
    fprintf(fp, "%g %g %g %g %d\n",
        x, h[i]/tot/dx, f, h[i], id);
  }
  fprintf(fp, "\n");
}


static void gaus_trimv(gaus_t *gaus, double *v)
{
  int i, n = gaus->n;
  double v0 = 0;

  for ( i = 0; i < n; i++ ) v0 += v[i];
  v0 /= n;
  for ( i = 0; i < n; i++ ) v[i] -= v0;
}



/* save histogram file */
__inline static int gaus_savehist(gaus_t *gaus, const char *fn)
{
  FILE *fp;
  int i, j, n = gaus->n, xn = gaus->xn;
  double xmin = gaus->xmin, dx = gaus->dx;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  gaus_trimv(gaus, gaus->lnz);
  for ( j = 0; j < xn; j++ )
    gaus->htot[j] = 0;
  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "# %g %g %g %g %g %g %d %g %g %g\n", gaus->ave[i], gaus->sig[i],
        gaus->beta1[i], gaus->beta2[i], gaus->lnz[i], gaus->cnt[i],
        gaus->mmwl[i].invt, mmwl_getalpha(gaus->mmwl+i), gaus->mmwl[i].fl[1], gaus->mmwl[i].fl[2]);
    savehistlow(i, gaus->hist + i * xn,
        xn, xmin, dx, gaus->ave[i], gaus->sig[i], fp);
    for ( j = 0; j < xn; j++ )
      gaus->htot[j] += gaus->hist[i*xn + j];
  }
  savehistlow(-1, gaus->htot, xn, xmin, dx, 0, 0, fp);
  fclose(fp);
  return 0;
}



/* update the parameters of the bias potential */
__inline static void gaus_update(gaus_t *gaus, int i, double x, double t)
{
  double ave = gaus->ave[i], sig = gaus->sig[i], y1, y2, amp;

  y1 = (x - ave) / sig;
  y2 = y1 * y1 - 1;
  amp = mmwl_getalpha(gaus->mmwl + i);
  gaus->beta1[i] += y1 / sig * amp;
  gaus->beta2[i] += y2 / (sig * sig) * amp;
  mmwl_add(gaus->mmwl + i, y1, y2);
  if ( fmod(t, 100) < 0.1 && mmwl_check(gaus->mmwl + i, 1, 0.05, 0.5)) {
    printf("t %g, id %d, new updating magnitude %g, fl %g, %g, c %g, %g, invt %d\n",
      t, i, gaus->mmwl[i].alpha, gaus->mmwl[i].fl[1], gaus->mmwl[i].fl[2], gaus->beta1[i], gaus->beta2[i], gaus->mmwl[i].invt);
  }
  /* disable update during update stage */
  if ( gaus->mmwl[i].invt ) {
    double alpha = gaus->invt ? gaus->n / (gaus->t + gaus->t0) : gaus->alpha;
    gaus->t += 1;
    gaus->cnt[i] += 1;
    gaus->lnz[i] += alpha;
  }
}



/* update the histogram */
__inline static void gaus_add(gaus_t *gaus, int i, int x, int acc)
{
  int j;
  gaus->acc[i] += acc;
  j = (x - gaus->xmin) / gaus->dx;
  gaus->hist[i * gaus->xn + j] += 1;
}



/* compute the histogram flatness (Wang-Landau version) */
__inline static double gaus_hflatness_wl(gaus_t *gaus)
{
  int i, n = gaus->n;
  double hmin, hmax, hi;

  hmin = hmax = gaus->cnt[0];
  for ( i = 1; i < n; i++ ) {
    hi = gaus->cnt[i];
    if ( hi < hmin ) {
      hmin = hi;
    } else if ( hi > hmax ) {
      hmax = hi;
    }
  }

  if ( hmax <= 0 ) return 100.0;
  return 2 * (hmax - hmin) / (hmax + hmin + 1e-16);
}


__inline static double gaus_hflatness(gaus_t *gaus)
{
  int i, n = gaus->n;
  double tot = 0, h, dx, s = 0;

  for ( i = 0; i < n; i++ ) tot += gaus->cnt[i];
  if ( tot <= 0 ) return 100.0;
  for ( i = 0; i < n; i++ ) {
    h = gaus->cnt[i] / tot;
    dx = h * n - 1;
    s += dx * dx;
  }
  return s / n;
}


__inline static int gaus_wlcheck(gaus_t *gaus,
    double fl, double magred, double t)
{
  int i, n = gaus->n;
  gaus->hflatness = gaus_hflatness_wl(gaus);
  if ( gaus->hflatness > fl || gaus->invt ) {
    return 0;
  }
  gaus->alpha *= magred;
  if ( gaus->alpha < gaus->n/(gaus->t + gaus->t0) ) {
    gaus->invt = 1;
    gaus->t0 = gaus->t;
  }
  fprintf(stderr, "alpha %g, %g, flatness %g, invt %d\n",
      gaus->alpha, gaus->n/gaus->t, gaus->hflatness, gaus->invt);
  gaus->t = 0;
  for ( i = 0; i < n; i++ ) {
    gaus->cnt[i] = 0;
  }
  for ( i = 0; i < n * gaus->xn; i++ )
    gaus->hist[i] = 0;
  return 1;
}



__inline static int gaus_bmove(gaus_t *gaus, double x, int *id)
{
  int acc = 0;
  int jd = (rand01() < 0.5) ? *id - 1 : *id + 1;
  double xi, xj, vi, vj, dv;

  if ( jd < 0 || jd >= gaus->n ) return 0;
  /* disable transition during WL stage */
  if ( !gaus->mmwl[*id].invt ) return 0;
  /* compute the acceptance probability */
  xi = x - gaus->ave[*id];
  xj = x - gaus->ave[jd];
  vi = gaus->beta1[*id] * xi + 0.5 * gaus->beta2[*id] * xi * xi + gaus->lnz[*id];
  vj = gaus->beta1[jd] * xj + 0.5 * gaus->beta2[jd] * xj * xj + gaus->lnz[jd];
  dv = vj - vi;
  acc = ( dv <= 0 || rand01() < exp(-dv) );
  if ( acc ) *id = jd;
  return acc;
}


#endif /* GAUS_H__ */
