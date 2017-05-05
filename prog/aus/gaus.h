#ifndef GAUS_H__
#define GAUS_H__


/* Generalized ensemble of adaptive umbrella sampling
 * */


#include "mmwl.h"
#include "cosmodes.h"


/* values for lnzmethod */
enum { LNZ_WL, LNZ_AVE };


typedef struct {
  int n;
  double *ave;
  double *sig;
  double *c1, *c2;
  double *lnz; /* partition function */
  int lnzmethod; /* WL or average */
  double alphawl; /* updating magnitude */
  double hflatness;
  double t, t0;
  int invt; /* using 1/t schedule */
  double alpha0; /* initial updating magnitude */
  double alphamm0; /* initial updating magnitude for moments */
  mmwl_t *mmwl;
  double *cnt;
  double *acc;
  int xmin, xmax, dx;
  int xn;
  double *hist, *htot;
  double *hfl; /* fluctuation of histogram modes */
  double *costab;
} gaus_t;



__inline static gaus_t *gaus_open(double xcmin, double xcmax,
    int n, double sig, int lnzmethod, double c1, double alpha0,
    int xmin, int xmax, int dx, int pbc)
{
  gaus_t *gaus;
  int i, xn;
  double delx;

  xnew(gaus, 1);
  gaus->n = n;
  gaus->lnzmethod = lnzmethod;
  gaus->alpha0 = alpha0;
  gaus->alphamm0 = alpha0;
  xnew(gaus->ave, n);
  xnew(gaus->sig, n);
  xnew(gaus->c1, n);
  xnew(gaus->c2, n);
  xnew(gaus->lnz, n);
  xnew(gaus->mmwl, n);
  xnew(gaus->cnt, n);
  xnew(gaus->acc, n);
  xnew(gaus->hfl, n);
  delx = (xcmax - xcmin) / (n - 1);
  for ( i = 0; i < n; i++ ) {
    gaus->ave[i] = xcmin + i * delx;
    gaus->sig[i] = sig;
    gaus->c1[i] = c1;
    gaus->c2[i] = 0;
    gaus->lnz[i] = 0;
    mmwl_init(gaus->mmwl + i, gaus->alphamm0);
    gaus->cnt[i] = 0;
    gaus->acc[i] = 0;
    gaus->hfl[i] = 0;
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
  gaus->alphawl = gaus->alpha0;
  gaus->t = 0;
  gaus->t0 = 1;
  gaus->invt = 0;
  gaus->costab = mkcostab(n, pbc);
  fprintf(stderr, "%d states Ec %g ... %g\n",
      n, gaus->ave[0], gaus->ave[n-1]);

  return gaus;
}



__inline static void gaus_close(gaus_t *gaus)
{
  free(gaus->ave);
  free(gaus->sig);
  free(gaus->c1);
  free(gaus->c2);
  free(gaus->lnz);
  free(gaus->mmwl);
  free(gaus->cnt);
  free(gaus->acc);
  free(gaus->hist);
  free(gaus->htot);
  free(gaus->costab);
  free(gaus);
}



/* save temperature dependent data */
__inline static void gaus_save(gaus_t *gaus, const char *fn)
{
  int i, n = gaus->n;
  FILE *fp;
  mmwl_t *mm;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fp = stderr;
  }
  fprintf(fp, "# %d %d %d %g %g\n", n, gaus->lnzmethod, gaus->invt, gaus->t, gaus->t0);
  for ( i = 0; i < n; i++ ) {
    mm = gaus->mmwl + i;
    fprintf(fp, "%4d %12.5f %12.5f %12.5f %12.5f %12.5f ", i,
        gaus->ave[i], gaus->sig[i], gaus->c1[i], gaus->c2[i], gaus->lnz[i]);
    fprintf(fp, "%10.0f %10.7f %10.7f %d\n", mm->mm[0], mm->fl[1], mm->fl[2], mm->invt);
  }
  if ( fp != stderr ) fclose(fp);
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
    } else {
      f = 0;
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
        gaus->c1[i], gaus->c2[i], gaus->lnz[i], gaus->cnt[i],
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



/* retrieve the updating magnitude */
__inline static double gaus_getalpha(const gaus_t *gaus, double *alphamm)
{
  double alpha = gaus->invt ? gaus->n / (gaus->t + gaus->t0) : gaus->alphawl;
  *alphamm = alpha;
  if ( *alphamm > gaus->alphamm0 ) *alphamm = gaus->alphamm0;
  return alpha;
}



/* update the parameters of the bias potential */
__inline static void gaus_add(gaus_t *gaus, int i, int x, int acc)
{
  double ave = gaus->ave[i], sig = gaus->sig[i], y1, y2, alphamm, alpha;
  int j;

  alpha = gaus_getalpha(gaus, &alphamm);

  y1 = (x - ave) / sig;
  y2 = y1 * y1 - 1;
  gaus->c1[i] += y1 * alphamm;
  gaus->c2[i] += y2 * alphamm;
  mmwl_add(gaus->mmwl + i, y1, y2);

  gaus->t += 1;
  gaus->cnt[i] += 1;
  gaus->lnz[i] += alpha;
  gaus->acc[i] += acc;
  j = (x - gaus->xmin) / gaus->dx;
  gaus->hist[i * gaus->xn + j] += 1;
}



/* compute the histogram flatness (Wang-Landau version) */
__inline static double gaus_hflatness_wl(gaus_t *gaus)
{
  int i, n = gaus->n;
  double hmin = 1e30, hmax = 0, hi;

  for ( i = 0; i < n; i++ ) {
    hi = gaus->mmwl[i].mm[i];
    if ( hi < hmin ) {
      hmin = hi;
    } else if ( hi > hmax ) {
      hmax = hi;
    }
  }

  if ( hmax <= 0 ) return 99.0;
  return 2 * (hmax - hmin) / (hmax + hmin + 1e-16);
}


__inline static double gaus_hflatness(gaus_t *gaus)
{
  int i, n = gaus->n;
  double tot = 0, h, dx, s = 0;

  for ( i = 0; i < n; i++ ) tot += gaus->mmwl[i].mm[0];
  if ( tot <= 0 ) return 100.0;
  for ( i = 0; i < n; i++ ) {
    h = gaus->mmwl[i].mm[0] / tot;
    dx = h * n - 1;
    s += dx * dx;
  }
  return s / n;
}


/* calculate fluctuation of histogram modes */
__inline static double gaus_calcfl(gaus_t *gaus)
{
  int i, k, n = gaus->n;
  double tot = 0, s, fl = 0;

  for ( i = 0; i < n; i++ )
    tot += gaus->mmwl[i].mm[0];
  if ( tot <= 0 ) return 99.0;

  for ( k = 1; k < n; k++ ) {
    s = 0;
    for ( i = 0; i < n; i++ ) {
      s += gaus->costab[k*n + i] * gaus->mmwl[i].mm[0] / tot;
    }
    gaus->hfl[k] = (s /= n);
    s = fabs(s);
    if ( s > fl ) fl = s;
  }
  return fl;
}


/* switch a stage */
__inline static void gaus_switch(gaus_t *gaus, double magred, int extended)
{
  int i, n = gaus->n;

  gaus->alphawl *= magred;
  if ( gaus->alphawl < gaus->n/(gaus->t + gaus->t0) ) {
    gaus->invt = 1;
    gaus->t0 = gaus->t;
  }
  fprintf(stderr, "alpha %g, %g, flatness %g, invt %d\n",
      gaus->alphawl, gaus->n/gaus->t, gaus->hflatness, gaus->invt);
  gaus->t = 0;
  for ( i = 0; i < n; i++ ) {
    //gaus->cnt[i] = 0;
    if ( extended ) {
      mmwl_init(gaus->mmwl + i, gaus->alphawl);
      gaus->mmwl[i].invt = gaus->invt;
      gaus->mmwl[i].t0 = gaus->t0;
    }
  }
  for ( i = 0; i < n * gaus->xn; i++ )
    gaus->hist[i] = 0;
}



__inline static int gaus_wlcheck(gaus_t *gaus,
    double fl, double magred)
{
  gaus_calcfl(gaus);
  if ( !gaus->invt && gaus->hflatness < fl ) {
    gaus_switch(gaus, magred, 0);
    return 1;
  } else {
    return 0;
  }
}



/* extensive check of histogram fluctuation */
__inline static int gaus_wlcheckx(gaus_t *gaus,
    double fl, double magred)
{
  int i, n = gaus->n;
  double f, f2;

  f = gaus_calcfl(gaus);
  for ( i = 0; i < n; i++ ) {
    f2 = mmwl_calcfl(gaus->mmwl + i, 1);
    if ( f2 > f ) f = f2;
  }
  gaus->hflatness = f;
  if ( !gaus->invt && f < fl ) {
    gaus_switch(gaus, magred, 1);
    return 1;
  } else {
    return 0;
  }
}




__inline static int gaus_bmove(gaus_t *gaus, double x, int *id)
{
  int acc = 0;
  int jd = (rand01() < 0.5) ? *id - 1 : *id + 1;
  double xi, xj, vi, vj, dv;

  if ( jd < 0 || jd >= gaus->n ) return 0;
  if ( gaus->lnzmethod == LNZ_AVE ) {
    /* disable transition during WL stage */
    if ( !gaus->mmwl[*id].invt ) return 0;
    /* approximate change of lnz */
    dv = 0.5 * (gaus->c1[*id] + gaus->c1[jd])
             * (gaus->ave[jd] - gaus->ave[*id]);
  } else {
    dv = gaus->lnz[jd] - gaus->lnz[*id];
  }
  /* compute the acceptance probability */
  xi = (x - gaus->ave[*id]) / gaus->sig[*id];
  xj = (x - gaus->ave[ jd]) / gaus->sig[ jd];
  vi = gaus->c1[*id] * xi + gaus->c2[*id] * 0.5 * xi * xi;
  vj = gaus->c1[ jd] * xj + gaus->c2[ jd] * 0.5 * xj * xj;
  dv += vj - vi;
  acc = ( dv <= 0 || rand01() < exp(-dv) );
  if ( acc ) *id = jd;
  return acc;
}


#endif /* GAUS_H__ */
