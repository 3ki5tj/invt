/* standalone histogram reweighting program */

#include "util.h"

#ifndef SQRT2
#define SQRT2 1.4142135623730951
#endif

const char *fndat = "pt2gaus.dat";
const char *fnhis = "pt2gaus.his";
int n; /* number of Gaussians */
double *ave, *sig, *c1, *c2, *c0, *cnt;
double xcmin, xcmax, delx;
double *hist, *htot, *lng, *lngs;
int xmin, xmax, dx, xn;

/* load the dat file */
static int getdatinfo(const char *fn)
{
  FILE *fp;
  int i, i1;
  char buf[256];

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fgets(buf, sizeof buf, fp);
  sscanf(buf, "# %d", &n);
  xnew(ave, n);
  xnew(sig, n);
  xnew(c1, n);
  xnew(c2, n);
  xnew(c0, n);
  xnew(cnt, n);
  for ( i = 0; i < n; i++ ) {
    fgets(buf, sizeof buf, fp);
    sscanf(buf, "%d%lf%lf%lf%lf%lf", &i1, &ave[i], &sig[i], &c1[i], &c2[i], &c0[i]);
    cnt[i] = 0;
  }
  xcmin = ave[0];
  delx = ave[1] - ave[0];
  xcmax = ave[n-1];
  fclose(fp);
  return 0;
}

/* get the number of histograms and number of bins, etc
 * from the tag line of the histogram file */
static int gethisinfo(const char *fn)
{
  FILE *fp;
  int i, ix, ln = 1;
  char buf[256];
  double x, y, z, h;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  fgets(buf, sizeof buf, fp);
  sscanf(buf, "# %d %d %d %d", &n, &xn, &xmin, &dx);
  xnew(hist, n * xn);
  xnew(htot, xn);
  xnew(lng, xn);
  xnew(lngs, xn);
  for ( ix = 0; ix < xn; ix++ ) htot[ix] = 0;

  for ( i = 0; i < n; i++ ) {
    fgets(buf, sizeof buf, fp);
    ln++;
    if ( buf[0] != '#' ) {
      fprintf(stderr, "cannot read histogram %d/%d from %s, line %d\n",
          i, n, fn, ln);
      fclose(fp);
      return -1;
    }
    cnt[i] = 0;
    while ( fgets(buf, sizeof buf, fp) ) {
      ln++;
      strstrip(buf);
      if (buf[0] == '\0') break;
      sscanf(buf, "%lf %lf %lf %lf", &x, &y, &z, &h);
      ix = (int)((x - xmin)/dx + 0.5);
      hist[i*xn + ix] = h;
      htot[ix] += h;
      cnt[i] += h;
    }
  }
  fclose(fp);
  return 0;
}

#define LN0 -1000000

/* ln(exp(x) + exp(y)) */
static double lnadd(double x, double y)
{
  if ( x > y ) {
    y -= x;
  } else {
    y = x - y;
    x = x - y; /* the old y */
  }
  return x + log(1 + exp(y));
}

/* non-iterative histogram reweighting */
static void reweight(void)
{
  int ix, i;
  double x, dxi, ui, lnden;

  for ( ix = 0; ix < xn; ix++ ) {
    /* compute the density of states
     * g[ix] = htot[ix] / Sum_i cnt[i] exp(-ui[ix]-c0[i]) */
    lng[ix] = lnden = LN0;
    x = xmin + ix * dx;
    if ( htot[ix] <= 0 ) continue;
    for ( i = 0; i < n; i++ ) {
      if ( cnt[i] <= 0 ) continue;
      dxi = (x - ave[i])/sig[i];
      ui = c0[i] + c1[i] * dxi + c2[i] * (dxi*dxi - 1) / SQRT2;
      //printf("ix %d, x %g, i %d, dxi %g, c1 %g, c2 %g, ui %g, lnden %g\n", ix, x, i, dxi, c1[i], c2[i], ui, lnden); // getchar();
      lnden = lnadd(lnden, log(cnt[i]) - ui);
    }
    lng[ix] = log(htot[ix]) - lnden;
    //printf("ix %d, lnden %g, htot %g, lng %g\n", ix, lnden, htot[ix], lng[ix]); // getchar();
  }
}

/* smooth array `arr` to `arrs` */
static void smooth(double *arrs, double *arr, int nb)
{
  int ix, jx, jmin, jmax;
  double sx, sy;

  for ( ix = 0; ix < xn; ix++ ) {
    if ( arr[ix] <= LN0 ) {
      arrs[ix] = LN0;
      continue;
    }
    jmin = ix - nb;
    if ( jmin < 0 ) jmin = 0;
    jmax = ix + nb + 1;
    if ( jmax > xn ) jmax = xn;
    sx = sy = 0;
    for ( jx = jmin; jx < jmax; jx++ ) {
      if ( arr[jx] <= LN0 ) continue;
      sy += arr[jx];
      sx += 1;
    }
    arrs[ix] = sy / sx;
  }
}

/* find the peak of arr[i] - x_i * bc in [ileft, iright) */
static int findpeak(const double *arr, double bc, int ileft, int iright, double *ym)
{
  double x, y;
  int i, imax = ileft;

  *ym = -1e100;
  for ( i = ileft; i < iright; i++ ) {
    if ( arr[i] <= LN0 ) continue;
    x = xmin + i * dx;
    y = arr[i] - x * bc;
    //printf("i %d, x %g, y %g, ymax %g\n", i, x, y, ymax); getchar();
    if ( y > *ym ) {
      *ym = y;
      imax = i;
    }
  }
  return imax;
}

/* find the critical inverse temperature */
static double seekcrit(const double *arr, int *x1, int *x2, double *shift)
{
  double bc = 0, y1, y2, y, lns, lns2;
  int i, ix, ix1, ix2, im, t, x;

  // since c1/sigma is roughly beta, we will estimate the critical temperature
  // at a point where c1/sigma rises quickest
  double dbmax = 0, b0, b1, db;
  for ( i = 0; i < n - 1; i++ ) {
    b0 = c1[i]/sig[i];
    b1 = c1[i+1]/sig[i+1];
    db = (b1 - b0) / (ave[i+1] - ave[i]);
    if ( db > dbmax ) {
      dbmax = db;
      bc = b0;
      im = (int) ((ave[i] - xmin) / dx);
    }
  }
  printf("estimated bc %g, mid %d\n", bc, im);

  /* adjust bc such that the heights are equal */
  for ( t = 0; t < 10; t++ ) {
    ix1 = findpeak(arr, bc, 0, im, &y1);
    ix2 = findpeak(arr, bc, im, xn, &y2);
    printf("%d: bc %g, x %d(%g) %d(%g) | %d(%d) %g\n", t, bc, xmin + ix1*dx, y1, xmin + ix2*dx, y2, xmin + im*dx, im, fabs(y2-y1));
    if ( fabs(y2 - y1) < 1e-10 ) break;
    bc += (y2 - y1) / ((ix2 - ix1)*dx);
  }
  *x1 = xmin + ix1*dx;
  *x2 = xmin + ix2*dx;
  printf("bc %.8f, T %.8f, x %d %d\n", bc, 1/bc, xmin + ix1*dx, xmin + ix2*dx);

  /* normalize lng at the critical point */
  lns = LN0;
  for ( ix = 0; ix < xn; ix++ ) {
    if ( arr[ix] <= LN0 ) continue;
    x = xmin + ix * dx;
    y = arr[ix] - x * bc;
    lns = lnadd(lns, y);
  }
  lns += log(dx);

  /* normalize c0 at the critical temperature */
  lns2 = LN0;
  for ( i = 0; i < n; i++ )
    lns2 = lnadd(lns2, c0[i] - c2[i]/SQRT2 - bc * ave[i]);
  lns2 += log(delx);
  printf("lns for normalization %g (lng), %g (c0hat)\n", lns, lns2);
  *shift = lns;
  return bc;
}

static int save(const char *fn, double bc, int x1, int x2)
{
  int ix, i1, i0;
  FILE *fp;
  double x, beta, y, ys;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d %d %g %d %d\n", xn, xmin, dx, bc, x1, x2);
  for ( ix = 0; ix < xn; ix++ ) {
    if ( htot[ix] <= 0 ) continue;
    x = xmin + ix * dx;
    i0 = ix - 1;
    if ( i0 < 0 || htot[i0] <= 0 ) i0 = ix;
    i1 = ix + 1;
    if ( i1 >= xn || htot[i1] <= 0 ) i1 = ix;
    if ( i0 < i1 ) beta = (lngs[i1] - lngs[i0])/(dx*(i1 - i0));
    else beta = 0;
    y = exp(lng[ix] - bc*x);
    ys = exp(lngs[ix] - bc*x);
    fprintf(fp, "%g %g %g %g %g %g %g\n",
        x, lng[ix], beta, y, htot[ix], lngs[ix], ys);
  }
  fclose(fp);
  fprintf(stderr, "saved lng to %s\n", fn);
  return 0;
}


int main(int argc, char **argv)
{
  double bc, lns;
  int ix, x1, x2, nb = 25;

  if ( argc > 1 ) fndat = argv[1];
  if ( argc > 2 ) fnhis = argv[2];
  if ( argc > 3 ) nb = atoi( argv[3] );

  if ( getdatinfo(fndat) != 0 ) return -1;
  if ( gethisinfo(fnhis) != 0 ) return -1;

  /* non-iterative histogram reweighting */
  reweight();

  /* smoothing lng to lngs */
  smooth(lngs, lng, nb);

  /* looking for the critical temperature */
  bc = seekcrit(lng, &x1, &x2, &lns);
  for ( ix = 0; ix < xn; ix++ ) {
    lng[ix] -= lns;
    lngs[ix] -= lns;
  }
  save("lng.dat", bc, x1, x2);
  return 0;
}
