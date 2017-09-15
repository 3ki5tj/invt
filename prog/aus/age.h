#ifndef AGE_H__
#define AGE_H__


/* Multiple Gaussian ensembles */


#include "mtrand.h"


#ifndef SQRT2
#define SQRT2 1.4142135623730951
#endif

typedef struct {
  int n;
  double *ave;
  double *sig;
  double *c0, *c1, *c2;
  double (*mm)[3];
  double alphawl; /* updating magnitude */
  double hfluc, flfr[3];
  double t, t0;
  int invt; /* using 1/t schedule */
  double alpha0; /* initial updating magnitude */
  double *cnt;
  double *acc;
  int xmin, xmax, dx;
  int xn;
  double *hist, *htot;
  double *hfl; /* fluctuation of histogram modes */
} age_t;



__inline static age_t *age_open(double xcmin, double xcmax,
    double delx, double sig, double c1, double alpha0,
    int xmin, int xmax, int dx)
{
  age_t *age;
  int i, xn, n;

  xnew(age, 1);
  n = (int)((xcmax - xcmin) / delx) + 1;
  age->n = n;
  age->alpha0 = alpha0;
  xnew(age->ave, n);
  xnew(age->sig, n);
  xnew(age->c0, n);
  xnew(age->c1, n);
  xnew(age->c2, n);
  xnew(age->mm, n);
  xnew(age->cnt, n);
  xnew(age->acc, n);
  xnew(age->hfl, n);
  for ( i = 0; i < n; i++ ) {
    age->ave[i] = xcmin + i * delx;
    age->sig[i] = sig;
    age->c0[i] = (i - (n - 1)*0.5) * c1 / sig * delx;
    age->c1[i] = c1;
    age->c2[i] = SQRT2*0.5;
    age->mm[i][0] = 0;
    age->mm[i][1] = 0;
    age->mm[i][2] = 0;
    age->cnt[i] = 0;
    age->acc[i] = 0;
    age->hfl[i] = 0;
  }
  age->xmin = xmin;
  age->dx = dx;
  age->xn = xn = (xmax - xmin) / age->dx + 1;
  age->xmax = age->xmin + age->dx * xn;
  xnew(age->hist, n * xn);
  xnew(age->htot, xn);
  for ( i = 0; i < n * xn; i++ ) {
    age->hist[i] = 0;
  }
  for ( i = 0; i < xn; i++ ) {
    age->htot[i] = 0;
  }
  age->alphawl = age->alpha0;
  age->t = 0;
  age->t0 = 1;
  age->invt = 0;
  fprintf(stderr, "%d states Ec %g, %g, ... %g\n",
      n, age->ave[0], age->ave[1], age->ave[n-1]);

  return age;
}



__inline static void age_close(age_t *age)
{
  free(age->ave);
  free(age->sig);
  free(age->c0);
  free(age->c1);
  free(age->c2);
  free(age->mm);
  free(age->cnt);
  free(age->acc);
  free(age->hist);
  free(age->htot);
  free(age);
}



static void age_trimv(age_t *age, double *v)
{
  int i, n = age->n;
  double v0 = 0;

  for ( i = 0; i < n; i++ ) v0 += v[i];
  v0 /= n;
  for ( i = 0; i < n; i++ ) v[i] -= v0;
}



/* retrieve the updating magnitude */
__inline static double age_getalpha(const age_t *age)
{
  return age->invt ? age->n / (age->t + age->t0) : age->alphawl;
}



/* save temperature dependent data */
__inline static void age_save(age_t *age, const char *fn)
{
  int i, n = age->n;
  FILE *fp;
  double alpha, tot = 0, acc = 0, mcnt;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fp = stderr;
  }
  alpha = age_getalpha(age);
  /* compute the overall acceptance ratio */
  for ( i = 0; i < n; i++ ) {
    acc += age->acc[i];
    tot += age->cnt[i];
  }
  fprintf(fp, "# %d 0 %d %g %g %14.7e %14.7e %.4f\n", n,
      age->invt, age->t, age->t0, alpha, alpha, 1.0*acc/tot);
  for ( i = 0; i < n; i++ ) {
    mcnt = age->mm[i][0] + 1e-14;
    fprintf(fp, "%4d %12.5f %12.5f %12.5f %12.5f %12.5f %10.0f ", i,
        age->ave[i], age->sig[i], age->c1[i], age->c2[i],
        age->c0[i], age->cnt[i]);
    fprintf(fp, "%10.0f %10.7f %10.7f %.4f\n",
        mcnt, age->mm[i][1]/mcnt, age->mm[i][2]/mcnt,
        age->acc[i]/(age->cnt[i]+1e-12));
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
    if ( sig > 0 ) { /* compute the reference distribution */
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


/* save histogram file */
__inline static int age_savehist(age_t *age, const char *fn)
{
  FILE *fp;
  int i, j, n = age->n, xn = age->xn;
  int xmin = age->xmin, dx = age->dx;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  /* print the header */
  fprintf(fp, "# %d %d %d %d\n", n, xn, xmin, dx);
  age_trimv(age, age->c0);
  for ( j = 0; j < xn; j++ )
    age->htot[j] = 0;
  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "# %g %g %g %g %g %g\n", age->ave[i], age->sig[i],
        age->c1[i], age->c2[i], age->c0[i], age->cnt[i]);
    savehistlow(i, age->hist + i * xn,
        xn, xmin, dx, age->ave[i], age->sig[i], fp);
    for ( j = 0; j < xn; j++ )
      age->htot[j] += age->hist[i*xn + j];
  }
  savehistlow(-1, age->htot, xn, xmin, dx, 0, 0, fp);
  fclose(fp);
  return 0;
}

/*
__inline static int age_loadhist(age_t *age, const char *fn)
{
  FILE *fp;
  char buf[256];
  int i, ix, n = age->n;
  double x, y, z, h;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    return -1;
  }
  fgets(buf, sizeof buf, fp);
  for ( i = 0; i < n; i++ ) {
    fgets(buf, sizeof buf, fp);
    if ( buf[0] != '#' ) {
      fprintf(stderr, "cannot read histogram %d/%d from %s\n",
          i, n, fn);
      fclose(fp);
      return -1;
    }
    age->cnt[i] = 0;
    while ( fgets(buf, sizeof buf, fp) ) {
      strstrip(buf);
      if (buf[0] == '\0') break;
      sscanf(buf, "%lf %lf %lf %lf", &x, &y, &z, &h);
      ix = (int)((x - age->xmin)/age->dx + 0.5);
      age->hist[i*age->xn + ix] = h;
      age->htot[ix] += h;
      age->cnt[i] += h;
    }
  }
  fclose(fp);
  return 0;
}
*/




/* update the parameters of the bias potential */
__inline static void age_add(age_t *age, int i, int x, int acc)
{
  double ave = age->ave[i], sig = age->sig[i], y1, y2, alpha;
  int j;

  alpha = age_getalpha(age);

  y1 = (x - ave) / sig;
  y2 = (y1 * y1 - 1) / SQRT2;
  age->c0[i] += alpha;
  age->c1[i] += y1 * alpha;
  age->c2[i] += y2 * alpha;
  age->mm[i][0] += 1;
  age->mm[i][1] += y1;
  age->mm[i][2] += y2;

  age->t += 1;
  age->cnt[i] += 1;
  age->acc[i] += acc;
  j = (x - age->xmin) / age->dx;
  age->hist[i * age->xn + j] += 1;
}



/* calculate the fluctuation of histogram modes (from the variance) */
__inline static double age_calcfl_rms(age_t *age)
{
  int i, k, n = age->n;
  double tot = 0, fl[3] = {0, 0, 0}, flsum = 0, s;

  for ( i = 0; i < n; i++ ) {
    for ( k = 0; k < 3; k++ )
      fl[k] += age->mm[i][k] * age->mm[i][k];
    tot += age->mm[i][0];
  }
  /* normalize the fluctuations */
  s = n/(tot*tot);
  for ( k = 0; k < 3; k++ ) {
    fl[k] *= s;
    if ( k == 0 ) fl[0] -= 1;
    flsum += fl[k];
  }
  for ( k = 0; k < 3; k++ )
    age->flfr[k] = fl[k]/flsum;
  return sqrt(flsum);
}


/* switch a stage */
__inline static void age_switch(age_t *age, double magred)
{
  int i, k, n = age->n;

  age->alphawl *= magred;
  if ( age->alphawl < age->n/(age->t + age->t0) ) {
    age->invt = 1;
    age->t0 = age->t;
  }
  fprintf(stderr, "alpha %g, %g, t %g, fluc %g, invt %d\n",
      age->alphawl, age->n/age->t, age->t, age->hfluc, age->invt);
  age->t = 0;
  for ( i = 0; i < n; i++ ) {
    for ( k = 0; k < 3; k++ )
      age->mm[i][k] = 0;
  }
  for ( i = 0; i < n * age->xn; i++ )
    age->hist[i] = 0;
}



/* extensive check of histogram fluctuation */
__inline static int age_wlcheckx(age_t *age,
    double fl, double magred)
{
  age->hfluc = age_calcfl_rms(age);
  if ( !age->invt && age->hfluc < fl ) {
    age_switch(age, magred);
    return 1;
  }
  return 0;
}




/* transition to a neighboring umbrella */
__inline static int age_move(age_t *age, double x, int *id, int local)
{
  int acc = 0, n = age->n;
  int jd;
  double xi, xj, vi, vj, dv, sigi, sigj;

  if ( local ) {
    // jump to a nearest neighbor
    jd = (rand01() < 0.5) ? *id - 1 : *id + 1;
    if ( jd < 0 || jd >= n ) return 0;
  } else {
    // jump any other umbrella
    jd = ( *id + 1 + (int) (rand01() * (n - 1)) ) % n;
  }
  sigi = age->sig[*id];
  sigj = age->sig[ jd];
  /* compute the acceptance probability */
  xi = (x - age->ave[*id]) / sigi;
  xj = (x - age->ave[ jd]) / sigj;
  vi = age->c0[*id] + age->c1[*id] * xi + age->c2[*id] * (xi * xi - 1) / SQRT2;
  vj = age->c0[ jd] + age->c1[ jd] * xj + age->c2[ jd] * (xj * xj - 1) / SQRT2;
  dv = vj - vi;
  //printf("id %d, jd %d, x %g, dv %g(%g), %g(%g), %g; local %d\n",
  //    *id, jd, x, vi, age->c0[*id], vj, age->c0[jd], dv, local);
  acc = metroacc(dv); // ( dv <= 0 || rand01() < exp(-dv) );
  if ( acc ) {
    *id = jd;
  }
  return acc;
}


#endif /* AGE_H__ */
