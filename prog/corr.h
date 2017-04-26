#ifndef CORR_H__
#define CORR_H__



/* module to compute the autocorrelation function
 * of the eigenmodes */



#include "util.h"



#define CORR_BLKSZ  1024



typedef struct {
  int n;
  int cnt; /* number of frames */
  int cap; /* maximal number of frames */
  double *arr; /* data array */
} corr_t;



/* open an correlation function object */
__inline static corr_t *corr_open(int n, int cap)
{
  corr_t *c;

  xnew(c, 1);
  c->n = n;
  c->cnt = 0;
  if ( cap <= 0 ) {
    cap = CORR_BLKSZ;
  }
  c->cap = cap;
  xnew(c->arr, c->cap * c->n);

  return c;
}



/* close the correlation function object */
__inline static void corr_close(corr_t *c)
{
  free(c->arr);
  free(c);
}



/* add a frame */
__inline static void corr_add(corr_t *c, double *u)
{
  int i, n = c->n;

  if ( c->cnt >= c->cap ) {
    c->cap += CORR_BLKSZ;
    /* make the array larger */
    xrenew(c->arr, c->cap * c->n);
  }

  /* add the new entry */
  for ( i = 0; i < n; i++ ) {
    c->arr[c->cnt * n + i] = u[i];
  }
  c->cnt += 1;
}



__inline static void corr_getave(corr_t *c, double *uave)
{
  int i, k, n = c->n;

  for ( i = 0; i < n; i++ ) {
    uave[i] = 0;
  }

  for ( k = 0; k < c->cnt; k++ ) {
    for ( i = 0; i < n; i++ ) {
      uave[i] += c->arr[k * n + i];
    }
  }

  for ( i = 0; i < n; i++ ) {
    uave[i] /= c->cnt;
  }
}



/* compute autocorrelation functions
 * for a separation of j frames */
__inline static int corr_compute(corr_t *c, double *uu, int j,
    const double *uave)
{
  int i, k, k2, n = c->n, tmax = c->cnt;
  double u1, u2;

  /* reset the accumulators */
  for ( i = 0; i < n; i++ ) {
    uu[i] = 0;
  }

  /* loop over the starting point of the window
   * to compute the autocorrelation functions
   * of k and k + j */
  for ( k = 0; k < tmax - j; k++ ) {
    k2 = k + j;
    for ( i = 0; i < n; i++ ) {
      u1 = c->arr[k * n + i];
      u2 = c->arr[k2 * n + i];
      if ( uave != NULL ) {
        u1 -= uave[i];
        u2 -= uave[i];
      }
      uu[i] += u1 * u2;
    }
  }

  /* normalize */
  for ( i = 0; i < n; i++ ) {
    uu[i] /= tmax - j;
  }

  return 0;
}



/* compute and save autocorrelation functions
 * `tol` is the tolerance level as a fraction of the peak value
 * of the correlation functions at time zero, below which
 * the correlation functions are assumed to be zero
 * if `subave` is true, the time-average value is subtracted
 * from c->arr[] before computing the correlation function
 * otherwise, zero is assumed as the average */
__inline static int corr_save(corr_t *c, int dt,
    double tmax, double tol,
    int subave, const char *fn)
{
  int i, j, jmax, n = c->n, err;
  int *stopped = NULL;
  double *uu0, *uu, *uave = NULL, *gam;
  FILE *fp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  xnew(uu0, n);
  xnew(uu, n);
  xnew(stopped, n);
  xnew(gam, n);

  for ( i = 0; i < n; i++ ) {
    stopped[i] = 0;
  }

  if ( subave ) {
    /* deduct the average values before computing correlations
     * in most cases, we test on model systems for which
     * the averages are zeroes, so there is no need to
     * deduct the averages to shift the baseline */
    xnew(uave, n);
    corr_getave(c, uave);
  }

  /* save the heading */
  fprintf(fp, "# %d %d %d\n", n, dt, c->cnt);

  /* try to determine the maximal separation
   * for correlation functions */
  if ( tmax <= 0 ) {
    jmax = c->cnt;
  } else {
    jmax = (int) ( tmax / dt );
  }

  /* loop over separations */
  for ( j = 0; j < jmax; j++ ) {
    /* compute autocorrelation functions
     * at a separation of j frames */
    corr_compute(c, uu, j, uave);

    if ( j == 0 ) {
      /* save the static fluctuation as a reference point */
      for ( i = 0; i < n; i++ ) {
        uu0[i] = uu[i];
        gam[i] = uu[i] * dt;
      }
    } else {
      /* decide if to stop the calculation */
      err = 0;
      for ( i = 0; i < n; i++ ) {
        if ( uu[i] <= uu0[i] * tol ) {
          stopped[i] = 1;
        }
        gam[i] += 2 * uu[i] * dt;
        if ( stopped[i] ) {
          err += 1;
          break;
        }
      }

      if ( err >= n ) break;
    }

    fprintf(fp, "%d", j * dt);
    for ( i = 0; i < n; i++ ) {
      fprintf(fp, "\t%g", uu[i]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "#gamma");
  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "\t%g", gam[i]);
  }
  fprintf(fp, "\n");

  fprintf(stderr, "autocorrelation functions saved in %s, %d frames\n",
      fn, j);
  fclose(fp);

  free(uave);
  free(uu0);
  free(uu);
  free(stopped);
  free(gam);

  return 0;
}



#if 0
/* the following function is now deprecated,
 * the code of the same functionality is integrated
 * at the end of `premeta()` in `invt.c` */

/* compute the integrals of the autocorrelation functions
 * `tol` is the tolerance level as a fraction of the peak value
 * of the correlation functions at time zero, below which
 * the correlation functions are assumed to be zero
 * otherwise, zero is assumed as the average */
__inline static int corr_getgamma(corr_t *c, double *gam,
    double alpha, const double *lam)
{
  int i, n = c->n;
  int *stopped;
  double *uave;

  xnew(uave, n);
  xnew(stopped, n);

  for ( i = 0; i < n; i++ ) {
    gam[i] = 0;
    stopped[i] = 0;
  }

  corr_getave(c, uave);

  /* compute the static fluctuation */
  corr_compute(c, gam, 0, uave);

  for ( i = 0; i < n; i++ ) {
    gam[i] *= 2 / alpha;
    if ( lam != NULL && fabs(lam[i]) > 0 ) {
      gam[i] /= lam[i];
    }
  }

  free(uave);
  free(stopped);

  return 0;
}
#endif



/* print the thermodynamic averages of the modes */
__inline static int corr_printfluc(corr_t *c, int subave,
    const double *xerr)
{
  int i, n = c->n;
  double *uu, *uave = NULL;

  xnew(uu, n);

  if ( subave ) {
    /* deduct the average values before computing correlations
     * in most cases, we test on model systems for which
     * the averages are zeroes, so there is no need to
     * deduct the averages to shift the baseline */
    xnew(uave, n);
    corr_getave(c, uave);
  }

  /* compute autocorrelation functions
   * at a separation of 0 frames */
  corr_compute(c, uu, 0, uave);

  for ( i = 0; i < n; i++ ) {
    printf("%4d\t%15.5e", i + 1, uu[i]);
    if ( xerr != NULL ) {
      printf("\t%15.5e", xerr[i]);
    }
    printf("\n");
  }

  free(uave);
  free(uu);

  return 0;
}



#endif /* CORR_H__ */
