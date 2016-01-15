#ifndef CORR_H__
#define CORR_H__



/* module to compute autocorrelation function */



#include "util.h"



#define CORR_BLKSZ  1024



typedef struct {
  int n;
  int cnt; /* number of frames */
  int cap; /* maximal number of frames */
  double *arr; /* data array */
} corr_t;



/* open an correction function object */
static corr_t *corr_open(int n, int cap)
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



/* close the correction function object */
static void corr_close(corr_t *c)
{
  free(c->arr);
  free(c);
}



/* add a frame */
static void corr_add(corr_t *c, double *u)
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



static void corr_getave(corr_t *c, double *uave)
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
static int corr_compute(corr_t *c, double *uu, int j,
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
      u1 = c->arr[k * n + i] - uave[i];
      u2 = c->arr[k2 * n + i] - uave[i];
      uu[i] += u1 * u2;
    }
  }

  /* normalize */
  for ( i = 0; i < n; i++ ) {
    uu[i] /= tmax - j;
  }

  return 0;
}



/* compute and save autocorrelation functions */
static int corr_save(corr_t *c, int dt, double tol, const char *fn)
{
  int i, j, n = c->n;
  double *uu0, *uu, *uave;
  FILE *fp;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  xnew(uave, n);
  xnew(uu0, n);
  xnew(uu, n);

  /* compute the averages */
  corr_getave(c, uave);

  fprintf(fp, "# %d %d %d\n", n, dt, c->cnt);

  for ( j = 0; j < c->cnt; j++ ) {
    /* compute correlation functions
     * at a separation of j frames */
    corr_compute(c, uu, j, uave);

    if ( j == 0 ) {
      /* save the static fluctuation as a reference point */
      for ( i = 0; i < n; i++ ) {
        uu0[i] = uu[i];
      }
    } else {
      int stop = 1;

      /* decide if to stop the calculation */
      for ( i = 0; i < n; i++ ) {
        if ( uu[i] > uu0[i] * tol ) {
          stop = 0;
          break;
        }
      }

      if ( stop ) break;
    }

    fprintf(fp, "%d", j * dt);
    for ( i = 0; i < n; i++ ) {
      fprintf(fp, "\t%g", uu[i]);
    }
    fprintf(fp, "\n");
  }

  fprintf(stderr, "autocorrelation functions saved in %s, %d frames\n", fn, j);
  fclose(fp);
  free(uave);
  free(uu0);
  free(uu);

  return 0;
}



#endif /* CORR_H__ */
