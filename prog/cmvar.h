#ifndef CMVAR_H__
#define CMVAR_H__


/* accumulator for computing the correlation integrals
 * of the eigenmodes from the variances of the bias potential */
#include "cosmodes.h"



typedef struct {
  int n;
  long cnt;
  double *u, *usum, *usqr, *uave, *uvar;
  double *costab;
} cmvar_t;

__inline static void cmvar_clear(cmvar_t *cm)
{
  int i;
  cm->cnt = 0;
  for ( i = 0; i < cm->n; i++ ) {
    cm->u[i] = 0;
    cm->usum[i] = 0;
    cm->usqr[i] = 0;
    cm->uave[i] = 0;
    cm->uvar[i] = 0;
  }
}

__inline static cmvar_t *cmvar_open(int n, int pbc)
{
  cmvar_t *cm;

  xnew(cm, 1);
  cm->n = n;
  cm->costab = mkcostab(n, pbc);
  xnew(cm->u, n);
  xnew(cm->usum, n);
  xnew(cm->usqr, n);
  xnew(cm->uave, n);
  xnew(cm->uvar, n);
  cmvar_clear(cm);
  return cm;
}

__inline static void cmvar_close(cmvar_t *cm)
{
  free(cm->u);
  free(cm->usum);
  free(cm->usqr);
  free(cm->uave);
  free(cm->uvar);
  free(cm->costab);
  free(cm);
}

/* deposit a vector into the mode */
__inline static void cmvar_add(cmvar_t *cm, double *v)
{
  double uave, du;
  long cnt = cm->cnt;
  int i, n = cm->n;

  /* cosine transform to the components of v */
  getcosmodes(v, n, cm->u, cm->costab);

  /* update accumulators for average and variance */
  for ( i = 0; i < n; i++ ) {
    uave = (cnt > 0) ? cm->usum[i] / cnt : 0;
    cm->usum[i] += cm->u[i]; /* update the average */
    if ( cnt > 0 ) { /* update the variance */
      du = cm->u[i] - uave;
      cm->usqr[i] += du * du * cnt / (cnt + 1);
    }
  }
  cm->cnt += 1;
}

/* compute the average and variance */
__inline static void cmvar_get(cmvar_t *cm)
{
  int i;

  if ( cm->cnt <= 0 ) return;
  for ( i = 0; i < cm->n; i++ ) {
    cm->uave[i] = cm->usum[i] / cm->cnt;
    cm->uvar[i] = cm->usqr[i] / cm->cnt;
  }
}



#endif /* CMVAR_H__ */
