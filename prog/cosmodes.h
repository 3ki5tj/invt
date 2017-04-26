#ifndef COSMODES_H__
#define COSMODES_H__



/* decomposition of histogram fluctuation into eigenmodes
 * for translationally-invariant (homogeneous) updating schemes */



/* compute the eigenvectors,
   costab = n * phi */
static double *mkcostab(int n, int pbc)
{
  int i, k;
  double *costab, norm = 1.0, ang, y;

  xnew(costab, n * n);

  /* first row, uniform distribution
   * costab(0, i) = n phi(0, i) = 1 */
  for ( i = 0; i < n; i++ ) {
    costab[i] = norm;
  }

  norm *= sqrt(2);
  if ( pbc ) {
    /* periodic case, for k = 1, ..., n - 1
     * in the paper, we use the exponential form
     *   n phi(k, i) = exp(I i k 2 Pi / n)
     * here, we use the equivalent cosine and sine form
     *   n phi(k, i) = sqrt(2) cos(i k 2 Pi / n)
     * for k <= n/2, or
     *   n phi(k, i) = sqrt(2) sin(i k 2 Pi / n)
     * for k > n/2 */
    for ( k = 1; k < n; k++ ) {
      ang = k * 2 * M_PI / n;
      for ( i = 0; i < n; i++ ) {
        if ( k <= n/2 ) {
          y = cos(i * ang);
        } else {
          y = sin(i * ang);
        }
        costab[k*n + i] = norm * y;
      }
    }
  } else {
    /* non-periodic case, for k = 1, ..., n - 1
     * n phi(k, i) = sqrt(2) cos[k (i+1/2) Pi / n] */
    for ( k = 1; k < n; k++ ) {
      ang = k * M_PI / n;
      for ( i = 0; i < n; i++ ) {
        costab[k*n + i] = norm * cos((i + 0.5) * ang);
      }
    }
  }

  return costab;
}



/* decompose histogram visit to id, into eigenmodes, `u`
 * the coefficients are pre-computed in costab */
__inline static void getcosmodesh(int id, int n, double *u,
    double *costab)
{
  int k;

  for ( k = 0; k < n; k++ )
    u[k] = costab[k*n + id] / n;
}



/* decompose a general fluctuation, `v`, into eigenmodes, `u`
 * the coefficients are pre-computed in costab */
__inline static void getcosmodes(double *v, int n, double *u,
    double *costab)
{
  int i, k;
  double s;

  /* this is a simple implementation
   * use FFT for performance */
  for ( k = 0; k < n; k++ ) {
    s = 0;
    for ( i = 0; i < n; i++ ) {
      s += costab[k*n + i] * v[i];
    }
    u[k] = s / n;
  }
}



/* compute the vector from the mode coefficients */
static void fromcosmodes(double *v, int n, double *u,
    double *costab)
{
  int i, k;

  for ( i = 0; i < n; i++ ) {
    v[i] = 0;
  }

  for ( k = 0; k < n; k++ ) {
    for ( i = 0; i < n; i++ ) {
      v[i] += costab[k*n + i] * u[k];
    }
  }
}


typedef struct {
  int n;
  long cnt;
  double *u, *usum, *usqr, *uave, *uvar;
  double *costab;
} cmvar_t;

__inline static cmvar_t *cmvar_open(int n, int pbc)
{
  cmvar_t *cm;
  int i;

  xnew(cm, 1);
  cm->n = n;
  cm->costab = mkcostab(n, pbc);
  cm->cnt = 0;
  xnew(cm->u, n);
  xnew(cm->usum, n);
  xnew(cm->usqr, n);
  xnew(cm->uave, n);
  xnew(cm->uvar, n);
  for ( i = 0; i < n; i++ ) {
    cm->u[i] = 0;
    cm->usum[i] = 0;
    cm->usqr[i] = 0;
    cm->uave[i] = 0;
    cm->uvar[i] = 0;
  }
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

  /* shift the baseline of the bias potential */
  for ( uave = 0, i = 0; i < n; i++ ) uave += v[i];
  uave /= n;
  for ( i = 0; i < n; i++ ) v[i] -= uave;

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


#endif /* COSMODES_H__ */
