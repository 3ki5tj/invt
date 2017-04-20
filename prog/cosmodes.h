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



/* decompose the histogram fluctuation, `v`, into eigenmodes, `u`
 * the coefficients are pre-computed in costab */
static void getcosmodes(double *v, int n, double *u,
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



#endif /* COSMODES_H__ */
