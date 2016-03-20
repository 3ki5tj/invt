#ifndef COSMODES_H__
#define COSMODES_H__


/* decomposition of updating scheme to cosine modes */



/* compute the cosine table
   costab = phi * n */
static double *mkcostab(int n)
{
  int i, k;
  double *costab, norm = 1.0, ang;

  xnew(costab, n * n);

  for ( i = 0; i < n; i++ ) {
    costab[i] = norm;
  }

  norm *= sqrt(2);
  ang = M_PI / n;
  for ( k = 1; k < n; k++ ) {
    for ( i = 0; i < n; i++ ) {
      costab[k*n + i] = norm * cos((i + 0.5) * k * ang);
    }
  }

  return costab;
}



/* compute the magnitude of the eigenmodes of the histogram
 * save them to `u` */
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



#if 0
/* compute the adjusted instantaneous histogram from the bin index i
   z_i = h_i / p_i - 1 */
__inline static void mkinsthist(double *ih, int n, int i)
{
  int j;

  for ( j = 0; j < n; j++ ) {
    ih[j] = -1;
  }
  ih[i] += n;
}
#endif


#endif /* COSMODES_H__ */
