/* decomposition of histogram fluctuation into eigenmodes
 * for translationally-invariant (homogeneous) updating schemes */


"use strict";


/* compute the eigenvectors,
   costab = n * phi */
function mkcostab(n, pbc)
{
  var i, k, costab, norm = 1.0, ang, y;

  costab = new Array(n * n);

  /* first row, uniform distribution
   * costab(0, i) = n phi(0, i) = 1 */
  for ( i = 0; i < n; i++ ) {
    costab[i] = norm;
  }

  norm *= Math.sqrt(2);
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
      ang = k * 2 * Math.PI / n;
      for ( i = 0; i < n; i++ ) {
        if ( k <= n/2 ) {
          y = Math.cos(i * ang);
        } else {
          y = Math.sin(i * ang);
        }
        costab[k*n + i] = norm * y;
      }
    }
  } else {
    /* non-periodic case, for k = 1, ..., n - 1
     * n phi(k, i) = sqrt(2) cos[k (i+1/2) Pi / n] */
    for ( k = 1; k < n; k++ ) {
      ang = k * Math.PI / n;
      for ( i = 0; i < n; i++ ) {
        costab[k*n + i] = norm * Math.cos((i + 0.5) * ang);
      }
    }
  }

  return costab;
}



/* decompose histogram visit to id, into eigenmodes, `u`
 * the coefficients are pre-computed in costab */
function getcosmodesh(id, n, u, costab)
{
  for ( var k = 0; k < n; k++ )
    u[k] = costab[k*n + id] / n;
}



/* decompose a general fluctuation, `v`, into eigenmodes, `u`
 * the coefficients are pre-computed in costab */
function getcosmodes(v, n, u, costab)
{
  /* this is a simple implementation
   * use FFT for performance */
  for ( var k = 0; k < n; k++ ) {
    var s = 0;
    for ( var i = 0; i < n; i++ ) {
      s += costab[k*n + i] * v[i];
    }
    u[k] = s / n;
  }
}



/* compute the vector from the mode coefficients */
function fromcosmodes(v, n, u, costab)
{
  var i, k;

  for ( i = 0; i < n; i++ ) {
    v[i] = 0;
  }

  for ( k = 0; k < n; k++ ) {
    for ( i = 0; i < n; i++ ) {
      v[i] += costab[k*n + i] * u[k];
    }
  }
}

