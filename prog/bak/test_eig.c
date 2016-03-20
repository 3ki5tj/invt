#include <time.h>
#include "util.h"
#include "eig.h"
#include "mtrand.h"



/* compute the eigenvalues of a random tridiagonal matrix */
static void test_eig(int n, double a)
{
  double *d, *e, *v;
  int i;

  xnew(d, n);
  xnew(e, n);
  xnew(v, n * n);

  /* initialize the tridiagonal matrix */
  for ( i = 0; i < n; i++ ) {
    d[i] = 1;
  }
  for ( i = 0; i < n - 1; i++ ) {
    e[i] = 0.5 * ( rand01() * (1 - a) + a );
    d[i] -= e[i];
    d[i+1] -= e[i];
  }
  e[n-1] = 0;

  for ( i = 0; i < n; i++ ) {
    fprintf(stderr, "diag %10.6f, off-diag %10.6f\n", d[i], e[i]);
  }

  /* compute the eigenvalues and eigenvectors */
  eigtriqr(d, e, n, v);
  eigsort(d, v, n);

  /* print out the eigenvalues */
  for ( i = 0; i < n; i++ ) {
    printf("%4d\t%10.6f\n", i, d[i]);
  }

  free(d);
  free(e);
  free(v);
}


int main(int argc, char **argv)
{
  double amp = 0.0;

  mtscramble( clock() );
  if ( argc > 1 ) {
    amp = atof( argv[1] );
  }
  test_eig(200, amp);
  return 0;
}
