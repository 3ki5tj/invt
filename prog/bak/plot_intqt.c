/* make a plot of the integrand versus q as an integral */
#include "invt.h"

static void plot(invtpar_t *m, int npt)
{
  int j, k;
  double *lambda, *gamma;
  double t, sum, y, q, qt;
  double *arr;

  lambda = trimwindow(m->n, &m->winn, m->win, 0);
  gamma = estgamma(m->n, m->sampmethod);

  t = (double) m->nsteps;
  //qt = intq_getqt(t, m->c, m->t0);
  qt = m->c * log( 1 + t / m->t0 );
  
  xnew(arr, npt + 1);

  for ( j = 0; j <= npt; j++ ) {
    q = qt * j / npt;
    sum = 0;
    for ( k = 0; k < m->n; k++ ) {
      y = lambda[k] * exp( lambda[k] * (q - qt) );
      sum += gamma[k] * y * y;
    }
    y = sqrt( sum );
    printf("%g %g\n", q, y);
  }

  free(lambda);
  free(gamma);
  free(arr);
}

int main(int argc, char **argv)
{
  invtpar_t m[1];
  
  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  
  plot(m, 100);

  invtpar_finish(m);

  return 0;
}
