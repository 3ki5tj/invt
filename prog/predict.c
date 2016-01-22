/* predict the error over the range */



/* allow additional parameters for c scan */
#define CSCAN
#include "invt.h"



static void invt_predict(invtpar_t *m)
{
  double c, t, t0;
  double err, err0, err1;

  t = (double) m->nsteps;
  for ( c = m->cmin; c < m->cmax + 0.001 * m->cdel; c += m->cdel ) {
    t0 = c / m->alpha0;

    /* c / (t + t0) */
    err = esterror_ez(c, t, t0,
        m->n, m->winn, m->win, m->sampmethod,
        0);

    /* initial equilibrium value */
    err0 = esterror0_ez(m->alpha0,
        m->n, m->winn, m->win, m->sampmethod,
        "initial", 0);

    /* final equilibrium value */
    err1 = esterror0_ez(c / (t0 + t),
        m->n, m->winn, m->win, m->sampmethod,
        "final", 0);

    printf("%g\t%g\t%g\t%g\n",
        c, err, err0, err1);
  }
}



int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  invt_predict(m);
  invtpar_finish(m);

  return 0;
}
