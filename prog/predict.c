/* predict the error over the range */



/* allow additional parameters for c scan */
#define SCAN
#include "invt.h"



static void invt_scanc(invtpar_t *m)
{
  double c, t, t0;
  double err, err0, err1;

  t = (double) m->nsteps;
  for ( c = m->cmin; c < m->cmax + 0.001 * m->cdel; c += m->cdel ) {
    t0 = c / m->alpha0;

    /* c / (t + t0) */
    err = esterror_ez(c, t, t0, m->alpha0,
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



static void invt_scannb(invtpar_t *m)
{
  double t, t0;
  double nb, c, err, errnorm;

  t = (double) m->nsteps;

  for ( nb = m->nbmin; nb < m->nbmax + 0.001 * m->nbdel; nb += m->nbdel ) {
    if ( fabs(nb) < DBL_EPSILON * 100 ) {
      nb = 0;
    }

    /* initialize the window */
    m->winn = 2;
    m->win[1] = nb;
    m->win[0] = 1 - nb;

    /* find the optimal c */
    c = estbestc(t, 0, m->alpha0, m->n, m->winn, m->win,
        m->sampmethod, 0, &err, 0);

    t0 = c / m->alpha0;
    errnorm = err * sqrt(t + t0);

    printf("%+8.5f\t%8.3f\t%g\t%g\n",
        nb, c, err, errnorm);
  }
}



static void invt_scansig(invtpar_t *m)
{
  double t, t0;
  double sig, c, err, errnorm;

  t = (double) m->nsteps;

  for ( sig = m->sigmin; sig < m->sigmax + 0.001 * m->sigdel; sig += m->sigdel ) {
    /* make the window */
    m->wingaus = sig;
    invtpar_mkgauswin(m);

    /* find the optimal c */
    c = estbestc(t, 0, m->alpha0, m->n, m->winn, m->win,
        m->sampmethod, 0, &err, 0);

    t0 = c / m->alpha0;
    errnorm = err * sqrt(t + t0);

    printf("%8.5f\t%8.3f\t%g\t%g\n",
        sig, c, err, errnorm);
  }
}



int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);

  if ( !m->cscan && !m->nbscan && !m->sigscan ) {
    double t, err1, err2;

    t = (double) m->nsteps;
    m->c = estbestc(t, 0, m->alpha0, m->n, m->winn, m->win,
        m->sampmethod, 0, &err1, 0);
    err2 = opterror_ez(m->c, t, m->alpha0,
        1000, "opta.dat",
        m->n, m->winn, m->win, m->sampmethod);

    printf("c %g, err %g, %g\n", m->c, err1, err2);
  }

  if ( m->cscan ) {
    invt_scanc(m);
  }

  if ( m->nbscan ) {
    invt_scannb(m);
  }

  if ( m->sigscan ) {
    invt_scansig(m);
  }

  invtpar_finish(m);

  return 0;
}
