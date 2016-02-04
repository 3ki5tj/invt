/* predict the error over the range */



/* allow additional parameters for c scan */
#define SCAN
#include "invt.h"
#include "intq.h"



static void invt_geterr(invtpar_t *m)
{
  double t, err1, err2;

  t = (double) m->nsteps;

  /* compute the optimal c and the error for 1/t formula */
  m->c = estbestc(t, 0, m->alpha0, m->n, m->winn, m->win,
      m->sampmethod, 0, &err1, 0);

  /* compute the exact minimal error under the same condition */
  err2 = opterror_ez(m->c, t, m->alpha0,
      m->alpha_nint, m->fnalpha,
      m->n, m->winn, m->win, m->sampmethod,
      NULL, 0);

  printf("c %g, err %g (invt), %g (exact)\n", m->c, err1, err2);
}



static void invt_scanc(invtpar_t *m)
{
  double c, t, t0;
  double err, err0, err1, err2;

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

    /* compute the exact minimal error under the same condition */
    err2 = opterror_ez(c, t, m->alpha0,
        m->alpha_nint, m->fnalpha,
        m->n, m->winn, m->win, m->sampmethod,
        NULL, 0);

    printf("%g\t%g\t%g\t%g\t%g\n",
        c, err, err0, err1, err2);
  }
}



static void invt_scannb(invtpar_t *m)
{
  double t, t0;
  double nb, c, err1, err1norm, err2, err2norm;

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
        m->sampmethod, 0, &err1, 0);

    t0 = c / m->alpha0;
    err1norm = err1 * sqrt(t + t0);

    /* compute the exact minimal error under the same condition */
    err2 = opterror_ez(c, t, m->alpha0,
        m->alpha_nint, m->fnalpha,
        m->n, m->winn, m->win, m->sampmethod,
        NULL, 0);

    err2norm = err2 * sqrt(t + t0);

    printf("%+8.5f\t%8.3f\t%g\t%g\t%g\t%g\n",
        nb, c, err1, err1norm, err2, err2norm);
  }
}



static void invt_scansig(invtpar_t *m)
{
  double t, t0;
  double sig, c, err1, err1norm, err2, err2norm;

  t = (double) m->nsteps;

  for ( sig = m->sigmin; sig < m->sigmax + 0.001 * m->sigdel; sig += m->sigdel ) {
    /* make the window */
    m->wingaus = sig;
    invtpar_mkgauswin(m);

    /* find the optimal c */
    c = estbestc(t, 0, m->alpha0, m->n, m->winn, m->win,
        m->sampmethod, 0, &err1, 0);

    t0 = c / m->alpha0;
    err1norm = err1 * sqrt(t + t0);

    /* compute the exact minimal error under the same condition */
    err2 = opterror_ez(c, t, m->alpha0,
        m->alpha_nint, m->fnalpha,
        m->n, m->winn, m->win, m->sampmethod,
        NULL, 0);

    err2norm = err2 * sqrt(t + t0);

    printf("%8.5f\t%8.3f\t%g\t%g\t%g\t%g\n",
        sig, c, err1, err1norm, err2, err2norm);
  }
}



int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);

  if ( !m->cscan && !m->nbscan && !m->sigscan ) {
    invt_geterr(m);
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
