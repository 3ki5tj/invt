/* predict the error over the range */



/* allow additional parameters for c scan */
#define SCAN
#include "invt.h"
#include "intq.h"



static void invt_geterr(invtpar_t *m,
    const double *lambda, const double *gamma)
{
  double t, err1, err2;

  t = (double) m->nsteps;

  /* compute the optimal c and the error for 1/t formula */
  m->c = estbestc_invt(t, m->alpha0, m->n, lambda, gamma,
      0, &err1, 0);

  /* compute the exact minimal error under the same condition */
  err2 = esterror_opt(t,
      intq_getqt(t, m->c, m->c / m->alpha0), m->alpha0,
      m->alpha_nint, NULL, m->n, NULL,
      lambda, gamma);

  printf("c %g, err %g (invt), %g (exact)\n", m->c, err1, err2);
}



static void invt_scanc(invtpar_t *m,
    const double *lambda, const double *gamma)
{
  double c, t, t0;
  double err, err0, err1, err2;

  t = (double) m->nsteps;
  for ( c = m->cmin; c < m->cmax + 0.001 * m->cdel; c += m->cdel ) {
    t0 = c / m->alpha0;

    /* c / (t + t0) */
    err = esterror_invt(t, c, m->alpha0, m->n, NULL,
        lambda, gamma);

    /* initial equilibrium value */
    err0 = esterror_eql(m->alpha0, m->n, NULL,
        lambda, gamma);

    /* final equilibrium value */
    err1 = esterror_eql(c / (t0 + t), m->n, NULL,
        lambda, gamma);

    /* compute the exact minimal error under the same condition */
    err2 = esterror_opt(t,
        intq_getqt(t, m->c, m->t0), m->alpha0,
        m->alpha_nint, NULL, m->n, NULL,
        lambda, gamma);

    printf("%g\t%g\t%g\t%g\t%g\n",
        c, err, err0, err1, err2);
  }
}



static void invt_scannb(invtpar_t *m,
    const double *lambda, const double *gamma)
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
    c = estbestc_invt(t, m->alpha0, m->n, lambda, gamma,
        0, &err1, 0);

    t0 = c / m->alpha0;
    err1norm = err1 * sqrt(t + t0);

    /* compute the exact minimal error under the same condition */
    err2 = esterror_opt(t,
        intq_getqt(t, c, c / m->alpha0), m->alpha0,
        m->alpha_nint, NULL, m->n, NULL,
        lambda, gamma);

    err2norm = err2 * sqrt(t + t0);

    printf("%+8.5f\t%8.3f\t%g\t%g\t%g\t%g\n",
        nb, c, err1, err1norm, err2, err2norm);
  }
}



static void invt_scansig(invtpar_t *m,
    const double *lambda, const double *gamma)
{
  double t, t0;
  double sig, c, err1, err1norm, err2, err2norm;

  t = (double) m->nsteps;

  for ( sig = m->sigmin; sig < m->sigmax + 0.001 * m->sigdel; sig += m->sigdel ) {
    /* make the window */
    m->wingaus = sig;
    invtpar_mkgauswin(m);

    /* find the optimal c */
    c = estbestc_invt(t, m->alpha0, m->n, lambda, gamma,
        0, &err1, 0);

    t0 = c / m->alpha0;
    err1norm = err1 * sqrt(t + t0);

    /* compute the exact minimal error under the same condition */
    err2 = esterror_opt(t,
        intq_getqt(t, c, c / m->alpha0), m->alpha0,
        m->alpha_nint, NULL, m->n, NULL,
        lambda, gamma);

    err2norm = err2 * sqrt(t + t0);

    printf("%8.5f\t%8.3f\t%g\t%g\t%g\t%g\n",
        sig, c, err1, err1norm, err2, err2norm);
  }
}



int main(int argc, char **argv)
{
  invtpar_t m[1];
  double *lambda = NULL, *gamma = NULL;

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);

  lambda = geteigvals(m->n, m->winn, m->win,
      0, NULL, 1);
  gamma = estgamma(m->n, m->sampmethod);

  if ( !m->cscan && !m->nbscan && !m->sigscan ) {
    invt_geterr(m, lambda, gamma);
  }

  if ( m->cscan ) {
    invt_scanc(m, lambda, gamma);
  }

  if ( m->nbscan ) {
    invt_scannb(m, lambda, gamma);
  }

  if ( m->sigscan ) {
    invt_scansig(m, lambda, gamma);
  }

  free(lambda);
  free(gamma);

  invtpar_finish(m);

  return 0;
}
