/* predict the error over the range */



/* allow additional parameters for c scan */
#define SCAN
#include "invt.h"
#include "intq.h"



static void invt_geterr(invtpar_t *m,
    const double *lambda, const double *gamma)
{
  double t, qt = 0, err1, err2;
  intq_t *intq;

  t = (double) m->nsteps;

  /* compute the optimal c and the error for 1/t formula */
  m->c = estbestc_invt(t, m->alpha0, m->n, lambda, gamma,
      0, &err1, 0);

  /* compute the exact minimal error under the same condition */
  err2 = esterror_opt(t, m->alpha0, &qt, m->qprec,
      m->alpha_nint, &intq, m->n, NULL,
      lambda, gamma, m->verbose);

  /* save the optimal schedule to file */
  intq_save(intq, m->c, m->c / m->alpha0, m->fnalpha);

  printf("c %g, err %g, sqr %g (invt), %g, sqr %g (exact), %s\n",
      m->c, err1, err1 * err1, err2, err2 * err2, m->fnalpha);

  intq_close(intq);
}



static void invt_scanc(invtpar_t *m,
    const double *lambda, const double *gamma)
{
  double c, qt = 0, t, t0;
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
    err2 = esterror_opt(t, m->alpha0, &qt, m->qprec,
        m->alpha_nint, NULL, m->n, NULL,
        lambda, gamma, m->verbose);

    printf("%g\t%g\t%g\t%g\t%g\n",
        c, err, err0, err1, err2);
  }
}



static void invt_scannb(invtpar_t *m,
    const double *gamma)
{
  double t, t0;
  double nb, c, qt = 0, err1, err1norm, err2, err2norm;
  double *lambda;

  t = (double) m->nsteps;

  for ( nb = m->nbmin; nb < m->nbmax + 0.001 * m->nbdel; nb += m->nbdel ) {
    if ( fabs(nb) < DBL_EPSILON * 100 ) {
      nb = 0;
    }

    /* initialize the window */
    m->winn = 2;
    m->win[1] = nb;
    m->win[0] = 1 - nb;
    lambda = geteigvals(m->n, m->winn, m->win, m->pbc,
        0, NULL, 1);

    /* find the optimal c */
    c = estbestc_invt(t, m->alpha0, m->n, lambda, gamma,
        0, &err1, 0);

    t0 = c / m->alpha0;
    err1norm = err1 * sqrt(t + t0);

    /* compute the exact minimal error under the same condition */
    err2 = esterror_opt(t, m->alpha0, &qt, m->qprec,
        m->alpha_nint, NULL, m->n, NULL,
        lambda, gamma, m->verbose);

    err2norm = err2 * sqrt(t + t0);

    printf("%+8.5f\t%10.6f\t%10.6f\t%g\t%10.6f\t%g\n",
        nb, c, err1, err1norm, err2, err2norm);
  }
}



static void invt_scansig(invtpar_t *m,
    const double *gamma)
{
  double t, t0;
  double sig, c, qt = 0, err1, err1norm, err2, err2norm;
  double *lambda;

  t = (double) m->nsteps;

  for ( sig = m->sigmin; sig < m->sigmax + 0.001 * m->sigdel; sig += m->sigdel ) {
    /* make the window */
    m->gaussig = sig;
    invtpar_mkgauswin(m);
    //lambda = geteigvals(m->n, m->winn, m->win, m->pbc, 0, NULL, 1);
    lambda = trimwindow(m->n, &m->winn, m->win, m->pbc, 0);

    /* find the optimal c, according to the inverse time schedule */
    c = estbestc_invt(t, m->alpha0, m->n, lambda, gamma,
        0, &err1, 0);

    /* compute the time-normalized error
     * of the optimized inverse-time schedule */
    t0 = c / m->alpha0;
    err1norm = err1 * sqrt(t + t0);

    /* compute the exact minimal error under the same condition */
    err2 = esterror_opt(t, m->alpha0, &qt, m->qprec,
        m->alpha_nint, NULL, m->n, NULL,
        lambda, gamma, m->verbose);

    /* compute the error of the optimal schedule */
    err2norm = err2 * sqrt(t + t0);

    printf("%8.5f\t%10.6f\t%10.6f\t%g\t%10.6f\t%g\n",
        sig, c, err1, err1norm, err2, err2norm);
    free(lambda);
  }
}



int main(int argc, char **argv)
{
  invtpar_t m[1];
  double *lambda = NULL, *gamma = NULL;

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);

  lambda = geteigvals(m->n, m->winn, m->win, m->pbc,
      0, NULL, 1);

  /* estimate or load the gamma values */
  gamma = estgamma(m->n, m->sampmethod, m->pbc, m->localg);
  if ( m->fngamma[0] != '\0' ) {
    loadgamma(m->n, gamma, m->fngamma);
  }

  if ( !m->cscan && !m->nbscan && !m->sigscan ) {
    invt_geterr(m, lambda, gamma);
  }

  if ( m->cscan ) {
    invt_scanc(m, lambda, gamma);
  }

  if ( m->nbscan ) {
    invt_scannb(m, gamma);
  }

  if ( m->sigscan ) {
    invt_scansig(m, gamma);
  }

  free(lambda);
  free(gamma);

  invtpar_finish(m);

  return 0;
}
