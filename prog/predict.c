/* predict the error over the range */



/* allow additional parameters for c scan */
#define SCAN
#include "invt.h"
#include "intq.h"



/* save the error components to file */
static int save_xerr(invtpar_t *m, const char *fn,
    const double *xerri, const double *xerrf,
    const double *xerrf_r, const double *xerrf_a,
    const double *lambda, const double *gamma,
    double qT, double a0)
{
  FILE *fp;
  int i, n = m->n;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d qT %.6f a0 %.8e nint %d\n",
      n, qT, a0, m->alpha_nint);
  for ( i = 1; i < n; i++ ) {
    fprintf(fp, "%4d %18.8e %18.8e %18.8e %18.8e %10.8f %18.8f\n",
        i, xerri[i], xerrf[i], xerrf_r[i], xerrf_a[i],
        lambda[i], gamma[i]);
  }

  fclose(fp);
  return 0;
}



static void invt_geterr(invtpar_t *m,
    const double *gamma)
{
  int n = m->n;
  double T, c0, t0, qT, inita, err0, err1, err2;
  double *lambda;
  double *xerri, *xerrf, *xerrf_r, *xerrf_a;
  intq_t *intq;

  xnew(xerri, n);
  xnew(xerrf, n);
  xnew(xerrf_r, n);
  xnew(xerrf_a, n);

  T = (double) m->nsteps;

  /* compute the eigenvalues of the updating matrix */
  lambda = geteigvals(m->n, m->winn, m->win, m->pbc,
      0, NULL, 1);
  /* save updating kernel or window function */
  if ( m->fnwin[0] != '\0' ) {
    savewin(m->winn, m->win, m->fnwin);
  }
  /* save the updating matrix */
  if ( m->fnwinmat[0] != '\0' ) {
    savewinmat(m->winn, m->win, m->n, m->pbc, m->fnwinmat);
  }

  /* initial errors */
  esterror_eql(m->alpha0, m->n, xerri, lambda, gamma);

  c0 = m->c;
  err0 = esterror_invt(T, m->c, m->alpha0, m->t0, m->n,
      NULL, lambda, gamma);

  /* compute the optimal c and the error
   * for the inverse-time formula */
  m->c = estbestc_invt(T, m->alpha0, 0, m->n, lambda, gamma,
      0, &err1, 0);

  /* compute the exact minimal error under the same condition */
  qT = m->qT;
  err2 = esterror_opt(T, m->alpha0, m->initalpha, &qT, m->qprec,
      m->alpha_nint, &intq, m->n, m->kcutoff, m->pbc,
      lambda, gamma, m->verbose);
  inita = intq_getinita(intq);

  /* save the optimal schedule to file */
  t0 = 2 / m->alpha0;
  intq_save(intq, m->c, t0, m->alpha_resample, m->fnalpha);

  printf("c %g, t0 %g, err %g, sqr %g (invt), "
         "opt. c %g, err %g, sqr %g (opt. invt), "
         "qT %g, a(0) %g, err %g, sqr %g (exact), %s\n",
      c0, t0, err0, err0 * err0,
      m->c, err1, err1 * err1,
      qT, inita, err2, err2 * err2, m->fnalpha);

  if ( m->fnxerr[0] != '\0' ) {
    intq_errcomp(intq, m->alpha0, qT, xerrf, xerrf_r, xerrf_a);
    save_xerr(m, m->fnxerr, xerri, xerrf, xerrf_r, xerrf_a,
        lambda, gamma, qT, inita);
  }

  intq_close(intq);
  free(lambda);
  free(xerri);
  free(xerrf);
  free(xerrf_r);
  free(xerrf_a);
}



static void invt_scanc(invtpar_t *m,
    const double *gamma)
{
  double c, T, t0, qT;
  double err, erri, errf, err2;
  double *lambda;

  T = (double) m->nsteps;

  /* compute the eigenvalues of the updating matrix */
  lambda = geteigvals(m->n, m->winn, m->win, m->pbc,
      0, NULL, 1);

  /* initial equilibrium error */
  erri = esterror_eql(m->alpha0, m->n, NULL,
      lambda, gamma);

  /* compute the exact minimal error under the same condition */
  qT = m->qT;
  err2 = esterror_opt(T, m->alpha0, m->initalpha, &qT, m->qprec,
      m->alpha_nint, NULL, m->n, m->kcutoff, m->pbc,
      lambda, gamma, m->verbose);

  /* print out a header */
  printf("# c     \t  final error\t  init. error\t  final equil\t  optimal error\n");

  for ( c = m->cmin; c < m->cmax + 0.001 * m->cdel; c += m->cdel ) {
    /* alpha(t) = c / (t + t0) */
    err = esterror_invt(T, c, m->alpha0, 0, m->n, NULL,
        lambda, gamma);

    /* final equilibrium value */
    t0 = 2 / m->alpha0;
    errf = esterror_eql(c / (t0 + T), m->n, NULL,
        lambda, gamma);

    printf("%8.5f\t%14.10f\t%14.10f\t%14.10f\t%10.6f\n",
        c, err, erri, errf, err2);
  }

  free(lambda);
}



/* scanning the initial updating magntidue */
static void invt_scania(invtpar_t *m,
    const double *gamma)
{
  double T, qT;
  double inita, iafac;
  double err, erri;
  double *lambda;
  intq_t *intq;

  T = (double) m->nsteps;

  /* compute the eigenvalues of the updating matrix */
  lambda = geteigvals(m->n, m->winn, m->win, m->pbc,
      0, NULL, 1);

  /* initial equilibrium error */
  erri = esterror_eql(m->alpha0, m->n, NULL,
      lambda, gamma);

  intq = intq_open(T, m->alpha_nint, m->n,
      m->kcutoff, m->pbc, lambda, gamma);

  /* print out a header */
  printf("# inita \t  final error\t  init. error\t   qT\n");

  iafac = pow(10.0, m->iadel);
  for ( inita = m->iamin; inita < m->iamax * iafac; inita *= iafac ) {
    /* compute the exact minimal error  */
    qT = intq_optqT(intq, inita, m->qprec, m->verbose);
    err = intq_geterr(intq, m->alpha0, qT);

    printf("%14.6e\t%14.10f\t%14.10f\t%14.6e\n",
        inita, err, erri, qT);
  }

  intq_close(intq);
  free(lambda);
}



/* possible scan types */
enum {
  SCAN_NB,
  SCAN_SIG,
  SCAN_OK
};



static void invt_scan(invtpar_t *m,
    const double *gamma, int scantype)
{
  double nb, sig;
  int ok, okmax = m->okmax, eigerr = 0;
  double T, t0, c, qT, err1, err1norm, err2, err2norm;
  double *lambda = NULL;

  T = (double) m->nsteps;

  /* initialize the scanning variable */
  if ( scantype == SCAN_NB )
  {
    nb = m->nbmin;
    printf("# nb    ");
  }
  else if ( scantype == SCAN_SIG )
  {
    sig = m->sigmin;
    printf("# sigma ");
  }
  else if ( scantype == SCAN_OK )
  {
    ok = m->okdel;
    printf("# Kmax  ");
  }
  /* print out a header */
  printf("\t  c-value  \tinvt error\t(normalized)\t  opt. error\t(normalized)\t"
      "   t0    \t   q(T)\n");

  for ( ; ; ) {
    if ( scantype == SCAN_NB )
    {
      /* avoid the round-off error in case
       * we start from a negative value of nb */
      if ( fabs(nb) < DBL_EPSILON * 100 ) {
        nb = 0;
      }

      /* create the updating window */
      m->winn = 2;
      m->win[1] = nb;
      m->win[0] = 1 - nb;
      lambda = geteigvals(m->n, m->winn, m->win, m->pbc,
          0, &eigerr, 1);
    }
    else if ( scantype == SCAN_SIG )
    {
      /* make the Gaussian window */
      m->gaussig = sig;
      invtpar_mkgauswin(m);
      lambda = trimwindow(m->n, &m->winn, m->win, m->pbc, 0, m->verbose);
    }
    else if ( scantype == SCAN_OK )
    {
      /* make the bandpass sinc window */
      m->okmax = ok;
      invtpar_mksinrwin(m);
      lambda = geteigvals(m->n, m->winn, m->win, m->pbc,
          0, &eigerr, 1);
    }

    if ( eigerr > 0 ) {
      fprintf(stderr, "negative eigenvalues encountered\n");
      free(lambda);
    }

    /* find the optimal c for the inverse-time schedule */
    c = estbestc_invt(T, m->alpha0, 0, m->n, lambda, gamma,
        0, &err1, 0);

    /* don't be smart about t0, for it affects
     * the normalization constant, better a constant */
    t0 = 2 / m->alpha0;

    //t0 = 2 * c / m->alpha0;
    err1norm = err1 * sqrt(T + t0);

    /* compute the exact minimal error under the same condition */
    qT = 0; /* use the optimal qT */
    err2 = esterror_opt(T, m->alpha0, m->initalpha, &qT, m->qprec,
        m->alpha_nint, NULL, m->n, m->kcutoff, m->pbc,
        lambda, gamma, m->verbose);

    // t0 = intq_estt0(T, qT);
    err2norm = err2 * sqrt(T + t0);

    /* print the scanning variable */
    if ( scantype == SCAN_NB )
    {
      printf("%+8.5f\t", nb);
    }
    else if ( scantype == SCAN_SIG )
    {
      printf("%+8.5f\t", sig);
    }
    else if ( scantype == SCAN_OK )
    {
      printf("%8d\t", ok);
    }

    printf("%10.6f\t%10.6f\t%10.4f\t%10.6f\t%10.4f\t%10g\t%10g\n",
        c, err1, err1norm, err2, err2norm, t0, qT);

    free(lambda);

    /* update the scanning variable and stop */
    if ( scantype == SCAN_NB )
    {
      nb += m->nbdel;
      if ( nb > m->nbmax + 0.01 * m->nbdel )
        break;
    }
    else if ( scantype == SCAN_SIG )
    {
      sig += m->sigdel;
      if ( sig > m->sigmax + 0.01 * m->sigdel )
        break;
    }
    else if ( scantype == SCAN_OK )
    {
      ok += m->okdel;
      if ( ok > okmax )
        break;
    }
  }
}



int main(int argc, char **argv)
{
  invtpar_t m[1];
  double *gamma = NULL;

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);

  /* estimate or load the gamma values */
  gamma = estgamma(m->n, m->sampmethod, m->pbc, m->localg);
  if ( m->fngamma[0] != '\0' ) {
    if ( m->sampmethod == SAMPMETHOD_MD ) {
      loadgamma(m->n, gamma, m->fngamma);
    } else {
      savegamma(m->n, gamma, m->fngamma);
    }
  }

  if ( !m->cscan
    && !m->iascan
    && !m->nbscan
    && !m->sigscan
    && !m->okscan ) {
    invt_geterr(m, gamma);
  }

  if ( m->cscan ) {
    invt_scanc(m, gamma);
  }

  if ( m->iascan ) {
    invt_scania(m, gamma);
  }

  if ( m->nbscan ) {
    invt_scan(m, gamma, SCAN_NB);
  }

  if ( m->sigscan ) {
    invt_scan(m, gamma, SCAN_SIG);
  }

  if ( m->okscan ) {
    invt_scan(m, gamma, SCAN_OK);
  }

  free(gamma);

  invtpar_finish(m);

  return 0;
}
