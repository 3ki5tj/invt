/* predict the error over the range */



/* allow additional parameters for c scan */
#include "invt.h"
#include "intq.h"


int cscan = 0;
double cmin = 0.1;
double cdel = 0.01;
double cmax = 5.0;

int iascan = 0;
double iamin = 5e-7;
double iadel = 0.2; /* ia is multiplied by 10^iadel each step */
double iamax = 1e-2;

int nbscan = 0;
double nbmin = 0.0;
double nbdel = 0.01;
double nbmax = 0.25;

int sigscan = 0;
double sigmin = 0.0;
double sigdel = 0.2;
double sigmax = 10.0;

int kcscan = 0;
int kcdel = 1;


static void predpar_help(void)
{
  fprintf(stderr, "  --cmin=:       set the minimal c in c-scan, default %g\n", cmin);
  fprintf(stderr, "  --cmax=:       set the maximal c in c-scan, default %g\n", cmax);
  fprintf(stderr, "  --dc=:         set the increment of c in c-scan, default %g\n", cdel);
  fprintf(stderr, "  --iamin=:      set the minimal initial updating magnitude in ia-scan, default %g\n", iamin);
  fprintf(stderr, "  --iamax=:      set the maximal initial updating magnitude in ia-scan, default %g\n", iamax);
  fprintf(stderr, "  --dia=:        set the increment of initial updating magnitude in ia-scan, default %g\n", iadel);
  fprintf(stderr, "  --nbmin=:      set the minimal nearest-neighbor updating magnitude nb in nb-scan, default %g\n", nbmin);
  fprintf(stderr, "  --nbmax=:      set the maximal nearest-neighbor updating magnitude nb in nb-scan, default %g\n", nbmax);
  fprintf(stderr, "  --dnb=:        set the increment of nb in nb-scan, default %g\n", nbdel);
  fprintf(stderr, "  --sigmin=:     set the minimal width of the Gaussian scheme in sig-scan, default %g\n", sigmin);
  fprintf(stderr, "  --sigmax=:     set the maximal width of the Gaussian scheme in sig-scan, default %g\n", sigmax);
  fprintf(stderr, "  --dsig=:       set the increment of the width in sig-scan, default %g\n", sigdel);
  fprintf(stderr, "  --kcscan=:     scan along kc, default %d\n", kcscan);
  fprintf(stderr, "  --dkc=:        set the increment of kc, default %d\n", kcdel);
}

static int predpar_keymatch(invtpar_t *m, const char *key, const char *val)
{
  /* c-scan paramerters */
  if ( strcmpfuzzy(key, "cmin") == 0 )
  {
    cmin = invtpar_getdouble(m, key, val);
    cscan = 1;
  }
  else if ( strcmpfuzzy(key, "cdel") == 0
         || strcmpfuzzy(key, "dc") == 0 )
  {
    cdel = invtpar_getdouble(m, key, val);
    cscan = 1;
  }
  else if ( strcmpfuzzy(key, "cmax") == 0 )
  {
    cmax = invtpar_getdouble(m, key, val);
    cscan = 1;
  }
  /* ia-scan (initial updating magnitude) parameter */
  else if ( strcmpfuzzy(key, "iamin") == 0 )
  {
    iamin = invtpar_getdouble(m, key, val);
    iascan = 1;
  }
  else if ( strcmpfuzzy(key, "iadel") == 0
         || strcmpfuzzy(key, "dia") == 0 )
  {
    iadel = invtpar_getdouble(m, key, val);
    iascan = 1;
  }
  else if ( strcmpfuzzy(key, "iamax") == 0 )
  {
    iamax = invtpar_getdouble(m, key, val);
    iascan = 1;
  }
  /* nb-scan parameters */
  else if ( strcmpfuzzy(key, "nbmin") == 0 )
  {
    nbmin = invtpar_getdouble(m, key, val);
    nbscan = 1;
  }
  else if ( strcmpfuzzy(key, "nbdel") == 0
         || strcmpfuzzy(key, "dnb") == 0 )
  {
    nbdel = invtpar_getdouble(m, key, val);
    nbscan = 1;
  }
  else if ( strcmpfuzzy(key, "nbmax") == 0 )
  {
    nbmax = invtpar_getdouble(m, key, val);
    nbscan = 1;
  }
  /* sig-scan parameters */
  else if ( strcmpfuzzy(key, "sigmin") == 0 )
  {
    sigmin = invtpar_getdouble(m, key, val);
    sigscan = 1;
  }
  else if ( strcmpfuzzy(key, "sigdel") == 0
         || strcmpfuzzy(key, "dsig") == 0 )
  {
    sigdel = invtpar_getdouble(m, key, val);
    sigscan = 1;
  }
  else if ( strcmpfuzzy(key, "sigmax") == 0 )
  {
    sigmax = invtpar_getdouble(m, key, val);
    sigscan = 1;
  }
  else if ( strcmpfuzzy(key, "kcscan") == 0
         || strcmpfuzzy(key, "Kscan") == 0 )
  {
    kcscan = 1;
  }
  else if ( strcmpfuzzy(key, "kcdel") == 0
         || strcmpfuzzy(key, "dk") == 0 )
  {
    kcdel = invtpar_getint(m, key, val);
    kcscan = 1;
  }
  else
  {
    return -1;
  }
  return 0;
}



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
  int n = m->n, winn;
  double T, c0, c1, t0, qT, inita, err0, err1, err2;
  double *win, *lambda;
  double *xerri, *xerrf, *xerrf_r, *xerrf_a;
  intq_t *intq;

  xnew(xerri, n);
  xnew(xerrf, n);
  xnew(xerrf_r, n);
  xnew(xerrf_a, n);

  T = (double) m->nsteps;

  xnew(win, n);
  xnew(lambda, n);
  /* prepare the window function and
   * compute the eigenvalues of the updating matrix */
  prepwin(lambda, m->n, m->win, m->winn, win, &winn,
      m->pbc, m->gaussig, m->kc,
      m->fnwin, m->fnwinmat, m->verbose);

  /* initial errors */
  esterror_eql(m->alpha0, m->n, xerri, lambda, gamma);

  c0 = m->c;
  err0 = esterror_invt(T, m->c, m->alpha0, m->t0, m->n,
      NULL, lambda, gamma);

  /* compute the optimal c and the error
   * for the inverse-time formula */
  c1 = estbestc_invt(T, m->alpha0, 0, m->n, lambda, gamma,
      0, &err1, 0);

  /* compute the minimal error from the optimal schedule
   * under the same condition */
  qT = m->qT;
  err2 = esterror_opt(T, m->alpha0, m->initalpha, &qT, m->qprec,
      m->alpha_nint, &intq, m->n, m->errkc, m->pbc,
      lambda, gamma, m->verbose);
  inita = intq_getinita(intq);

  /* save the optimal schedule to file */
  t0 = 2 / m->alpha0;
  intq_save(intq, c1, t0, m->alpha_resample, m->fnalpha);

  printf("c %g, t0 %g, err %g, sqr %g (invt), "
         "opt. c %g, err %g, sqr %g (opt. invt), "
         "qT %g, a(0) %g, err %g, sqr %g (exact), %s\n",
      c0, t0, err0, err0 * err0,
      c1, err1, err1 * err1,
      qT, inita, err2, err2 * err2, m->fnalpha);

  /* print out the mass distribution */
  //intq_getmint(intq, qT, NULL, "mass.dat");

  if ( m->fnxerr[0] != '\0' ) {
    if ( m->opta ) {
      /* optimal schedule */
      intq_errcomp(intq, m->alpha0, qT, xerrf, xerrf_r, xerrf_a);
    } else {
      /* inverse-time schedule */
      err0 = esterror_invt_x(T, m->c, m->alpha0, m->t0, m->n,
          NULL, NULL, xerrf, xerrf_r, xerrf_a, lambda, gamma);
    }
    save_xerr(m, m->fnxerr, xerri, xerrf, xerrf_r, xerrf_a,
        lambda, gamma, qT, inita);
  }

  intq_close(intq);
  free(win);
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
  xnew(lambda, m->n);
  geteigvals(lambda, m->n, m->win, m->winn, m->pbc,
      0, NULL, 1);

  /* initial equilibrium error */
  erri = esterror_eql(m->alpha0, m->n, NULL,
      lambda, gamma);

  /* compute the exact minimal error under the same condition */
  qT = m->qT;
  err2 = esterror_opt(T, m->alpha0, m->initalpha, &qT, m->qprec,
      m->alpha_nint, NULL, m->n, m->errkc, m->pbc,
      lambda, gamma, m->verbose);

  /* print out a header */
  printf("# c     \t  final error\t  init. error\t  final equil\t  optimal error\n");

  for ( c = cmin; c < cmax + 0.001 * cdel; c += cdel ) {
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
  xnew(lambda, m->n);
  geteigvals(lambda, m->n, m->win, m->winn, m->pbc,
      0, NULL, 1);

  /* initial equilibrium error */
  erri = esterror_eql(m->alpha0, m->n, NULL,
      lambda, gamma);

  intq = intq_open(T, m->alpha_nint, m->n,
      m->errkc, m->pbc, lambda, gamma);

  /* print out a header */
  printf("# inita \t  final error\t  init. error\t   qT\n");

  iafac = pow(10.0, iadel);
  for ( inita = iamin; inita < iamax * iafac; inita *= iafac ) {
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
  SCAN_KC
};



static void invt_scan(invtpar_t *m,
    const double *gamma, int scantype)
{
  double nb, sig;
  int kc, kcmax = m->kc, eigerr = 0;
  double T, t0, c, qT, err1, err1norm, err2, err2norm;
  double *lambda;

  T = (double) m->nsteps;
  xnew(lambda, m->n);

  /* initialize the scanning variable */
  if ( scantype == SCAN_NB )
  {
    nb = nbmin;
    printf("# nb    ");
  }
  else if ( scantype == SCAN_SIG )
  {
    sig = sigmin;
    printf("# sigma ");
  }
  else if ( scantype == SCAN_KC )
  {
    kc = kcdel;
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
      geteigvals(lambda, m->n, m->win, m->winn, m->pbc,
          0, &eigerr, 1);
    }
    else if ( scantype == SCAN_SIG )
    {
      /* make the Gaussian window */
      m->gaussig = sig;
      mkgauswin(sig, m->n, m->pbc, m->win, &m->winn);
      stablizewin(lambda, m->n, m->win, &m->winn, m->pbc, 0.0, m->verbose);
    }
    else if ( scantype == SCAN_KC )
    {
      /* make the bandpass sinc window */
      m->kc = kc;
      mksincwin(kc, m->n, m->pbc, m->win, &m->winn);
      geteigvals(lambda, m->n, m->win, m->winn, m->pbc,
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
        m->alpha_nint, NULL, m->n, m->errkc, m->pbc,
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
    else if ( scantype == SCAN_KC )
    {
      printf("%8d\t", kc);
    }

    printf("%10.6f\t%10.6f\t%10.4f\t%10.6f\t%10.4f\t%10g\t%10g\n",
        c, err1, err1norm, err2, err2norm, t0, qT);

    free(lambda);

    /* update the scanning variable and stop */
    if ( scantype == SCAN_NB )
    {
      nb += nbdel;
      if ( nb > nbmax + 0.01 * nbdel ) break;
    }
    else if ( scantype == SCAN_SIG )
    {
      sig += sigdel;
      if ( sig > sigmax + 0.01 * sigdel ) break;
    }
    else if ( scantype == SCAN_KC )
    {
      kc += kcdel;
      if ( kc > kcmax ) break;
    }
  }
}



int main(int argc, char **argv)
{
  invtpar_t m[1];
  double *gamma;

  invtpar_init(m);
  m->userhelp = &predpar_help;
  m->usermatch = &predpar_keymatch;
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);

  /* estimate or load the gamma values */
  xnew(gamma, m->n);
  estgamma(gamma, m->n, m->sampmethod, m->pbc, m->mvsize);
  if ( m->fngamma[0] != '\0' ) {
    if ( m->sampmethod == SAMPMETHOD_MD ) {
      loadgamma(m->n, gamma, m->fngamma);
    } else {
      savegamma(m->n, gamma, m->fngamma);
    }
  }

  if ( !cscan && !iascan && !nbscan && !sigscan && !kcscan ) {
    invt_geterr(m, gamma);
  }

  if ( cscan ) {
    fprintf(stderr, "c-scan window: %g:%g:%g\n", cmin, cdel, cmax);
    invt_scanc(m, gamma);
  }
  if ( iascan ) {
    fprintf(stderr, "ia-scan window: %g:%g:%g\n", iamin, iadel, iamax);
    invt_scania(m, gamma);
  }
  if ( nbscan ) {
    fprintf(stderr, "nb-scan window: %g:%g:%g\n", nbmin, nbdel, nbmax);
    invt_scan(m, gamma, SCAN_NB);
  }
  if ( sigscan ) {
    fprintf(stderr, "sig-scan window: %g:%g:%g\n", sigmin, sigdel, sigmax);
    invt_scan(m, gamma, SCAN_SIG);
  }
  if ( kcscan ) {
    fprintf(stderr, "kc-scan window: 0:%d:%d\n", kcdel, m->kc);
    invt_scan(m, gamma, SCAN_KC);
  }

  free(gamma);

  invtpar_finish(m);

  return 0;
}
