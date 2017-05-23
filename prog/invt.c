/* main driver
 * One-dimensional model system */



#include "invt.h"
#include "cmvar.h" /* eigenmode decomposition */
#include "intq.h"
#include "invtsamp.h" /* MD and MC sampling of the 1D test system */
#include "corr.h"
#include "ave.h"



/* compute the updating magntidue */
static double getalpha(const invtpar_t *m, double t, intq_t *intq)
{
  double a = 0;

  if ( intq != NULL ) {
    a = intq_interpa(intq, t);
    //printf("a %g (a_invt: %g), id %d, t %g, %g\n", a, m->c/(t+m->t0), id, t, intq->tarr[id]);
  } else {
    a = m->c / (t + m->t0);
    //printf("a %g, t %g, t0 %g, c %g\n", a, t, m->t0, m->c); getchar();
  }
  return a;
}



/* multiple-bin update */
__inline static void mbin_update(double *v, int n, int i,
    double a, const double *win, int winn, int pbc)
{
  int j, k;

  v[i] += a * win[0];
  /* update the neighbors */
  for ( j = 1; j < winn; j++ ) {

    k = i - j;
    if ( k < 0 ) {
      k = pbc ? k + n : - k - 1;
    }
    v[k] += a * win[ j ];

    if ( j * 2 == n && pbc ) {
      continue;
    }

    k = i + j;
    if ( k >= n ) {
      k = pbc ? k - n : 2 * n - 1 - k;
    }

    v[k] += a * win[j];
  }
}



/* preparation metadynamics run to compute the gamma values */
__inline static void premeta(const invtpar_t *m, double *gamma)
{
  int i, n = m->n;
  long t;
  double a, *gamma0;

  /* data for sampling */
  invtsamp_t *its;

  /* accumulator for computing the variance of the bias potential */
  cmvar_t *cm;

  /* open an object for sampling */
  its = invtsamp_open(m);

  cm = cmvar_open(n, m->pbc);

  xnew(gamma0, n);

  for ( i = 0; i < n; i++ ) {
    gamma0[i] = gamma[i];
    gamma[i] = 0;
  }

  i = 0;
  for ( t = 0; t < m->gam_nsteps; t++ ) {
    /* sampling */
    if ( m->tcorr <= 0 || rand01() * m->tcorr < 1 ) {
      invtsamp_step(its, &i);
    }

    /* compute the updating magnitude */
    a = m->alpha0 / m->p[i];

    /* update the bias potential (WL) way */
    its->v[i] += a;

    /* accumulate data for variance */
    if ( (t + 1) % m->gam_nstave == 0 ) {
      cmvar_add(cm, its->v);
    }
  }

  /* compute the integrals of the autocorrelation functions */
  cmvar_get(cm);
  fprintf(stderr, "mode    old-gamma    new-gamma\n");
  for ( i = 1; i < n; i++ ) {
    gamma[i] = 2 * cm->uvar[i] / m->alpha0;
    fprintf(stderr, "%4d %12.6f %12.6f\n", i, gamma0[i], gamma[i]);
  }

  if ( m->fngamma[0] != '\0' )
    savegamma(n, gamma, m->fngamma);

  invtsamp_close(its);
}



/* simulate a metadynamics process
 * return the root-mean-squared error of the inverse time scheme */
static double simulmeta(const invtpar_t *m, intq_t *intq,
    const double *win, int winn,
    double *err0, const double *xerr0)
{
  double a, err;
  int i, n = m->n, prod = 0;
  long t;

  /* data for sampling */
  invtsamp_t *its;

  /* data for computing correlation functions */
  corr_t *corr = NULL;
  int nfrcorr = 0;


  /* open an object for sampling */
  its = invtsamp_open(m);

  /* open an object for correlation functions */
  if ( m->docorr ) {
    nfrcorr = m->nsteps / m->nstcorr + 1;
    corr = corr_open(n - 1, nfrcorr);
  }

  i = 0;
  for ( t = 0; t < m->nsteps + m->nequil; t++ ) {
    /* sampling */
    if ( m->tcorr <= 0 || rand01() * m->tcorr < 1 ) {
      invtsamp_step(its, &i);
    }

    /* try to turn on production after equilibration */
    if ( !prod && t >= m->nequil ) {
      prod = 1;
      /* compute the initial error */
      *err0 = geterror(its->v, n, m->p);
      if ( m->verbose >= 1 ) {
        fprintf(stderr, "starting production at step %ld, t0 %g, err %g\n",
            t, m->t0, *err0);
      }
    }

    /* compute the updating magnitude */
    if ( !prod || m->fixa ) {
      a = m->alpha0;
    } else {
      a = getalpha(m, (double) (t - m->nequil), intq);
    }
    a /= m->p[i];

    /* update the bias potential  */
    mbin_update(its->v, n, i, a, win, winn, m->pbc);

    /* accumulate data for correlation functions */
    if ( prod && corr != NULL && (t + 1) % m->nstcorr == 0 ) {
      shift(its->v, n);
      /* Fourier transform to get the modes `u` */
      getcosmodes(its->v, n, its->u, its->costab);
      //getcosmodesh(i, n, its->u, its->costab);
      /* the first mode should always be zero,
       * so we start from the second mode, u + 1 */
      corr_add(corr, its->u + 1);
    }
  }

  /* compute the error */
  err = geterror(its->v, n, m->p);

  /* print out the error of the modes */
  if ( m->verbose >= 2 ) {
    getcosmodes(its->v, n, its->u, its->costab);
    for ( i = 0; i < n; i++ ) {
      printf("  %+11.8f", its->u[i]);
    }
    printf("\n");
  }

  if ( corr != NULL ) {
    /* print out the thermodynamic fluctuations */
    corr_printfluc(corr, 0, xerr0 + 1);

    /* compute the correlation functions and save them to file
     * the maximal span, 100/alpha, should be large enough */
    corr_save(corr, m->nstcorr, 100 / m->alpha0,
        m->corrtol, 0, m->fncorr);
  }

  invtsamp_close(its);

  if ( corr != NULL ) {
    corr_close(corr);
  }

  return err;
}



typedef struct {
  int n;
  double T;
  int winn;
  double *win;
  double *lambda;
  double *gamma;
  intq_t *intq;
  double *xerr0, *xerr;
  double err0ref, errref, err1ref; /* initial, final, final saturated */
} invtdata_t;


static invtdata_t *invt_open(invtpar_t *m)
{
  invtdata_t *invt;

  xnew(invt, 1);
  invt->n = m->n;
  invt->T = (double) m->nsteps;
  xnew(invt->win, invt->n);
  xnew(invt->lambda, invt->n);
  xnew(invt->gamma, invt->n);
  invt->intq = NULL;

  /* prepare the window function, compute lambda's */
  prepwin(invt->lambda, invt->n, m->win, m->winn,
      invt->win, &invt->winn, m->pbc, m->gaussig, m->kc,
      m->fnwin, m->fnwinmat, m->verbose);

  /* estimate the integrals of the autocorrelation functions
   * of the eigenmodes for the updating scheme */
  estgamma(invt->gamma, m->n, m->sampmethod, m->pbc, m->mvsize);

  xnew(invt->xerr0, invt->n);
  xnew(invt->xerr, invt->n);
  invt->errref = 0;
  invt->err1ref = 0;

  /* theoretical estimate of the initial saturated error */
  invt->err0ref = esterror_eql(m->alpha0, m->n, invt->xerr0, invt->lambda, invt->gamma);

  return invt;
}

static void invt_close(invtdata_t *invt)
{
  if ( invt->intq != NULL ) {
    intq_close( invt->intq );
  }
  free(invt->win);
  free(invt->lambda);
  free(invt->gamma);
  free(invt->xerr0);
  free(invt->xerr);
  free(invt);
}

/* prepare the schedule */
static void invt_getalpha(invtdata_t *invt, invtpar_t *m)
{
  double optc;
  double alphaf; /* final updating magnitude */
  double errmin = 0; /* predicted minimal error */
  /* reference values */

  if ( m->fixa ) return;

  /* compute the theoretically optimal c */
  optc = estbestc_invt(invt->T, m->alpha0, 0, m->n, invt->lambda, invt->gamma,
      0, &errmin, m->verbose - 1);

  /* use the theoretically optimal c */
  if ( m->optc ) {
    m->c = optc;
    fprintf(stderr, "use the optimal c = %g\n", m->c);
  }

  /* don't be smart about t0, we want a constant t0 as
   * a part of the normalization factor
   * even if the schedule its inverse-time, don't multiply c */
  m->t0 = 2 / m->alpha0;

  if ( m->opta ) {
    /* compute the theoretically optimal schedule */
    invt->errref = esterror_opt(invt->T, m->alpha0, m->initalpha, &m->qT, m->qprec,
        m->alpha_nint, &invt->intq, m->n, m->errkc, m->pbc,
        invt->lambda, invt->gamma, m->verbose);
    errmin = invt->errref;

    //m->t0 = intq_estt0(invt->T, m->qT);
    //if ( m->initalpha > 0 ) m->t0 = 1 / m->initalpha;
    fprintf(stderr, "q(T) %g, t0 = %g\n", m->qT, m->t0);

    /* save the optimal schedule to file */
    intq_save(invt->intq, optc, m->t0,
        m->alpha_resample, m->fnalpha);

    alphaf = invt->intq->aarr[invt->intq->m - 1];
  } else {
    /* compute the optimal error from the inverse-time formula */
    invt->errref = esterror_invt(invt->T, m->c, m->alpha0, m->t0, m->n,
        invt->xerr, invt->lambda, invt->gamma);

    printf("predicted optimal c %g, err %g, sqr: %g\n",
        optc, errmin, errmin * errmin);

    alphaf = m->c / (invt->T + m->t0);
  }

  /* theoretical estimate of the final saturated error */
  invt->err1ref = esterror_eql(alphaf, m->n, NULL, invt->lambda, invt->gamma);
  printf("estimated final saturated error %g, sqr: %g\n",
      invt->err1ref, invt->err1ref * invt->err1ref);

  printf("estimated final error %g, sqr: %g, norm. sqr: %g\n",
      invt->errref, invt->errref * invt->errref,
      invt->errref * invt->errref * (invt->T + m->t0));

  if ( m->verbose ) {
    if ( m->opta ) {
      intq_errcomp(invt->intq, m->alpha0, m->qT, invt->xerr, NULL, NULL);
    }
    dumperror(m->n, invt->lambda, invt->gamma, 2, invt->xerr0, invt->xerr);
  }
}


static double invt_run(invtpar_t *m)
{
  invtdata_t *invt;
  double err = 0, averr, stde;  /* final */
  double err0, averr0, stde0; /* initial */
  ave_t ei[1], ef[1];
  long ntr = m->ntrials;
  long i;

  /* clock() its probably better than time(NULL) */
  mtscramble( clock() );

  invt = invt_open(m);

  ave_clear(ei);
  ave_clear(ef);

  /* do a trial run to compute the gamma values */
  if ( m->gammethod == GAMMETHOD_LOAD ) {
    loadgamma(invt->n, invt->gamma, m->fngamma);
  } else if ( m->gammethod != GAMMETHOD_NONE ) {
    premeta(m, invt->gamma);
  }

  if ( m->docorr ) {
    /* compute the correlation functions for a single run */
    err = simulmeta(m, NULL, invt->win, invt->winn, &err0, invt->xerr0);
  } else {
    /* do multiple runs to compute the average error */
    printf("estimated initial error %g, sqr: %g\n",
        invt->err0ref, invt->err0ref * invt->err0ref);

    /* prepare the schedule */
    invt_getalpha(invt, m);

    /* repeat `ntr` independent simulations */
    for ( i = 0; i < ntr; i++ ) {
      /* run simulation */
      err = simulmeta(m, invt->intq, invt->win, invt->winn, &err0, invt->xerr0);

      /* accumulators for the initial and final errors */
      ave_add(ei, err0 * err0);
      averr0 = sqrt( ei->ave );
      ave_add(ef, err * err);
      averr  = sqrt( ef->ave );

      printf("%4ld: err %10.8f -> %10.8f, ave %10.8f -> %10.8f, sqr %e -> %e\n",
          i, err0, err, averr0, averr, ei->ave, ef->ave);
    }

    /* statistics for the final square error */
    averr = sqrt( ef->ave );
    stde = sqrt( ef->var );
    if ( ntr > 1 ) stde /= sqrt(ntr - 1.0) ;

    /* statistics for the initial square error */
    averr0 = sqrt( ei->ave );
    stde0 = sqrt( ei->var );
    if ( ntr > 1 ) stde0 /= sqrt(ntr - 1.0);

    printf("average error: %10.8f -> %10.8f, sqr %e -> %e, "
        "norm. sqr %e, stdsqr %e -> %e\n",
        averr0, averr, ei->ave, ef->ave,
        ef->ave * (invt->T + m->t0), stde0, stde);
    printf("predicted val: %10.8f -> %10.8f, sqr %e -> %e, "
        "norm. sqr %e\n",
        invt->err0ref, invt->errref, invt->err0ref * invt->err0ref, invt->errref * invt->errref,
        invt->errref * invt->errref * (invt->T + m->t0));
    printf("saturated val: %10.8f -> %10.8f, sqr %e -> %e\n",
        invt->err0ref, invt->err1ref, invt->err0ref * invt->err0ref, invt->err1ref * invt->err1ref);
  }

  invt_close(invt);
  return err;
}



int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  invt_run(m);
  invtpar_finish(m);

  return 0;
}
