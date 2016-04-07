/* main driver */



#include "invt.h"
#include "cosmodes.h"
#include "intq.h"
#include "invtsamp.h"
#include "corr.h"



/* compute the updating magntidue */
static double getalpha(const invtpar_t *m, double t, intq_t *intq)
{
  double a = 0;
  static int id = 0;

  if ( intq != NULL ) {
    a = intq_interpa(intq, t, &id);
    //printf("a %g (a_invt: %g), id %d, t %g, %g\n", a, m->c/(t+m->t0), id, t, intq->tarr[id]);
  } else {
    a = m->c / (t + m->t0);
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

    v[k] += a * win[ abs(j) ];
  }
}



/* preparation metadynamics run to compute the gamma values */
__inline static void premeta(const invtpar_t *m, double *gamma)
{
  int i, j, n = m->n;
  long t;
  double a, *gamma0;
  double ucnt = 0, *uave, *uvar;

  /* data for sampling */
  invtsamp_t *is;

  /* open an object for sampling */
  is = invtsamp_open(m);

  xnew(gamma0, n);
  xnew(uave, n);
  xnew(uvar, n);

  for ( i = 0; i < n; i++ ) {
    gamma0[i] = gamma[i];
    gamma[i] = 0;
    uave[i] = 0;
    uvar[i] = 0;
  }

  i = 0;
  for ( t = 0; t < m->gam_nsteps; t++ ) {
    /* sampling */
    if ( m->tcorr <= 0 || rand01() * m->tcorr < 1 ) {
      invtsamp_step(is, &i);
    }

    /* compute the updating magnitude */
    a = m->alpha0 / m->p[i];

    /* update the bias potential (WL) way */
    is->v[i] += a;

    /* accumulate data for variance */
    if ( (t + 1) % m->gam_nstave == 0 ) {
      shift(is->v, n);
      /* Fourier transform to get the modes `u` */
      getcosmodes(is->v, n, is->u, is->costab);
      /* accumulate data for average and variance */
      for ( j = 0; j < n; j++ ) {
        double y = 0, du;
        if ( ucnt > 0 ) {
          y = uave[j] / ucnt;
        }
        uave[j] += is->u[j]; /* update the average */
        if ( ucnt > 0 ) { /* update the variance */
          du = is->u[j] - y;
          uvar[j] += du * du * ucnt / (ucnt + 1);
        }
      }
      ucnt += 1;
    }
  }

  /* compute the integrals of the autocorrelation functions */
  fprintf(stderr, "mode    old-gamma    new-gamma\n");
  for ( i = 1; i < n; i++ ) {
    uave[i] /= ucnt;
    uvar[i] /= ucnt;
    gamma[i] = 2 * uvar[i] / m->alpha0;
    fprintf(stderr, "%4d %12.6f %12.6f\n", i, gamma0[i], gamma[i]);
  }

  if ( m->pregamma == 1 && m->fngamma[0] != '\0' ) {
    savegamma(n, gamma, m->fngamma);
  }

  invtsamp_close(is);
}



/* simulate a metadynamics process
 * return the root-mean-squared error of the inverse time scheme */
static double simulmeta(const invtpar_t *m, intq_t *intq,
    double *err0, const double *xerr0)
{
  double a, err;
  int i, n = m->n, prod = 0;
  long t;

  /* data for sampling */
  invtsamp_t *is;

  /* data for computing correlation functions */
  corr_t *corr = NULL;
  int nfrcorr = 0;


  /* open an object for sampling */
  is = invtsamp_open(m);

  /* open an object for correlation functions */
  if ( m->docorr ) {
    nfrcorr = m->nsteps / m->nstcorr + 1;
    corr = corr_open(n - 1, nfrcorr);
  }

  i = 0;
  for ( t = 0; t < m->nsteps + m->nequil; t++ ) {
    /* sampling */
    if ( m->tcorr <= 0 || rand01() * m->tcorr < 1 ) {
      invtsamp_step(is, &i);
    }

    /* try to turn on production after equilibration */
    if ( !prod && t >= m->nequil ) {
      prod = 1;
      /* compute the initial error */
      *err0 = geterror(is->v, n, m->p);
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
    mbin_update(is->v, n, i, a, m->win, m->winn, m->pbc);

    /* accumulate data for correlation functions */
    if ( prod && corr != NULL && (t + 1) % m->nstcorr == 0 ) {
      shift(is->v, n);
      /* Fourier transform to get the modes `u` */
      getcosmodes(is->v, n, is->u, is->costab);
      /* the first mode is always zero,
       * so we start from the second mode, u + 1 */
      corr_add(corr, is->u + 1);
    }
  }

  /* compute the error */
  err = geterror(is->v, n, m->p);

  /* print out the error of the modes */
  if ( m->verbose >= 2 ) {
    getcosmodes(is->v, n, is->u, is->costab);
    for ( i = 0; i < n; i++ ) {
      printf("  %+11.8f", is->u[i]);
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

  invtsamp_close(is);

  if ( corr != NULL ) {
    corr_close(corr);
  }

  return err;
}



static double invt_run(invtpar_t *m)
{
  double err = 0,  e,  se = 0,  see = 0,  ave,  averr,  stde;  /* final */
  double err0, e0, se0 = 0, see0 = 0, ave0, averr0, stde0; /* initial */
  /* reference values */
  double errref = 0, err0ref, err1ref = 0; /* final, initial, final saturated */
  double optc, errmin = 0; /* optimal c, predicted minimal error */
  double T;
  double alphaf; /* final updating magnitude */
  double *lambda = NULL, *gamma = NULL;
  double *xerr0 = NULL, *xerr = NULL;
  intq_t *intq = NULL;
  long ntr = m->ntrials;
  long i;

  /* clock() is probably better than time(NULL) */
  mtscramble( clock() );

  /* estimate the eigenvalues of the w matrix,
   * for the updating scheme */
  //lambda = geteigvals(m->n, m->winn, m->win, m->pbc, 0, NULL, 1);
  /* modify the window function such that all eigenvalues
   * lambda[i] are positive-definite */
  lambda = trimwindow(m->n, &m->winn, m->win, m->pbc, 0, m->verbose);
  if ( m->fnwin[0] != '\0' ) {
    savewin(m->winn, m->win, m->fnwin);
  }
  if ( m->fnwinmat[0] != '\0' ) {
    savewinmat(m->winn, m->win, m->n, m->pbc, m->fnwinmat);
  }

  /* estimate the integrals of correlations
   * of the eigenmodes for the updating scheme */
  gamma = estgamma(m->n, m->sampmethod, m->pbc, m->localg);

  xnew(xerr0, m->n);
  xnew(xerr, m->n);

  /* theoretical estimate of the initial saturated error */
  err0ref = esterror_eql(m->alpha0, m->n, xerr0, lambda, gamma);

  /* do a trial run to compute the gamma values */
  if ( m->pregamma ) {
    premeta(m, gamma);
  }

  if ( m->docorr ) {
    /* compute the correlation functions for a single run */
    err = simulmeta(m, NULL, &err0, xerr0);
  } else {
    /* do multiple runs to compute the average error */
    T = (double) m->nsteps;

    printf("estimated initial error %g, sqr: %g\n",
        err0ref, err0ref * err0ref);

    if ( !m->fixa ) {
      /* compute the theoretically optimal c */
      optc = estbestc_invt(T, m->alpha0, m->n, lambda, gamma,
          0, &errmin, m->verbose - 1);

      /* use the theoretically optimal c */
      if ( m->optc ) {
        m->c = optc;
        fprintf(stderr, "use the optimal c = %g\n", m->c);
      }

      if ( m->opta ) {
        double qT = 0;

        /* compute the theoretically optimal schedule */
        errref = esterror_opt(T, m->alpha0, &qT, m->qprec,
            m->alpha_nint, &intq, m->n, xerr,
            lambda, gamma, m->verbose);
        errmin = errref;

        m->t0 = intq_estt0(T, qT);
        fprintf(stderr, "q(T) %g, estimated t0 = %g\n", qT, m->t0);

        /* save the optimal schedule to file */
        intq_save(intq, optc, 2 * optc / m->alpha0, m->fnalpha);

        alphaf = intq->aarr[intq->m - 1];
      } else {
        /* compute the optimal error from the inverse-time formula */
        errref = esterror_invt(T, m->c, m->alpha0, m->n, xerr,
            lambda, gamma);

        printf("predicted optimal c %g, err %g, sqr: %g\n",
            optc, errmin, errmin * errmin);

        alphaf = m->c / (T + m->t0);
      }

      /* theoretical estimate of the final saturated error */
      err1ref = esterror_eql(alphaf, m->n, NULL, lambda, gamma);
      printf("estimated final saturated error %g, sqr: %g\n",
          err1ref, err1ref * err1ref);

      printf("estimated final error %g, sqr: %g, norm. sqr: %g\n",
          errref, errref * errref, errref * errref * (T + m->t0));

      if ( m->verbose ) {
        dumperror(m->n, lambda, gamma, 2, xerr0, xerr);
      }
    }

    /* repeat `ntr` independent simulations */
    for ( i = 0; i < ntr; i++ ) {
      /* run simulation */
      err = simulmeta(m, intq, &err0, xerr0);

      /* accumulators for the final error */
      e = err * err;
      se += e;
      see += e * e;
      ave = se / (i + 1);
      averr = sqrt( ave );

      /* accumulators for the initial error */
      e0 = err0 * err0;
      se0 += e0;
      see0 += e0 * e0;
      ave0 = se0 / (i + 1);
      averr0 = sqrt( ave0 );

      printf("%4ld: err %10.8f -> %10.8f, "
                   "ave %10.8f -> %10.8f, "
                   "sqr %e -> %e\n",
          i, err0, err,
          averr0, averr,
          ave0, ave);
    }

    /* statistics for the final square error */
    ave = se / ntr;
    averr = sqrt( ave );
    stde = sqrt( see / ntr - ave * ave );
    if ( ntr > 1 ) {
      stde /= sqrt(ntr - 1.0) ;
    }

    /* statistics for the initial square error */
    ave0 = se0 / ntr;
    averr0 = sqrt( ave0 );
    stde0 = sqrt( see0 / ntr - ave0 * ave0 );
    if ( ntr > 1 ) {
      stde0 /= sqrt(ntr - 1.0);
    }

    printf("average error: %10.8f -> %10.8f, sqr %e -> %e, "
        "norm. sqr %e, stdsqr %e -> %e\n",
        averr0, averr, ave0, ave,
        ave * (T + m->t0), stde0, stde);
    printf("predicted val: %10.8f -> %10.8f, sqr %e -> %e, "
        "norm. sqr %e\n",
        err0ref, errref, err0ref * err0ref, errref * errref,
        errref * errref * (T + m->t0));
    printf("saturated val: %10.8f -> %10.8f, sqr %e -> %e\n",
        err0ref, err1ref, err0ref * err0ref, err1ref * err1ref);
  }

  if ( intq != NULL ) {
    intq_close( intq );
  }
  free(lambda);
  free(gamma);
  free(xerr0);
  free(xerr);

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
