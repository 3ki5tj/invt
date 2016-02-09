/* main driver */



#include "invt.h"
#include "invtsamp.h"
#include "corr.h"
#include "cosmodes.h"
#include "intq.h"



/* compute the updating magntidue */
static double getalpha(const invtpar_t *m, double t,
    intq_t *intq, int fixa)
{
  double a = 0;
  static int id = 0;

  if ( fixa || m->fixa ) {
    return m->alpha0;
  }

  if ( intq != NULL ) {
    a = intq_interpa(intq, t, &id);
    //printf("a %g (a_invt: %g), id %d, t %g, %g\n", a, m->c/(t+m->t0), id, t, intq->tarr[id]);
  } else {
    a = m->c / (t + m->t0);
  }
  return a;
}



/* multiple-bin update */
static void mbin_update(double *v, int n, int i,
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

    k = i + j;
    if ( k >= n ) {
      k = pbc ? k - n : 2 * n - 1 - k;
    }

    v[k] += a * win[ abs(j) ];
  }
}



/* simulate a metadynamics process
 * return the root-mean-squared error of the inverse time scheme */
static double simulmeta(const invtpar_t *m, intq_t *intq,
    double *err0, const double *xerr0)
{
  double *v = NULL, a, err;
  int i, n = m->n, prod = 0;
  long t;

  /* data for sampling method */
  double *vac = NULL; /* for the heatbath algorithm */
  ouproc_t *ou = NULL; /* for Ornstein-Uhlenbeck process */
  invtmd_t invtmd[1]; /* for molecular dynamics */

  /* data for computing correlation functions */
  corr_t *corr = NULL;
  double *u = NULL, *costab = NULL;
  int nfrcorr = 0;


  xnew(v, n);
  xnew(u, n);

  /* compute the coefficients for Fourier transform
   * used in decomposition of modes */
  costab = mkcostab(n);

  /* initially randomize the error */
  /* initialize the magnitude of the Fourier modes
   * up to the cutoff wave number */
  for ( i = 0; i < m->kcutoff; i++ ) {
    u[i] = randgaus();
  }
  /* combine the modes */
  fromcosmodes(v, n, u, costab);
  normalize(v, n, m->initrand, m->p);

  /* initialize data for the sampling methods */
  if ( m->sampmethod == SAMPMETHOD_HEATBATH ) {
    /* allocated space for the accumulative distribution function
     * used for heatbath algorithm */
    xnew(vac, n + 1);
  } else if ( m->sampmethod == SAMPMETHOD_OU ) {
    ou = ouproc_open(v, n, 1.0 / ( m->tcorr + 1e-16 ) );
  } else if ( m->sampmethod == SAMPMETHOD_MD ) {
    invtmd_init(invtmd, n,
        m->mddt, m->tp, m->thermdt, v);
  }

  /* open an object for correlation functions */
  if ( m->docorr ) {
    nfrcorr = m->nsteps / m->nstcorr + 1;
    corr = corr_open(n - 1, nfrcorr);
  }

  i = 0;
  for ( t = 0; t < m->nsteps + m->nequil; t++ ) {
    /* sampling */
    if ( m->tcorr <= 0 || rand01() * m->tcorr < 1 )
    {
      if ( m->sampmethod == SAMPMETHOD_METROGLOBAL ) {
        i = mc_metro_g(v, n, i);
      } else if ( m->sampmethod == SAMPMETHOD_METROLOCAL ) {
        i = mc_metro_l(v, n, i, m->pbc);
      } else if ( m->sampmethod == SAMPMETHOD_HEATBATH ) {
        i = mc_heatbath(v, vac, n);
      } else if ( m->sampmethod == SAMPMETHOD_OU ) {
        i = ouproc_step(ou);
      } else if ( m->sampmethod == SAMPMETHOD_MD ) {
        i = invtmd_vv(invtmd);
      }
    }

    /* try to turn on production after equilibration */
    if ( !prod && t >= m->nequil ) {
      prod = 1;
      /* compute the initial error */
      *err0 = geterror(v, n, m->p);
      if ( m->verbose >= 1 ) {
        fprintf(stderr, "starting production at step %ld, t0 %g, err %g\n",
            t, m->t0, *err0);
      }
    }

    /* compute the updating magnitude */
    a = getalpha(m, (double) (t - m->nequil), intq, !prod);
    a /= m->p[i];

    /* update the bias potential */
    mbin_update(v, n, i, a, m->win, m->winn, m->pbc);

    /* accumulate data for correlation functions */
    if ( prod && corr != NULL && (t + 1) % m->nstcorr == 0 ) {
      shift(v, n);
      /* Fourier transform to get the modes `u` */
      getcosmodes(v, n, u, costab);
      /* the first mode is always zero,
       * so we start from the second mode, u + 1 */
      corr_add(corr, u + 1);
    }
  }

  /* compute the error */
  err = geterror(v, n, m->p);

  /* print out the error of the modes */
  if ( m->verbose >= 2 ) {
    getcosmodes(v, n, u, costab);
    for ( i = 0; i < n; i++ ) {
      printf("  %+11.8f", u[i]);
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

  free(v);
  free(vac);
  free(u);
  free(costab);

  if ( m->sampmethod == SAMPMETHOD_OU ) {
    ouproc_close(ou);
  }

  if ( corr != NULL ) {
    corr_close(corr);
  }

  return err;
}



static double invt_run(invtpar_t *m)
{
  double err,  e,  se = 0,  see = 0,  ave,  averr,  stde;  /* final */
  double err0, e0, se0 = 0, see0 = 0, ave0, averr0, stde0; /* initial */
  /* reference values */
  double errref, err0ref, err1ref = 0; /* final, initial, final saturated */
  double optc, errmin = 0; /* optimal c, predicted minimal error */
  double t;
  double *lambda = NULL, *gamma = NULL;
  double *xerr0 = NULL, *xerr = NULL;
  intq_t *intq = NULL;
  long ntr = m->ntrials;
  long i;

  /* clock() is probably better than time(NULL) */
  mtscramble( clock() );

  /* estimate the eigenvalues of the w matrix,
   * for the updating scheme */
  //lambda = geteigvals(m->n, m->winn, m->win, NULL, 1);
  lambda = trimwindow(m->n, &m->winn, m->win, 0);

  /* estimate the correlation integrals
   * of the eigenmodes of the w matrix,
   * for the updating scheme */
  gamma = estgamma(m->n, m->sampmethod);

  xnew(xerr0, m->n);
  xnew(xerr, m->n);

  /* initial saturated error */
  err0ref = esterror_eql(m->alpha0, m->n, xerr0, lambda, gamma);

  if ( m->docorr ) {
    /* compute correlation functions
     * for a single run */
    err = simulmeta(m, NULL, &err0, xerr0);
  } else {
    /* do multiple runs to compute the average error */
    t = (double) m->nsteps;

    printf("estimated initial error %g, sqr: %g\n",
        err0ref, err0ref * err0ref);

    if ( !m->fixa ) {
      err1ref = esterror_eql(m->c / (m->t0 + t), m->n, NULL,
          lambda, gamma);

      printf("estimated final saturated error %g, sqr: %g\n",
          err1ref, err1ref * err1ref);

      /* compute the optimal c */
      optc = estbestc_invt(t, m->alpha0, m->n, lambda, gamma,
          0, &errmin, m->verbose - 1);

      if ( m->opta ) {
        /* compute the optimal schedule */
        errref = esterror_opt(t,
            intq_getqt(t, m->c, m->t0), m->alpha0,
            m->alpha_nint, &intq, m->n, xerr,
            lambda, gamma);

        /* save the optimal schedule to file */
        intq_save(intq, m->c, m->c / m->alpha0, m->fnalpha);

        /* compute the optimal schedule for the optimal c */
        errmin = esterror_opt(t,
            intq_getqt(t, optc, optc / m->alpha0), m->alpha0,
            m->alpha_nint, NULL, m->n, NULL,
            lambda, gamma);
      } else {
        /* compute the optimal error from the inverse-time formula */
        errref = esterror_invt(t, m->c, m->alpha0, m->n, xerr,
            lambda, gamma);
      }

      printf("estimated final error %g, sqr: %g\n",
          errref, errref * errref);

      printf("predicted optimal c %g, err %g, sqr: %g\n", optc, errmin, errmin * errmin);

      if ( m->verbose ) {
        dumperror(m->n, lambda, gamma,
            2, xerr0, xerr);
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

    /* statistics for the final error */
    ave = se / ntr;
    averr = sqrt( ave );
    stde = sqrt( see / ntr - ave * ave );
    if ( ntr > 1 ) {
      stde /= sqrt(ntr - 1.0) ;
    }

    /* statistics for the initial error */
    ave0 = se0 / ntr;
    averr0 = sqrt( ave0 );
    stde0 = sqrt( see0 / ntr - ave0 * ave0 );
    if ( ntr > 1 ) {
      stde0 /= sqrt(ntr - 1.0);
    }

    printf("average error: %10.8f -> %10.8f, sqr %e -> %e, stdsqr %e -> %e\n",
        averr0, averr, ave0, ave, stde0, stde);
    printf("predicted val: %10.8f -> %10.8f, sqr %e -> %e\n",
        err0ref, errref, err0ref * err0ref, errref * errref);
    printf("saturated val: %10.8f -> %10.8f, sqr %e -> %e\n",
        err0ref, err1ref, err0ref * err0ref, err1ref * err1ref);
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
