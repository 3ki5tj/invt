/* main driver */



#include "invt.h"
#include "invtsamp.h"
#include "corr.h"



/* simulate a metadynamics process
 * return the root-mean-squared error of the inverse time scheme */
static double simulmeta(const invtpar_t *m, double *err0)
{
  double *v = NULL, *vac = NULL, a, err;
  int i, n = m->n, prod = 0;
  long t;

  ouproc_t *ou = NULL;
  invtmd_t invtmd[1];

  corr_t *corr = NULL;
  double *u = NULL, *costab = NULL;
  int nfrcorr = 0;


  xnew(v, n);
  xnew(u, n);
  costab = makecostab(n);

  /* initially randomize the error */
  for ( i = 0; i < m->kcutoff; i++ ) {
    u[i] = randgaus();
  }
  fromcosmodes(v, n, u, costab);
  normalize(v, n, m->initrand, m->p);

  if ( m->sampmethod == SAMPMETHOD_OU ) {
    ou = ouproc_open(v, n, m->tcorr);
  } else if ( m->sampmethod == SAMPMETHOD_MD ) {
    invtmd_init(invtmd, n,
        m->mddt, m->tp, m->thermdt, v);
  }

  /* space for the accumulative distribution function */
  xnew(vac, n + 1);

  /* open an object for correlation functions */
  if ( m->docorr ) {
    nfrcorr = m->nsteps / m->nstcorr + 1;
    corr = corr_open(n - 1, nfrcorr);
  }

  i = 0;
  for ( t = 0; t < m->nsteps + m->nequil; t++ ) {
    /* MCMC sampling */
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

    /* try to start production */
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
    if ( prod && !m->fixa ) {
      a = m->c / (t - m->nequil + m->t0);
    } else {
      a = m->alpha0;
    }

    a /= m->p[i];

    if ( m->winn > 1 ) {
      mbin_update(v, n, i, a, m->win, m->winn, m->pbc);
    } else {
      v[i] += a;
    }

    if ( prod && corr != NULL && (t + 1) % m->nstcorr == 0 ) {
      shift(v, n);
      getcosmodes(v, n, u, costab);
      /* the first mode is always zero,
       * so we start from u + 1 */
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
  long ntr = m->ntrials;
  long i;

  /* clock() is probably better than time(NULL) */
  //mtscramble( clock() );

  if ( m->docorr ) {
    /* compute correlation functions
     * for a single run */
    err = simulmeta(m, &err0);
  } else {
    /* do multiple runs to compute the average error */
    double t = (double) m->nsteps;

    /* compute the prediction from the analytical result */
    /* initial saturated error */
    err0ref = esterror0_ez(m->alpha0,
        m->n, m->winn, m->win, m->sampmethod,
        "initial", m->verbose + 1);

    if ( !m->fixa ) {
      err1ref = esterror0_ez(m->c / (m->t0 + t),
          m->n, m->winn, m->win, m->sampmethod,
          "final", m->verbose + 1);

      errref = esterror_ez(m->c, t, m->t0,
          m->n, m->winn, m->win, m->sampmethod,
          m->verbose + 1);

      /* compute the optimal c */
      optc = estbestc(t, m->t0,
          m->n, m->winn, m->win, m->sampmethod,
          1e-8, &errmin, m->verbose);
      printf("predicted optimal c %g, err %g\n", optc, errmin);
    }

    for ( i = 0; i < ntr; i++ ) {
      err = simulmeta(m, &err0);

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
