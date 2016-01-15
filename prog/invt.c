

#include "invt.h"
#include "corr.h"



/* return the root-mean-squared error of the inverse time scheme */
static double comperr(const invtpar_t *m, double *err0)
{
  double *v = NULL, *vac = NULL, a, err;
  int i, n = m->n, prod = 0;
  long t;

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
  normalize(v, n, m->initrand);

  /* space for the accumulative distribution function */
  xnew(vac, n + 1);

  /* open an object for correlation functions */
  if ( m->nstcorr > 0 ) {
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
        i = mc_metro_l(v, n, i);
      } else if ( m->sampmethod == SAMPMETHOD_HEATBATH ) {
        i = mc_heatbath(v, vac, n);
      }
    }

    /* try to start production */
    if ( !prod && t >= m->nequil ) {
      prod = 1;
      *err0 = geterror(v, n);
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

    if ( m->nbn > 1 ) {
      mbin_update(v, n, i, a, m->nbs, m->nbn);
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
  err = geterror(v, n);

  /* print out the error */
  if ( m->verbose >= 2 ) {
    for ( i = 0; i < n; i++ ) {
      printf("  %+11.8f", v[i]);
    }
    printf("\n");
  }

  if ( corr != NULL ) {
    /* save the correlation functions to file */
    corr_save(corr, m->nstcorr, m->corrtol, m->fncorr);
  }

  free(v);
  free(vac);
  free(u);
  free(costab);

  if ( corr != NULL ) {
    corr_close(corr);
  }

  return err;
}



static double invt_run(invtpar_t *m)
{
  double err, see = 0;
  double err0, see0 = 0;
  int i;

  mtscramble( time(NULL) );

  if ( m->fixa ) {
    err = comperr(m, &err0);
  } else {
    for ( i = 0; i < m->ntrials; i++ ) {
      err = comperr(m, &err0);
      see += err * err;
      see0 += err0 * err0;
      printf("%4d: err %10.8f -> %10.8f, "
                  "ave %10.8f -> %10.8f, "
                  "sqr %e -> %e\n",
          i, err0, err,
          sqrt(see0/(i+1)), sqrt(see/(i+1)),
          see0/(i+1), see/(i+1));
    }

    err = see / m->ntrials;
    err0 = see0 / m->ntrials;
    printf("average error: %10.8f -> %10.8f, sqr %e -> %e\n",
        sqrt(err0), sqrt(err), err0, err);
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
