

#include "invt.h"



/* return the root squared error of the inverse time scheme
 * */
static double comperr(const invtpar_t *m, double *err0)
{
  double *v, *vac, a, err;
  int i, n = m->n, prod = 0;
  long t;

  xnew(v, n);
  for ( i = 0; i < n; i++ ) {
    v[i] = 0;
  }

  /* space for the accumulative distribution function */
  xnew(vac, n + 1);

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
    if ( !prod && t >= m->nequil + m->c / m->alpha0 - m->t0 ) {
      prod = 1;
      if ( err0 != NULL ) {
        *err0 = geterror(v, n);
        if ( m->verbose >= 1 ) {
          fprintf(stderr, "starting production at step %ld, t0 %g, err %g\n",
              t, m->t0, *err0);
        }
      }
    }

    /* compute the updating magnitude */
    if ( prod ) {
      a = m->c / (t - m->nequil + m->t0);
    } else {
      a = m->alpha0;
    }


    /* the distribution density is p = 1/n
     * this is why we have to multiply by n */
    a *= n;

    if ( m->nbn > 1 ) {
      mbin_update(v, n, i, a, m->nbs, m->nbn);
    } else {
      v[i] += a;
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

  free(v);
  free(vac);

  return err;
}



static double invt_test(invtpar_t *m)
{
  double err, see = 0;
  double err0, see0 = 0;
  int i;

  mtscramble( time(NULL) );
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
  return err;
}



int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  invt_test(m);

  return 0;
}
