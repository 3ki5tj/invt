

#include "invt.h"



/* return the root squared error of the inverse time scheme
 * */
static double comperr(const invtpar_t *m)
{
  double *v, *vac, a, err;
  int i, n = m->n;
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

    /* compute the updating magnitude */
    if ( t >= m->nequil + m->c / m->alpha0 - m->t0 ) {
      a = m->c / (t - m->nequil + m->t0);
    } else {
      a = m->alpha0;
    }

    /* the distribution density is p = 1/n
     * this is why we have to multiply by n */
    a *= n;

    if ( m->nbn > 1 ) {
      double rem = 1, z;

      if ( i > 0 ) {
        z = m->nbs[1];
        v[i - 1] += a * z;
        rem -= z;
      }
      if ( i < n - 1 ) {
        z = m->nbs[1];
        v[i + 1] += a * z;
        rem -= z;
      }
      v[i] += a * rem;

    } else {
      v[i] += a;
    }
  }

  /* normalize */
  normalize(v, n);

  /* compute the error */
  err = geterror(v, n);

  /* print out the error */
  if ( m->verbose ) {
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
  int i;

  mtscramble( time(NULL) );
  for ( i = 0; i < m->ntrials; i++ ) {
    err = comperr(m);
    see += err * err;
    printf("%4d: err %10.8f, ave %10.8f, sqr %e\n",
        i, err, sqrt(see/(i+1)), see/(i+1));
  }

  err = see / m->ntrials;
  printf("average error: %10.8f, sqr %e\n", sqrt(err), err);
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
