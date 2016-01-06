#include "mtrand.h"
#include "util.h"
#include "invtpar.h"



/* return the root squared error of the inverse time scheme
 * */
static double comperr(const invtpar_t *p)
{
  double *v, a, dv, err, t0;
  int i, j, n = p->n;
  long t;

  xnew(v, n);
  for ( i = 0; i < n; i++ ) {
    v[i] = 0;
  }

  /* constant */
  t0 = p->c / p->alpha0;

  for ( t = 0; t < p->nsteps + p->nequil; t++ ) {
    /* MCMC sampling */
    if ( i < 0 || p->tcorr <= 0 || rand01() * p->tcorr < 1 ) {
      j = (int) ( n * rand01() );
      dv = v[j] - v[i];
      if ( dv < 0 || rand01() < exp(-dv) ) {
        i = j;
      }
    }

    /* compute the updating magnitude
     * the distribution density is p = 1/n
     * this is why we have to multiply by n */
    if ( t >= p->nequil ) {
      /* the constant t0 makes the transition
       * of alpha at t = t0 smooth */
      a = p->c * n / (t - p->nequil + t0);
    } else {
      a = p->alpha0 * n;
    }

    if ( p->nbn > 1 ) {
      double rem = 1, z;

      if ( i > 0 ) {
        z = p->nbs[1];
        v[i - 1] += a * z;
        rem -= z;
      }
      if ( i < n - 1 ) {
        z = p->nbs[1];
        v[i + 1] += a * z;
        rem -= z;
      }
      v[i] += a * rem;

    } else {
      v[i] += a;
    }
  }

  /* normalize */
  for ( a = 0, i = 0; i < n; i++ ) {
    a += v[i];
  }
  a /= n;
  for ( i = 0; i < n; i++ ) {
    v[i] -= a;
  }

  /* compute the error */
  err = 0;
  for ( i = 0; i < n; i++ ) {
    err += v[i] * v[i];
  }

  /* print out the error */
  if ( p->verbose ) {
    for ( i = 0; i < n; i++ ) {
      printf("  %+11.8f", v[i]);
    }
    printf("\n");
  }

  free(v);

  return sqrt(err / n);
}



static double invt_test(invtpar_t *p)
{
  double err, se = 0;
  int i;

  mtscramble( time(NULL) );
  for ( i = 0; i < p->ntrials; i++ ) {
    err = comperr(p);
    se += err;
    printf("%4d: err %10.8f, ave %10.8f\n", i, err, se/(i+1));
  }

  err = se / p->ntrials;
  printf("average error: %g\n", err);
  return err;
}



int main(int argc, char **argv)
{
  invtpar_t p[1];

  invtpar_init(p);
  invtpar_doargs(p, argc, argv);
  invtpar_dump(p);
  invt_test(p);

  return 0;
}
