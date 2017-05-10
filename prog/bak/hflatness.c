/* compare two definitions of histogram flatness */
#include "../util.h"
#include "../cosmodes.h"
#include "../mtrand.h"

int main(int argc, char **argv)
{
  int n = 100, pbc = 0;
  int i, j, k, km;
  double *u, *v, *hist, *costab;
  double um, vmax, vmin, fl, vsq;
  double rcnt, sum, sfl, svsq;
  long t, nsteps = 100000000L, period;
  FILE *fp;

  mtscramble(clock());

  if ( argc > 1 ) n = atoi(argv[1]);
  if ( argc > 2 ) nsteps = atol(argv[2]);
  if ( argc > 3 ) pbc = atoi(argv[3]);
  nsteps *= n;
  period = n * 10000;

  xnew(u, n);
  xnew(v, n);
  xnew(hist, n);
  costab = mkcostab(n, pbc);

  /* do sampling */
  i = 0;
  sum = sfl = 0;
  for ( t = 1; t <= nsteps; t++ ) {
    //i = (int) (rand01() * n);
    if ( pbc ) {
      i = (i + ((rand01() > 0.5) ? 1 : -1) + n) % n;
    } else {
      if ( rand01() > 0.5 ) { /* increase i */
        if ( i < n - 1 ) i++;
      } else { /* decrease i */
        if ( i > 0 ) i--;
      }
    }
    hist[i] += 1;

    if ( t % period == 0 ) {
      /* normalize the histogram and
       * compute the histogram flatness */
      vmax = 0;
      vmin = n;
      um = 0;
      vsq = 0;
      //FILE *fq = fopen("ah.dat", "w");
      for ( j = 0; j < n; j++ ) {
        v[j] = n*hist[j]/period;
        if ( v[j] > vmax ) vmax = v[j];
        if ( v[j] < vmin ) vmin = v[j];
        um += costab[1*n+j] * v[j];
        vsq += v[j] * v[j];
        hist[j] = 0;
        //fprintf(fq, "%d %g\n", j, v[j]-1);
      }
      //fclose(fq); getchar();
      fl = (vmax - vmin)/(vmax + vmin);
      um /= n;
      vsq /= n;
#if 0
      /* convert them to modes */
      getcosmodes(v, n, u, costab);
      /* compute the maximal modes */
      um = 0;
      for ( k = 1; k < n; k++ ) {
        if ( fabs(u[k]) > um )
          um = fabs(u[km = k]);
        //printf("%d %g\n", k, u[k]);
      }
      //printf("km %d, um %g\n", km, um);getchar();
#endif
      rcnt += 1;
      sum += um * um * period;
      sfl += fl * fl * period;
      svsq += vsq * vsq * period;
      if ( t % 10000000L == 0 )
        printf("%ld, fl %g(%g), um %g(%g), vsq %g(%g), ratio %g\n",
            t, fl, sqrt(sfl/rcnt), um, sqrt(sum/rcnt),
            vsq, sqrt(svsq/rcnt), sqrt(sum/sfl));
    }
  }

  free(u);
  free(v);
  free(hist);
  free(costab);

  return 0;
}
