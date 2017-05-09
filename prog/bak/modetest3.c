/* test the residue error due to mode truncation */
#include "../util.h"
#include "../cosmodes.h"
#include "../mtrand.h"

int main(int argc, char **argv)
{
  int n = 100, pbc = 0;
  int i, j, k, km;
  double *u, *v, *hist, *costab;
  double um, vmax, vmin, fl, rcnt = 0, srat = 0;
  long t, nsteps = 100000000L;
  FILE *fp;

  mtscramble(clock());

  if ( argc > 1 ) n = atoi(argv[1]);
  if ( argc > 2 ) nsteps = atol(argv[2]);
  if ( argc > 3 ) pbc = atoi(argv[3]);
  nsteps *= n;

  xnew(u, n);
  xnew(v, n);
  xnew(hist, n);
  costab = mkcostab(n, pbc);

  /* do sampling */
  i = 0;
  for ( t = 1; t <= nsteps; t++ ) {
    //i = (int) (rand01() * n);
    if ( pbc ) {
      i = (i + ((rand01() > 0.5) ? 1 : -1) + n) % n;
    } else {
      if ( rand01() > 0.5 ) {
        if ( i < n - 1) i++;
      } else {
        if ( i > 0) i--;
      }
    }
    hist[i] += 1;

    if ( t % n == 0 ) {
      for ( j = 0; j < n; j++ ) {
        v[j] = n*hist[j]/t;
      }
      /* convert them to modes */
      getcosmodes(v, n, u, costab);
      /* compute the maximal modes */
      um = 0;
      for ( k = 1; k < n; k++ ) {
        if ( fabs(u[k]) > um )
          um = fabs(u[km = k]);
        //printf("%d %g\n", k, u[k]);
      }
      //getchar();
      /* compute the histogram flatness */
      vmax = 0;
      vmin = n;
      for ( j = 0; j < n; j++ ) {
        if ( v[j] > vmax ) vmax = v[j];
        if ( v[j] < vmin ) vmin = v[j];
      }
      fl = 2*(vmax - vmin)/(vmax + vmin);
      rcnt += 1;
      srat += um*um/(fl*fl);
      if ( t % 10000000L == 0 )
        printf("%ld, fl %g(%g), um %g(%g), ratio %g\n",
            t, fl, fl*fl*t, um, um*um*t, sqrt(srat/rcnt));
    }
  }

  free(u);
  free(v);
  free(hist);
  free(costab);

  return 0;
}
