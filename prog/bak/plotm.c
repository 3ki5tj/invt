/* plot the mass function against qbar
 *
 *   icc plotm.c && ./a.out --fngamma=../lj/dr0.01/rho0.1/gamma.dat --sig=20
 *
 * Note: `sig` should be in the number of bins.
 * */

#include "../invt.h"
#include "../intq.h"

static int loadn(const char *fn)
{
  FILE *fp;
  int n;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return 0;
  }
  fscanf(fp, "# %d", &n);
  fclose(fp);
  return n;
}

static int plotm(int n, const double *lambda, const double *gamma,
    const char *fnmass)
{
  int i, k;
  double q, m;
  FILE *fp;

  if ( (fp = fopen(fnmass, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fnmass);
    return -1;
  }
  for ( i = 0; i < 100000; i++ ) {
    q = i * 0.1;
    for ( m = 0, k = 1; k < n; k++ ) {
      m += gamma[k] * lambda[k] * lambda[k] * exp(-2*lambda[k]*q);
    }
    fprintf(fp, "%g %g\n", q, sqrt(m));
  }
  fclose(fp);
  return -1;
}


int main(int argc, char **argv)
{
  invtpar_t m[1];
  double *gamma, *lambda, *win;
  int n, winn;

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  n = loadn(m->fngamma);
  xnew(gamma, n);
  loadgamma(n, gamma, m->fngamma);
  xnew(lambda, n);
  xnew(win, n);
  prepwin(lambda, n, m->win, m->winn, win, &winn,
      m->pbc, m->gaussig, m->kc,
      m->fnwin, m->fnwinmat, m->verbose);

  plotm(n, lambda, gamma, "mass.dat");
  //savegamma(n, lambda, "lam.dat");

  invtpar_finish(m);
  free(gamma);
  free(lambda);
  free(win);
  return 0;
}
