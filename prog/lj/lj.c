#include "lj.h"
#include <time.h>

int n = 55;
double rho = 0.8;
double rcdef = 100.0;
long nsteps = 10000000;
double tp = 1.5;

static void savehist(const char *fn,
    const double *vb, const double *his,
    int xn, double dx)
{
  FILE *fp;
  int i;
  double y, s = 0;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(fp, "cannot open %s\n", fn);
    return;
  }
  for ( i = 0; i < xn; i++ ) s += his[i];
  s = 1./(dx*s);

  for ( i = 0; i < xn; i++ ) {
    y = his[i];
    fprintf(fp, "%g %g %g\n", (i+.5)*dx, y*s, vb[i]);
  }
  fclose(fp);
}

static int foo(void)
{
  lj_t *lj;
  long t;
  int ir, acc;
  double dx[D], dr, l, invl, amp;
  double sacc = 0, lnf = 0.01;
  double *vb, xn, xdel = 0.02, *his;

  mtscramble(clock());

  lj = lj_open(n, rho, rcdef);
  xn = (int) (lj->l*0.5 / xdel);
  xnew(vb, xn);
  xnew(his, xn);
  lj_energy(lj);
  invl = 1/(l = lj->l);
  amp = 0.1/rho;
  dr = sqrt( lj_pbcdist2(dx, lj->x[0], lj->x[1], l, invl) );
  printf("%g %g\n", dr, xn * xdel);
  for ( t = 1; t <= nsteps; t++ ) {
    acc = lj_metro(lj, amp, 1/tp, vb, xn, xdel);
    sacc += acc;
    if ( t % 10 == 0 ) {
      dr = sqrt( lj_pbcdist2(dx, lj->x[0], lj->x[1], l, invl) );
      ir = (int) (dr / xdel);
      if ( ir < xn ) {
        vb[ir] += lnf;
        his[ir] += 1;
      }
    }
  }
  savehist("x.his", vb, his, xn, xdel);
  printf("tp %g, acc %g%%\n", tp, 100.*sacc/nsteps);
  lj_close(lj);
  free(vb);
  free(his);
  return 0;
}

int main(void)
{
  foo();
  return 0;
}
