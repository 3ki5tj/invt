/* test program for simple Langevin equation */



#include <stdio.h>
#include <math.h>
#include "mtrand.h"
#include "corr.h"



int gausdev = 0;
double alpha = 0.0001;
long nsteps = 100000000L;
int block = 10000;



/* simulate a simple Langevin equation */
static void langevin(void)
{
  double x = alpha, dx, u;
  double suu = 0;
  long t, nblk = 0;
  corr_t *corr = corr_open(1, nsteps / block + 1);

  for ( t = 1; t <= nsteps; t++ ) {
    dx = x;
    if ( gausdev ) {
      dx += randgaus();
    } else {
      dx += (rand01() > 0.5 ? 1 : -1);
    }
    x -= alpha * dx;
    if ( t % block == 0 ) {
      u = x;
      corr_add(corr, &u);
      suu += u * u;
      nblk++;
    }
  }
  printf("suu %g\n", suu/nblk);

  corr_save(corr, block, block * 100, 1e-4, 0, "langcorr.dat");
  corr_close(corr);
}



int main(void)
{
  langevin();
  return 0;
}
