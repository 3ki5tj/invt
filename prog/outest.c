/* Ornstein-Uhlenbeck process */



#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mtrand.h"
//#include "corr.h"



double gamdt = 0.01;
long nsteps = 100000000L;
int block = 100;



/* simulate the Ornstein-Uhlenbeck process */
static void ou(double dt)
{
  double x, sxx = 0;
  long t, nblk = 0;
  //corr_t *corr = corr_open(1, nsteps / block + 1);

  double expndt = exp(-0.5 * dt);
  double fluc = sqrt(2 * dt);

  x = randgaus();
  for ( t = 1; t <= nsteps; t++ ) {

    x *= expndt;
    x += fluc * randgaus();
    x *= expndt;

    sxx += x * x;

    if ( t % block == 0 ) {
      //u = x;
      //corr_add(corr, &u);
      nblk++;
    }
  }
  printf("sxx %g\n", sxx/nsteps);

  //corr_save(corr, block, block * 100, 1e-4, 0, "langcorr.dat");
  //corr_close(corr);
}



int main(void)
{
  mtscramble( clock() );
  ou( gamdt );
  return 0;
}
