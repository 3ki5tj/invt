/* Ornstein-Uhlenbeck process */



#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mtrand.h"



double gamdt = 0.01;
long nsteps = 100000000L;
int block = 100;



/* simulate the Ornstein-Uhlenbeck process */
static void ou(double dt)
{
  double x, sxx = 0;
  long t, nblk = 0;

  double expndt = exp(-0.5 * dt);
  double fluc = sqrt(2 * dt);

  x = randgaus();
  for ( t = 1; t <= nsteps; t++ ) {

    x *= expndt;
    x += fluc * randgaus();
    x *= expndt;

    sxx += x * x;

    if ( t % block == 0 ) {
      nblk++;
    }
  }
  printf("sxx %g\n", sxx/nsteps);
}



int main(void)
{
  mtscramble( clock() );
  ou( gamdt );
  return 0;
}
