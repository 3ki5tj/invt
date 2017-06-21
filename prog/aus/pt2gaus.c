#ifndef POTTS2_LB
#define POTTS2_LB 5
#endif
#include "potts2.h"
#include "util.h"
#include "gaus.h"



const char *fnhis = "pt2gaus.his";
const char *fndat = "pt2gaus.dat";



enum { SAMP_METROPOLIS, SAMP_WOLFF };



/* modified Metropolis algorithm */
__inline static int potts2_metro_mod(potts2_t *pt,
    double c1, double c2, double Eave)
{
  int id, h, sn, enew;
  double edev, dv;

  POTTS2_PICK(pt, id, h, sn);
  enew = pt->E + h;
  /* compute the change of bias potential */
  edev = (enew + pt->E)*.5 - Eave;
  dv = h * (c1 + c2 * edev);
  if ( dv <= 0 || rand01() < exp(-dv) ) {
    POTTS2_FLIP(pt, id, h, sn);
    return 1;
  }
  return 0;
}



/* modified Wolff algorithm */
__inline static int potts2_wolff_mod(potts2_t *pt,
    double c1, double c2, double Eave)
{
  int l = pt->l, n = pt->n, i, ix, iy, id, so, sn, cnt = 0, h = 0;
  double padd, enew, dv;

  padd = 1 - exp(-c1);

  /* randomly selected a seed */
  id = (int) ( rand01() * n );
  so = pt->s[id];
  sn = (so + 1 + (int) (rand01() * pt->q)) % pt->q;
  pt->queue[ cnt++ ] = id;
  for ( i = 0; i < n; i++ )
    pt->used[i] = 0;
  pt->used[id] = (char) 1;

  /* go through spins in the queue */
  for ( i = 0; i < cnt; i++ ) {
    id = pt->queue[i];
    pt->s[id] = sn;
    /* add neighbors of i with the same spins */
    ix = id % l;
    iy = id - ix;
    h += potts2_addtoqueue(pt, iy + (ix + 1) % l,     so, sn, padd, &cnt);
    h += potts2_addtoqueue(pt, iy + (ix + l - 1) % l, so, sn, padd, &cnt);
    h += potts2_addtoqueue(pt, (iy + l) % n + ix,     so, sn, padd, &cnt);
    h += potts2_addtoqueue(pt, (iy + n - l) % n + ix, so, sn, padd, &cnt);
  }

  enew = pt->E + h;
  dv = h * c2 * ((enew + pt->E)*.5 - Eave);
  if ( dv <= 0 || rand01() < exp(-dv) ) {
    pt->E += h;
    return 1;
  } else { /* reject */
    for ( i = 0; i < cnt; i++ ) {
      id = pt->queue[i];
      pt->s[id] = so;
    }
    return 0;
  }
}


/* equilibrate the system to raise the energy above ene */
static void potts2_equil(potts2_t *pt, double ene)
{
  int id, h, sn;
  long t;

  for ( t = 1; ; t++ ) {
    POTTS2_PICK(pt, id, h, sn);
    POTTS2_FLIP(pt, id, h, sn);
    if ( pt->E > ene ) break;
  }
}



static void potts2_gaus(potts2_t *pt,
    int sampmethod, int lnzmethod, long nsteps)
{
  long t, nstsave;
  int m, id, acc;
  double ecmin, ecmax, esig, beta1, beta2, fl, alpha0;
  const double beta_c = 1.4;
  gaus_t *gaus;

  //mtscramble(clock());

  ecmin = -1.7 * pt->n;
  ecmax = -0.9 * pt->n;
  esig = pt->l;
  m = (int) ((ecmax - ecmin) / esig) + 1;
  //esig *= 0.3; // make the Gaussian width narrower to test the correction
  potts2_equil(pt, ecmin);
  if (lnzmethod == LNZ_WL) {
    fl = 0.2;
    alpha0 = 0.01;
  } else {
    fl = 0.05;
    alpha0 = 0.001;
  }
  gaus = gaus_open(ecmin, ecmax, m, esig, lnzmethod,
      beta_c * esig, alpha0, -2*pt->n, 0, 1, 0);

  id = 0;
  if (sampmethod == SAMP_WOLFF) {
    nstsave = 100000;
  } else {
    nstsave = 10000000;
  }

  for ( t = 1; t <= nsteps; t++ ) {
    beta1 = gaus->c1[id] / esig;
    beta2 = gaus->c2[id] / (esig * esig);
    if ( sampmethod == SAMP_WOLFF ) { /* cluster algorithm */
      acc = potts2_wolff_mod(pt, beta1, beta2, gaus->ave[id]);
    } else { /* Metropolis algorithm */
      acc = potts2_metro_mod(pt, beta1, beta2, gaus->ave[id]);
    }
    gaus_add(gaus, id, pt->E, acc);
    gaus_move(gaus, pt->E, &id);
    if ( t % 100 == 0 ) {
      gaus_wlcheckx(gaus, fl, 0.5);
    }
    if ( t % nstsave == 0 ) {
      double alpha, alphamm;
      alpha = gaus_getalpha(gaus, id, &alphamm);
      gaus_save(gaus, fndat);
      gaus_savehist(gaus, fnhis);
      printf("t %ld/%g, flatness %g, alpha %g/%g, id %d, invt %d\n",
          t, gaus->t, gaus->hflatness, alpha, alphamm, id, gaus->invt);
    }
  }
  gaus_close(gaus);
}



int main(int argc, char **argv)
{
  potts2_t *pt;
  int sampmethod = SAMP_METROPOLIS, lnzmethod = LNZ_WL, q = 10;
  long nsteps = 0;

  if ( argc > 1 ) sampmethod = atoi( argv[1] );
  if ( argc > 2 ) lnzmethod  = atoi( argv[2] );
  if ( argc > 3 ) q          = atoi( argv[3] );
  if ( argc > 4 ) nsteps     = atol( argv[4] );
  if ( nsteps <= 0 )
    nsteps = (sampmethod == 0) ? 100000000L : 100000L;

  pt = potts2_open(POTTS2_L, q);
  nsteps *= pt->n;
  potts2_gaus(pt, sampmethod, lnzmethod, nsteps);
  potts2_close(pt);
  return 0;
}

