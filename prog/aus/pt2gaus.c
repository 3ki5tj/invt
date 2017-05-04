#define POTTS2_LB 5
#include "potts2.h"
#include "util.h"
#include "gaus.h"




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
    int wolff, long nsteps)
{
  long t;
  int m, id, acc;
  double ecmin, ecmax, esig;
  const double tp = 1.4;
  gaus_t *gaus;

  //mtscramble(clock());

  ecmin = -1.7 * pt->n;
  ecmax = -0.9 * pt->n;
  //ecmin = -1.3 * pt->n;
  //ecmax = -1.25 * pt->n;
  //ecmin = -0.5 * pt->n;
  //ecmax = -0.45 * pt->n;
  esig = pt->l;
  m = (int) ((ecmax - ecmin) / esig) + 1;
  potts2_equil(pt, ecmin);
  gaus = gaus_open(ecmin, ecmax, m, esig, tp,
      -2*pt->n, 0, 1);

  id = 0;
  for ( t = 1; t <= nsteps; t++ ) {
    if ( wolff ) { /* cluster algorithm */
      acc = potts2_wolff_mod(pt, gaus->beta1[id], gaus->beta2[id], gaus->ave[id]);
    } else { /* Metropolis algorithm */
      acc = potts2_metro_mod(pt, gaus->beta1[id], gaus->beta2[id], gaus->ave[id]);
    }
    gaus_update(gaus, id, pt->E, (double) t);
    gaus_add(gaus, id, pt->E, acc);
    gaus_bmove(gaus, pt->E, &id);
    if ( t % 100 == 0 ) {
      gaus_wlcheck(gaus, 0.1, 0.5, (double) t);
    }
    if ( t % 100000 == 0 ) {
      printf("t %ld, flatness %g, alpha %g, id %d, invt %d\n",
          t, gaus->hflatness, gaus->alpha, id, gaus->invt);
      gaus_savehist(gaus, "pt2gaus.his");
    }
  }
  gaus_close(gaus);
}



int main(int argc, char **argv)
{
  potts2_t *pt;
  int method = 0, q = 10;
  long nsteps = 0;

  if ( argc > 1 ) method = atoi( argv[1] );
  if ( argc > 2 ) q      = atoi( argv[2] );
  if ( argc > 3 ) nsteps = atol( argv[3] );
  if ( nsteps <= 0 )
    nsteps = (method == 0) ? 100000000L : 5000000L;

  pt = potts2_open(POTTS2_L, q);
  potts2_gaus(pt, method, nsteps);
  potts2_close(pt);
  return 0;
}

