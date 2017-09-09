#ifndef L
#ifndef POTTS2_LB
#define POTTS2_LB 5
#endif
#endif
#include "potts2.h"
#include "util.h"
#include "gaus.h"



const char *fnhis = "pt2gaus.his";
const char *fndat = "pt2gaus.dat";



enum { SAMP_METROPOLIS, SAMP_WOLFF };
const char *samp_names[] = {"Metropolis", "Wolff"};
double localmovef = 1.0;


/* modified Metropolis algorithm
 * p(E) ~ exp(-b1*(E-Eave)-b2*(E-Eave)^2/2) */
__inline static int potts2_metro_mod(potts2_t *pt,
    double b1, double b2, double Eave)
{
  int id, sn, h, enew;
  double edev, dv;

#ifndef L
  POTTS2_PICK(pt, id, sn, h);
#else
  id = potts2_pick(pt, &sn, &h);
#endif
  enew = pt->E + h;
  /* compute the change of the bias potential */
  edev = (enew + pt->E)*.5 - Eave;
  dv = h * (b1 + b2 * edev);
  if ( dv <= 0 || rand01() < exp(-dv) ) {
#ifndef L
    POTTS2_FLIP(pt, id, sn, h);
#else
    potts2_flip(pt, id, sn, h);
#endif
    return 1;
  }
  return 0;
}



/* modified Wolff algorithm
 * p(E) ~ exp(-b1*(E-Eave)-b2*(E-Eave)^2/2) */
__inline static int potts2_wolff_mod(potts2_t *pt,
    double b1, double b2, double Eave)
{
  int l = pt->l, q = pt->q, n = pt->n, i, ix, iy, id, so, sn, cnt = 0, h = 0, enew;
  double padd, dv;

  padd = 1 - exp(-b1);

  /* randomly selected a seed */
  id = (int) ( rand01() * n );
  so = pt->s[id];
  sn = (so + 1 + (int) (rand01() * (q - 1))) % q;
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
  //printf("spin %d, cnt %d, E %d %d %d\n", id0, cnt, pt->E, h, enew);
  dv = h * b2 * ((enew + pt->E)*.5 - Eave);
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
  int id, sn, h;
  long t;

  for ( t = 1; ; t++ ) {
    id = potts2_pick(pt, &sn, &h);
    potts2_flip(pt, id, sn, h);
    if ( pt->E > ene ) break;
  }
}



static void potts2_gaus(potts2_t *pt,
    int sampmethod, int lnzmethod, long nsteps, double sig)
{
  long t, nstsave, ntrips = 0;
  int id, acc, sgn = 0;
  double ecmin, ecmax, esig, espacing, beta1, beta2, fl, alpha0;
  const double beta_c = 1.4;
  const double magred = 0.5;
  gaus_t *gaus;

  //mtscramble(clock());

  ecmin = -(int) (1.8 * pt->n + 0.5);
  ecmax = -(int) (0.8 * pt->n + 0.5);
  esig = sig * pt->l;
  espacing = esig;
  potts2_equil(pt, ecmin);
  if (lnzmethod == LNZ_WL) {
    fl = 0.2;
    alpha0 = 0.01;
  } else {
    fl = 0.05;
    alpha0 = 0.001;
  }
  gaus = gaus_open(ecmin, ecmax, espacing, esig, lnzmethod,
      beta_c * esig, alpha0, -2*pt->n, 0, 1);

  fprintf(stderr, "samping-method %s, lnz-method %s, sig %g, spacing %g, localmove %g%%\n",
      samp_names[sampmethod], lnz_names[lnzmethod], esig, espacing, localmovef*100.0);

  id = 0;
  sgn = 0;
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
    gaus_move(gaus, pt->E, &id, rand01() <= localmovef);
    gaus_add(gaus, id, pt->E, acc);

    /* update the number of round trips */
    if ( (sgn == 0 && id == gaus->n - 1)
      || (sgn == 1 && id == 0) ) {
      sgn = !sgn;
      ntrips++;
    }

    if ( t % 10000 == 0 ) {
      gaus_wlcheckx(gaus, fl, magred);
      //printf("t %d, fl %g, id %d, E %d, c1 %g, c2 %g\n", t, gaus->hfluc, id, pt->E, gaus->c1[id]/esig, gaus->c2[id]); if ( t % 100 == 0 ) getchar();
    }
    if ( t % nstsave == 0 ) {
      double alpha, alphamm;
      alpha = gaus_getalpha(gaus, id, &alphamm);
      gaus_save(gaus, fndat);
      gaus_savehist(gaus, fnhis);
      printf("t %ld/%g, fluc %g(%.0f%%/%.0f%%/%.0f%%), alpha %g/%g, id %d, invt %d, trips %ld\n",
          t, gaus->t, gaus->hfluc, 100*gaus->flfr[0], 100*gaus->flfr[1], 100*gaus->flfr[2],
          alpha, alphamm, id, gaus->invt, ntrips);
    }
  }
  gaus_close(gaus);
}



int main(int argc, char **argv)
{
  potts2_t *pt;
  int sampmethod = SAMP_METROPOLIS, lnzmethod = LNZ_WL, q = 10;
  long nsteps = 0;
  double sig = 1.0;
#ifndef L
  int l = POTTS2_L;
#else
  int l = L;
#endif

  if ( argc > 1 ) sampmethod = atoi( argv[1] );
  if ( argc > 2 ) lnzmethod  = atoi( argv[2] );
  if ( argc > 3 ) nsteps     = atol( argv[3] );
  if ( argc > 4 ) sig        = atof( argv[4] );
  if ( argc > 5 ) localmovef = atof( argv[5] );
  if ( argc > 6 ) q          = atoi( argv[6] );
  if ( nsteps <= 0 )
    nsteps = (sampmethod == 0) ? 100000000L : 100000L;

  pt = potts2_open(l, q);
  nsteps *= pt->n;
  potts2_gaus(pt, sampmethod, lnzmethod, nsteps, sig);
  potts2_close(pt);
  return 0;
}

