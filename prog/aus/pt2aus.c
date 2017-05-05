#define POTTS2_LB 5
#include "potts2.h"
#include "util.h"
#include "mmwl.h"



static void savehist(double *h, int n,
    double emin, double de, const char *fn)
{
  FILE *fp;
  int i;
  double tot = 0;

  fp = fopen(fn, "w");
  for ( i = 0; i < n; i++ )
    tot += h[i];
  for ( i = 0; i < n; i++ ) {
    if ( h[i] <= 0 ) continue;
    fprintf(fp, "%g %g %g\n", emin+i*de, h[i]/tot/de, h[i]);
  }
  fclose(fp);
}



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


static void potts2_aus(potts2_t *pt, double Eave, double Esig,
    int wolff, long nsteps)
{
  int id, h, sn, dovar = 1;
  long t, nacc = 0, nstrep;
  double c1 = 0, c2 = 0, y1, y2, amp, beta1, beta2;
  double *his;
  mmwl_t mmwl[1];

  dovar = ( Esig > 0 );
  Esig = fabs(Esig);

  mtscramble(clock());
  xnew(his, pt->n * 2 + 1);

  /* equilibrate the system to raise the energy */
  for ( t = 1; ; t++ ) {
    POTTS2_PICK(pt, id, h, sn);
    POTTS2_FLIP(pt, id, h, sn);
    if ( pt->E > Eave ) break;
  }
  c1 = 1.4245 * Esig;
  mmwl_init(mmwl, 1e-3);

  nstrep = wolff ? 50000 : 10000000;
  for ( t = 1; t <= nsteps; t++ ) {
    beta1 = c1 / Esig;
    beta2 = c2 / (Esig * Esig);
    if ( wolff ) { /* cluster algorithm */
      nacc += potts2_wolff_mod(pt, beta1, beta2, Eave);
    } else { /* Metropolis algorithm */
      nacc += potts2_metro_mod(pt, beta1, beta2, Eave);
    }
    amp = mmwl_getalpha(mmwl);
    y1 = (pt->E - Eave) / Esig;
    c1 += y1 * amp;
    y2 = y1 * y1 - 1;
    if ( dovar ) c2 += y2 * amp;
    his[pt->E + 2 * pt->n] += 1;
    mmwl_add(mmwl, y1, y2);
    /* control the updating magnitude */
    if ( t % 100 == 0 && mmwl_check(mmwl, dovar, 0.05, 0.5) ) {
      printf("t %ld, new updating magnitude %g, fl %g, %g, c %g, %g, invt %d\n",
          t, mmwl_getalpha(mmwl), mmwl->fl[1], mmwl->fl[2], c1, c2, mmwl->invt);
    }
    if ( t % nstrep == 0 ) {
      printf("t %ld, c1 %g, c2 %g, E %d, Eave %g, y %g(%+g), %g(%+g), amp %g, acc %g%%, invt %d\n",
          t, c1, c2, pt->E, Eave, y1, mmwl->fl[1], y2, mmwl->fl[2], amp, 100.*nacc/t, mmwl->invt);
      savehist(his, pt->n * 2 + 1, -2*POTTS2_N, 1, "pt2aus.his");
      //getchar();
    }
  }
  free(his);
}



int main(int argc, char **argv)
{
  potts2_t *pt;
  int method = 0, q = 10;
  double eave = -1.3, edev = 1.0;
  long nsteps = 0;

  if ( argc > 1 ) method = atoi( argv[1] );
  if ( argc > 2 ) eave   = atof( argv[2] );
  if ( argc > 3 ) edev   = atof( argv[3] );
  if ( argc > 4 ) nsteps = atol( argv[4] );
  if ( nsteps <= 0 )
    nsteps = (method == 0) ? 1000000000L : 500000000L;

  pt = potts2_open(POTTS2_L, q);
  potts2_aus(pt, pt->n * eave, pt->l * edev, method, nsteps);
  potts2_close(pt);
  return 0;
}

