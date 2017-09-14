#define IS2_LB 5
#define L IS2_L
#include "is2.h"
#include "util.h"



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



/* modified Wolff algorithm */
__inline static int is2_wolff_mod(is2_t *is,
    double c1, double c2, double Eave)
{
  int l = is->l, n = is->n, i, ix, iy, id, s, cnt = 0, h = 0;
  double padd, enew, dv;

  padd = 1 - exp(-2*c1);

  /* randomly selected a seed */
  id = (int) ( rand01() * n );
  s = is->s[id];
  is->queue[ cnt++ ] = id;
  for ( i = 0; i < n; i++ )
    is->used[i] = 0;
  is->used[id] = (char) 1;

  /* go through spins in the queue */
  for ( i = 0; i < cnt; i++ ) {
    id = is->queue[i];
    is->s[id] = -s;
    /* add neighbors of i with the same spins */
    ix = id % l;
    iy = id - ix;
    h += is2_addtoqueue(is, iy + (ix + 1) % l,     s, padd, &cnt);
    h += is2_addtoqueue(is, iy + (ix + l - 1) % l, s, padd, &cnt);
    h += is2_addtoqueue(is, (iy + l) % n + ix,     s, padd, &cnt);
    h += is2_addtoqueue(is, (iy + n - l) % n + ix, s, padd, &cnt);
  }

  enew = is->E + 2 * s * h;
  dv = 2 * s * h * c2 * ((enew + is->E)/2 - Eave);
  //if ( dv <= 0 || rand01() < exp(-dv) ) {
  if ( metroacc(dv) ) {
    is->E += 2 * s * h;
    is->M -= 2 * s * cnt;
    return 1;
  } else { /* reject */
    for ( i = 0; i < cnt; i++ ) {
      id = is->queue[i];
      is->s[id] = s;
    }
    return 0;
  }
}


static void is2_aus(is2_t *is, double Eave, double Esig,
    int wolff, long nsteps)
{
  int id, h;
  int enew;
  long t, nacc = 0;
  double edev, c1 = 0, c2 = 0, dv, y1, y2;
  double *his;

  mtscramble(clock());
  xnew(his, is->n + 1);

  /* equilibrate the system to raise the energy */
  for ( t = 1; ; t++ ) {
    IS2_PICK(is, id, h);
    IS2_FLIP(is, id, h);
    if ( is->E > Eave ) break;
  }
  c1 = 0.44;

  for ( t = 1; t <= nsteps; t++ ) {
    if ( wolff ) {
      /* cluster algorithm */
      nacc += is2_wolff_mod(is, c1, c2, Eave);
    } else {
      /* Metropolis way */
      IS2_PICK(is, id, h);
      enew = is->E + 2 * h;
      /* compute the change of bias potential */
      edev = (enew + is->E)/2 - Eave;
      dv = 2 * h * (c1 + c2 * edev);
      //if ( dv <= 0 || rand01() < exp(-dv) ) {
      if ( metroacc(dv) ) {
        IS2_FLIP(is, id, h);
        nacc++;
      }
    }
    y1 = (is->E - Eave) / Esig;
    c1 += y1 / Esig / t;
    y2 = (y1 * y1 - 1);
    c2 += y2 / (Esig * Esig) / t;
    if ( t % 10000 == 0 )
      printf("t %ld, c1 %g, c2 %g, E %d, Eave %g, y %g, %g, acc %g%%\n",
          t, c1, c2, is->E, Eave, y1, y2, 100.*nacc/t);
    his[(is->E + 2 * is->n)/4] += 1;
  }
  savehist(his, is->n + 1, -2*IS2_N, 4, "is2aus.his");
  free(his);
}



int main(int argc, char **argv)
{
  is2_t *is;
  int method = 0;
  double Eave = -1468, Edev = 50;
  long nsteps = 0;

  if ( argc > 1 ) method = atoi( argv[1] );
  if ( argc > 2 ) Eave   = atof( argv[2] );
  if ( argc > 3 ) Edev   = atof( argv[3] );
  if ( argc > 4 ) nsteps = atol( argv[4] );
  if ( nsteps <= 0 )
    nsteps = (method == 0) ? 100000000L : 500000L;

  is = is2_open(L);
  is2_aus(is, Eave, Edev, method, nsteps);
  is2_close(is);
  return 0;
}

