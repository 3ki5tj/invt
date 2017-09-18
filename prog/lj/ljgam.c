#include "lj.h"
#include "../ave.h"
#include <time.h>

int np = 100;
double rho = 0.1;
double rcdef = 100.0;
double tp = 3.0;
int mcblk = 5;
int hblksz = 0;
double delr = 0.01; /* spacing */
long nstsave = 1000000; /* interval of saving data */
const char *fnvbias = "vbias.dat";
const char *fnvref = "vref.dat";
const char *fnlog = "verr.log";
const char *fnxerr = "xerr.dat";

static void ljpar_help(void)
{
  fprintf(stderr, "  --np=:      number of Lennard Jones particles, default %d\n", np);
  fprintf(stderr, "  --rho=:     density of Lennard Jones particles, default %g\n", rho);
  fprintf(stderr, "  --rc=:      distance cutoff, default %g\n", rcdef);
  fprintf(stderr, "  --tp=:      temperature, default %g\n", tp);
  fprintf(stderr, "  --dr=:      bin width, default %g\n", delr);
  fprintf(stderr, "  --blk=:     number of steps in a Monte Carlo, default %d\n", mcblk);
  fprintf(stderr, "  --hblk=:    number of steps in a histogram block, default %d\n", hblksz);
  fprintf(stderr, "  --vref=:    reference bias potential, default %s\n", fnvref);
}

static int ljpar_keymatch(invtpar_t *m, const char *key, const char *val)
{
  if ( strcmpfuzzy(key, "np") == 0 ) {
    np = invtpar_getint(m, key, val);
  } else if ( strcmpfuzzy(key, "rho") == 0 ) {
    rho = invtpar_getdouble(m, key, val);
  } else if ( strcmpfuzzy(key, "rc") == 0
           || strcmpfuzzy(key, "rcutoff") == 0 ) {
    rcdef = invtpar_getdouble(m, key, val);
  } else if ( strcmpfuzzy(key, "tp") == 0
           || strcmpfuzzy(key, "temp") == 0 ) {
    tp = invtpar_getdouble(m, key, val);
  } else if ( strcmpfuzzy(key, "dr") == 0
           || strcmpfuzzy(key, "delr") == 0 ) {
    delr = invtpar_getdouble(m, key, val);
  } else if ( strcmpfuzzy(key, "blk") == 0
           || strcmpfuzzy(key, "mcblk") == 0 ) {
    mcblk = invtpar_getint(m, key, val);
  } else if ( strcmpfuzzy(key, "vref") == 0
           || strcmpfuzzy(key, "fnvref") == 0 ) {
    fnvref = val;
  } else if ( strcmpfuzzy(key, "hblk") == 0
           || strcmpfuzzy(key, "hblksz") == 0 ) {
    hblksz = invtpar_getint(m, key, val);
  }
  return 0;
}

/* return the index for distance between the first two atoms */
static int dist01(metad_t *metad, lj_t *lj, double *pdr)
{
  static double invl = 0, dx[D], dr;
  if ( invl <= 0 ) invl = 1/lj->l; /* assuming l doesn't change */
  dr = sqrt( lj_pbcdist2(dx, lj->x[0], lj->x[1], lj->l, invl) );
  if ( pdr != NULL ) *pdr = dr;
  return metad_getindexf(metad, dr);
}

/* a block of Metropolis move */
static double lj_metroblk(lj_t *lj, metad_t *metad)
{
  int it, i, sacc = 0;
  double amp = 0.1/rho;

  /* move for non-special */
  for ( it = 0; it < mcblk - 1; it++ ) {
    i = 2 + (int) (rand01() * (lj->n - 2));
    sacc += lj_metro(lj, i, amp, 1/tp, metad);
  }
  i = (int) (rand01() * 2);
  sacc += lj_metro(lj, i, amp, 1/tp, metad);
  return 1. * sacc / mcblk;
}

/* constant updating magnitude run
 * to compute the gamma values */
static int eqlrun(metad_t *metad, lj_t *lj, long nsteps)
{
  long t;
  int ir0, ir;
  double dr, sacc = 0, *hist;

  metad->a = 0;

  if ( hblksz <= 0 ) {
    hblksz = (int)(0.1/metad->a);
    if ( hblksz <= 0 ) hblksz = 0;
  }

  ir0 = dist01(metad, lj, &dr);
  xnew(hist, metad->n);
  for ( ir = 0; ir < metad->n; ir++ ) hist[ir] = 0;

  fprintf(stderr, "starting equilibirium run of %ld steps (updating magnitude %g)...%20s\n",
      nsteps, metad->a, "");
  metad_hblk_clear(metad);
  metad_clearav(metad); /* clear the average correction data */
  for ( t = 1; t <= nsteps; t++ ) {
    sacc += lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    hist[ir] += 1;
    metad->tmat[ir*metad->n + ir0] += 1;
    metad_hblk_add(metad, ir);
    ir0 = ir;
    if ( t % hblksz == 0 ) {
      metad_hblk_dump(metad);
    }
    if ( t % 10000 == 0 ) fprintf(stderr, "t %ld/%ld = %5.2f%%, acc %.2f%% \r", t, nsteps, 100.*t/nsteps, 100*sacc/t);
    if ( t % nstsave == 0 || t == nsteps ) {
      metad_getgamma_tmat(metad, 1, metad->tgamma, "tgamma.dat");
      metad_getgamma_hblk(metad, metad->hgamma, "hgamma.dat");
      {
        int i;
        FILE *fp = fopen("ljhist.dat", "w");
        for ( i = 0; i < metad->n; i++ )
          fprintf(fp, "%d %g\n", i, hist[i]);
        fclose(fp);
      }
      metad_save(metad, fnvbias);
    }
  }
  fprintf(stderr, "\n");
  return 0;
}


static int work(invtpar_t *m)
{
  lj_t *lj;
  int ir, n;
  double dr;
  metad_t *metad;

  mtscramble(clock());

  lj = lj_open(np, rho, rcdef);
  lj_energy(lj);

  metad = metad_openf(0, lj->l*0.5, delr, m->pbc,
      METAD_SHIFT_TAIL, m->gaussig, m->kc, m->win, m->winn);
  n = metad->n;
  /* load the reference bias potential */
  metad_load(metad, metad->vref, fnvref);
  metad_load(metad, metad->v, fnvref);
  ir = dist01(metad, lj, &dr);
  fprintf(stderr, "n %d, rmax %g, r %d/%g\n", n, metad->xmax, ir, dr);

  eqlrun(metad, lj, m->gam_nsteps);

  lj_close(lj);
  metad_close(metad);
  return 0;
}

int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  m->n = 0;
  m->alpha0 = 0;
  m->gam_nsteps = 10000000L;
  m->gam_nstave = 100;
  m->nequil = 10000000L;
  m->nsteps = 10000000L;
  m->ntrials = 1000;
  //m->alpha_nint = 20000;
  m->pbc = 0;
  m->gammethod = GAMMETHOD_NONE;
  m->sampmethod = SAMPMETHOD_GAUSS;
  strcpy(m->fnalpha, "alpha.dat");
  m->userhelp = ljpar_help;
  m->usermatch = ljpar_keymatch;
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  work(m);
  invtpar_finish(m);
  return 0;
}
