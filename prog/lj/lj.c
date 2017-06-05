#include "lj.h"
#include "../ave.h"
#include <time.h>

int np = 100;
double rho = 0.8;
double rcdef = 100.0;
double tp = 3.0;
int mcblk = 5;
double delr = 0.01; /* spacing */
const char *fnvbias = "vbias.dat";
const char *fnvref = "vref.dat";
const char *fnlog = "verr.log";

static void ljpar_help(void)
{
  fprintf(stderr, "  --np=:      number of Lennard Jones particles, default %d\n", np);
  fprintf(stderr, "  --rho=:     density of Lennard Jones particles, default %g\n", rho);
  fprintf(stderr, "  --rc=:      distance cutoff, default %g\n", rcdef);
  fprintf(stderr, "  --tp=:      temperature, default %g\n", tp);
  fprintf(stderr, "  --dr=:      bin width, default %g\n", delr);
  fprintf(stderr, "  --blk=:     number of steps in a Monte Carlo, default %d\n", mcblk);
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
static int gammrun(invtpar_t *m, metad_t *metad, lj_t *lj)
{
  long t;
  int ir0, ir;
  double dr, sacc = 0;

  metad->a = m->alpha0;
  ir0 = dist01(metad, lj, &dr);

  fprintf(stderr, "starting constant magnitude %g metadynamics run of %ld steps...\n",
      metad->a, m->gam_nsteps);
  for ( t = 1; t <= m->gam_nsteps; t++ ) {
    sacc += lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    metad_updatev(metad, ir);
    metad->tmat[ir*metad->n + ir0] += 1;
    ir0 = ir;
    if ( t % 10 == 0 )
      metad_add_varv(metad);
    if ( t % 10000 == 0 ) fprintf(stderr, "t %ld/%ld = %5.2f%%, acc %.2f%% \r", t, m->gam_nsteps, 100.*t/m->gam_nsteps, 100*sacc/t);
    if ( t % 1000000 == 0 || t == m->gam_nsteps ) {
      metad_getgamma_varv(metad, m->alpha0, "gamma.dat");
      metad_getgamma_tmat(metad, 1, "tgamma.dat");
    }
  }
  metad_save(metad, fnvbias);
  return 0;
}


/* production run */
static double prodrun(invtpar_t *m, metad_t *metad, lj_t *lj,
    int prod, long nsteps, const char *fn, double *hfl, int *ntrip)
{
  int ir, n = metad->n, sgn = 0;
  long t;
  double t0, dr;

  metad->a = m->alpha0;
  t0 = 2/metad->a;
  metad_clearh(metad);
  *ntrip = 0;
  for ( t = 1; t <= nsteps; t++ ) {
    lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    if ( prod ) {
      if ( m->opta ) {
        //double a0 = intq_interpa(metad->intq, (double) t);
        metad->a = intq_evala(metad->intq, (double) t);
        //printf("t %ld, %g %g\n", t, a0, metad->a); if (t %10000 == 0) getchar();
      } else {
        metad->a = 1.0/(t + t0);
      }
    }
    /* update number of round trips */
    if ( ir == 0 ) {
      if ( sgn > 0 ) *ntrip += 1;
      sgn = -1;
    } else if ( ir == n - 1 ) {
      if ( sgn < 0 ) *ntrip += 1;
      sgn = 1;
    }
    if ( t % 10000 == 0 ) {
      *hfl = metad_hfl(metad, -1);
      fprintf(stderr, "%s t %8ld/%8ld = %5.2f%%, a %.3e, err %12.6e, hist. fl %12.6e/%5.3f%% ntrip %d %20s\r",
          (prod ? "prod." : "prep."), t, nsteps, 100.*t/nsteps,
          metad->a, metad_geterror(metad), (*hfl) * (*hfl), 100 * (*hfl), *ntrip, " ");
    }
    if ( t % 1000000 == 0 ) metad_save(metad, fn);
    metad_updatev(metad, ir);
  }
  metad_save(metad, fn);
  *hfl *= *hfl;
  return metad_geterror(metad);
}

static int work(invtpar_t *m)
{
  lj_t *lj;
  int ir, kcerr, ntrip;
  double dr, errtrunc, hfl;
  metad_t *metad;

  //mtscramble(clock());

  lj = lj_open(np, rho, rcdef);
  lj_energy(lj);

  metad = metad_openf(0, lj->l*0.5, delr, m->pbc,
      METAD_SHIFT_TAIL, m->gaussig, m->kc, m->win, m->winn);
  metad->a = m->alpha0;
  /* load the reference bias potential */
  metad_load(metad, metad->vref, fnvref);
  errtrunc = metad_errtrunc(metad, metad->vref, &kcerr, "vtrunc.dat");
  ir = dist01(metad, lj, &dr);
  fprintf(stderr, "n %d, rmax %g, r %d/%g, err trunc %g, mode %d\n", metad->n, metad->xmax, ir, dr, errtrunc, kcerr);

  /* constant magnitude run */
  if ( m->gammethod != GAMMETHOD_NONE
    && m->gammethod != GAMMETHOD_LOAD ) {
    prodrun(m, metad, lj, 0, m->nequil, "vg.dat", &hfl, &ntrip); /* equilibration */
    gammrun(m, metad, lj);
  }

  /* compute the optimal schedule and the error */
  metad_getalphaerr(metad, m->opta, (double) m->nsteps, m->gammethod,
      m->fngamma, m->sampmethod, m->alpha0, (double) m->nequil, &m->qT,
      m->qprec, m->alpha_nint, m->fnalpha);

  /* multiple production runs */
  {
    int i, itr, ntrip0, ntrip1;
    ave_t ei[1], ef[1], fi[1], ff[1];
    double erri, errf, hfli, hflf; // , *v0;
    FILE *fplog;

    ave_clear(ei);
    ave_clear(ef);
    ave_clear(fi);
    ave_clear(ff);
    fprintf(stderr, "starting %ld metadynamics run of %ld/%ld steps..., a %g, err %g -> %g\n",
        m->ntrials, m->nequil, m->nsteps, metad->a, metad->eiref, metad->efref);
    fplog = fopen(fnlog, "w");
    metad_saveheader(metad, fplog);
    fprintf(fplog, "# %ld %ld %ld %g %g %g\n",
        m->ntrials, m->nequil, m->nsteps, metad->a, metad->eiref, metad->efref);
    fclose(fplog);
    for ( itr = 0; itr < m->ntrials; itr++ ) {
      for ( i = 0; i < metad->n; i++ ) metad->v[i] = 0; // v0[i];
      erri = prodrun(m, metad, lj, 0, m->nequil, "vi.dat", &hfli, &ntrip0);
      ave_add(ei, erri);
      ave_add(fi, hfli);
      errf = prodrun(m, metad, lj, 1, m->nsteps, "vf.dat", &hflf, &ntrip1);
      ave_add(ef, errf);
      ave_add(ff, hflf);
      printf("%4d: %9.7f %9.7f(%5.2f) %9.7f %9.7f(%5.2f) | %9.7f %9.7f %9.7f %9.7f | %d %d\n",
          itr, errf, ef->ave, ef->ave/metad->efref, erri, ei->ave, ei->ave/metad->eiref, hflf, ff->ave, hfli, fi->ave, ntrip1, ntrip0);
      fplog = fopen(fnlog, "a");
      fprintf(fplog, "%d %.8f %.8f %.8f %.8f %d %d\n", itr, errf, erri, hflf, hfli, ntrip1, ntrip0);
      fclose(fplog);
    }
  }

  lj_close(lj);
  metad_close(metad);
  return 0;
}

int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  m->n = 0;
  m->alpha0 = 1e-4; /* updating magnitude for equilibration and gamma run */
  m->gam_nsteps = 10000000L;
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
