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

#if 0
/* run with decreasing magnitude
 * until it falls under m->alpha0
 * cannot be used for Gaussian updating scheme */
__inline static void decmagrun(invtpar_t *m, metad_t *metad, lj_t *lj)
{
  int ir;
  long t;
  double dr;

  ir = dist01(metad, lj, &dr);
  for ( t = 1; ; t++ ) {
    lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    metad_updatev(metad, ir);
    if ( t % 1000 == 0 ) {
      int sacc = metad_wlcheck(metad, m->fluc, m->magred);
      if ( sacc && metad->a < m->alpha0 )
        break;
    }
    if ( t % 1000000 == 0 ) {
      printf("t %ld, fl %g\n", t, metad->hfl);
      metad_save(metad, fnvbias);
    }
  }
  metad_save(metad, fnvbias);
}
#endif

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
    int prod, long nsteps, const char *fn, double *hfl)
{
  int ir;
  long t;
  double t0, dr;

  metad->a = m->alpha0;
  t0 = 2/metad->a;
  metad_clearh(metad);
  for ( t = 1; t <= nsteps; t++ ) {
    lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    if ( prod ) {
      if ( m->opta ) {
        metad->a = intq_interpa(metad->intq, (double) t);
      } else {
        metad->a = 1.0/(t + t0);
      }
    }
    //if ( prod ) {
    //  printf("t %ld, a %g\n", t, metad->a); getchar();
    //}
    if ( t % 10000 == 0 ) {
      *hfl = metad_hfl(metad, -1);
      fprintf(stderr, "%s t %8ld/%8ld = %5.2f%%, a %.3e, err %12.6e, hist. fl %12.6e/%5.3f%% %20s\r",
          (prod ? "prod." : "prep."), t, nsteps, 100.*t/nsteps,
          metad->a, metad_geterror(metad), (*hfl) * (*hfl), 100 * (*hfl), " ");
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
  int ir, kcerr;
  double dr, errtrunc;
  metad_t *metad;

  mtscramble(clock());

  lj = lj_open(np, rho, rcdef);
  lj_energy(lj);

  metad = metad_openf(0, lj->l*0.5, delr, m->pbc,
      METAD_SHIFT_TAIL, m->gaussig, m->kc, m->win, m->winn);
  /* load the reference bias potential */
  metad_load(metad, metad->vref, fnvref);
  errtrunc = metad_errtrunc(metad, metad->vref, &kcerr, "vtrunc.dat");
  ir = dist01(metad, lj, &dr);
  fprintf(stderr, "n %d, rmax %g, r %d/%g, err trunc %g, mode %d\n", metad->n, metad->xmax, ir, dr, errtrunc, kcerr);

  /* reduce the updating magnitude until it falls below m->alpha0 */
  //decmagrun(m, metad, lj);

  /* constant magnitude run */
  if ( m->gammethod != GAMMETHOD_NONE
    && m->gammethod != GAMMETHOD_LOAD ) {
    double hfl;
    prodrun(m, metad, lj, 0, m->nequil, "vg.dat", &hfl); /* equilibration */
    gammrun(m, metad, lj);
  }

  /* compute the optimal schedule and the error */
  metad_getalphaerr(metad, m->opta, (double) m->nsteps, m->gammethod,
      m->fngamma, m->sampmethod, m->alpha0, &m->qT,
      m->qprec, m->alpha_nint, m->fnalpha);

  /* multiple production runs */
  {
    int i, itr;
    ave_t ei[1], ef[1], fi[1], ff[1];
    double erri, errf, hfli, hflf; // , *v0;
    FILE *fplog;

    //xnew(v0, metad->n);
    //for ( i = 0; i < metad->n; i++ ) v0[i] = metad->v[i];
    ave_clear(ei);
    ave_clear(ef);
    ave_clear(fi);
    ave_clear(ff);
    fprintf(stderr, "starting production metadynamics run of %ld/%ld steps..., a %g, err %g\n", m->nequil, m->nsteps, metad->a, metad->errref);
    fplog = fopen(fnlog, "a");
    metad_saveheader(metad, fplog);
    fclose(fplog);
    for ( itr = 0; itr < m->ntrials; itr++ ) {
      for ( i = 0; i < metad->n; i++ ) metad->v[i] = 0; // v0[i];
      erri = prodrun(m, metad, lj, 0, m->nequil, "vi.dat", &hfli);
      ave_add(ei, erri);
      ave_add(fi, hfli);
      errf = prodrun(m, metad, lj, 1, m->nsteps, "vf.dat", &hflf);
      ave_add(ef, errf);
      ave_add(ff, hflf);
      printf("%4d: %14.10f %14.10f %14.10f %14.10f | %14.10f %14.10f %14.10f %14.10f\n",
          itr, errf, ef->ave, erri, ei->ave, hflf, ff->ave, hfli, fi->ave);
      fplog = fopen(fnlog, "a");
      fprintf(fplog, "%d %12.6f %12.6f %12.6f %12.6f\n", itr, errf, erri, hfli, hflf);
      fclose(fplog);
    }
    //free(v0);
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
  m->nequil = 1000000L;
  m->nsteps = 10000000L;
  m->ntrials = 1000;
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
