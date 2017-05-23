#include "metad.h"
#include "lj.h"
#include "ave.h"
#include <time.h>

int n = 100;
double rho = 0.8;
double rcdef = 100.0;
double tp = 3.0;
int mcblk = 5;
double delr = 0.005; /* spacing */
const char *fnvbias = "vbias.dat";
const char *fnvref = "vref.dat";

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

/* run with decreasing magnitude
 * until it falls under m->alpha0 */
static void decmagrun(invtpar_t *m, metad_t *metad, lj_t *lj)
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
    metad_add_varv(metad);
    if ( t % 10000 == 0 ) fprintf(stderr, "t %ld/%ld = %5.2f%%, acc %.2f%% \r", t, m->gam_nsteps, 100.*t/m->gam_nsteps, 100*sacc/t);
  }
  metad_save(metad, fnvbias);

  /* estimate gamma */
  metad_getgamma_varv(metad, m->alpha0, "gamma.dat");
  metad_getgamma_tmat(metad, 1, "tgamma.dat");
  return 0;
}


/* production run */
static double prodrun(invtpar_t *m, metad_t *metad, lj_t *lj,
    int prod, long nsteps)
{
  int ir;
  long t;
  double t0, dr;

  metad->a = m->alpha0;
  t0 = 2/metad->a;
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
      fprintf(stderr, "%s t %8ld/%8ld = %5.2f%%, a %.3e, err %g %20s\r",
          (prod ? "prod." : "prep."), t, nsteps, 100.*t/nsteps,
          metad->a, metad_geterror(metad), " ");
    }
    metad_updatev(metad, ir);
  }
  return metad_geterror(metad);
}

static int work(invtpar_t *m)
{
  lj_t *lj;
  int ir;
  double dr;
  metad_t *metad;

  mtscramble(clock());

  lj = lj_open(n, rho, rcdef);
  lj_energy(lj);

  metad = metad_openf(0, lj->l*0.5, delr,
      m->pbc, m->gaussig, m->kc, m->win, m->winn);
  ir = dist01(metad, lj, &dr);
  fprintf(stderr, "n %d, rmax %g, r %d/%g\n", metad->n, metad->xmax, ir, dr);

  /* reduce the updating magnitude until it falls below m->alpha0 */
  decmagrun(m, metad, lj);

  /* constant magnitude run */
  if ( m->gammethod != GAMMETHOD_NONE
    && m->gammethod != GAMMETHOD_LOAD )
    gammrun(m, metad, lj);

  /* compute the optimal schedule */
  if ( m->opta )
    metad_getalpha(metad, (double) m->nsteps, m->gammethod,
        m->fngamma, m->sampmethod, m->alpha0, &m->qT,
        m->qprec, m->alpha_nint, m->fnalpha);

  /* multiple production runs */
  {
    int i, itr;
    ave_t ei[1], ef[1];
    double erri, errf, *v0;

    metad_load(metad, metad->vref, fnvref);
    xnew(v0, metad->n);
    for ( i = 0; i < metad->n; i++ ) v0[i] = metad->v[i];
    ave_clear(ei);
    ave_clear(ef);
    fprintf(stderr, "starting production metadynamics run of %ld/%ld steps..., a %g, err %g\n", m->nequil, m->nsteps, metad->a, metad->errref);
    for ( itr = 0; itr < m->ntrials; itr++ ) {
      for ( i = 0; i < metad->n; i++ ) metad->v[i] = v0[i];
      erri = prodrun(m, metad, lj, 0, m->nequil);
      ave_add(ei, erri);
      metad_save(metad, "vi.dat");
      errf = prodrun(m, metad, lj, 1, m->nsteps);
      ave_add(ef, errf);
      metad_save(metad, "vf.dat");
      printf("%4d: %14g %14g %14g %14g %20s\n", itr, errf, ef->ave, erri, ei->ave, " ");
    }
    free(v0);
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
  m->gam_nsteps = 10000000L;
  m->nequil = 100000L;
  m->nsteps = 10000000L;
  m->alpha0 = 1e-4;
  m->pbc = 0;
  m->gammethod = GAMMETHOD_NONE;
  m->sampmethod = SAMPMETHOD_GAUSS;
  strcpy(m->fnalpha, "alpha.dat");
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  work(m);
  invtpar_finish(m);
  return 0;
}
