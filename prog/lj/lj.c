#include "metad.h"
#include "lj.h"
#include "ave.h"
#include <time.h>

int n = 108;
double rho = 0.8;
double rcdef = 100.0;
double tp = 1.5;
int mcblk = 5;

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
  int ir, it;
  long t;
  double dx[D], l = lj->l, invl = 1/l, dr, amp = 0.1/rho;

  ir = dist01(metad, lj, &dr);
  for ( t = 1; ; t++ ) {
    lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    metad_updatev(metad, ir);
    if ( t % 1000 == 0 ) {
      int sacc = metad_wlcheck(metad, m->flatness, m->magred);
      if ( sacc && metad->a < m->alpha0 )
        break;
    }
    if ( t % 1000000 == 0 ) {
      printf("t %ld, fl %g\n", t, metad->hflatness);
      metad_save(metad, "vbias.dat");
    }
  }
}


/* constant updating magnitude run
 * to compute the gamma values */
static int cmagrun(invtpar_t *m, metad_t *metad, lj_t *lj)
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
    if ( t % 1000 == 0 ) fprintf(stderr, "t %ld/%ld, acc %g%% \r", t, m->gam_nsteps, 100*sacc/t);
  }

  /* estimate gamma */
  metad_getgamma_varv(metad, m->alpha0, "gamma.dat");
  metad_getgamma_tmat(metad, 1, "tgamma.dat");
  return 0;
}


/* production run */
static double prodrun(invtpar_t *m, metad_t *metad, lj_t *lj)
{
  int ir;
  long t;
  double t0, dr;

  metad->a = m->alpha0;
  t0 = 2/metad->a;
  for ( t = 1; t <= m->nsteps; t++ ) {
    lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    if ( m->opta ) {
      metad->a = intq_interpa(metad->intq, (double) t);
    } else {
      metad->a = 1.0/(t + t0);
    }
    metad_updatev(metad, ir);
  }
  return metad_geterror(metad);
}


static int work(invtpar_t *m)
{
  lj_t *lj;
  long t;
  int it, ir, acc;
  double dx[D], dr, l, invl, amp;
  double sacc = 0, xdel = 0.02;
  metad_t *metad;

  mtscramble(clock());

  lj = lj_open(n, rho, rcdef);
  lj_energy(lj);
  invl = 1/(l = lj->l);
  amp = 0.1/rho;

  metad = metad_openf(0, lj->l*0.5, xdel,
      m->pbc, m->gaussig, m->okmax, m->win, m->winn);
  ir = dist01(metad, lj, &dr);
  fprintf(stderr, "n %d, rmax %g, r %d/%g\n", metad->n, metad->xmax, ir, dr);

  /* reduce the updating magnitude until it falls below m->alpha0 */
  decmagrun(m, metad, lj);

  /* constant magnitude run */
  if ( m->gammethod != GAMMETHOD_NONE )
    cmagrun(m, metad, lj);

  /* compute the optimal schedule */
  if ( m->opta )
    metad_getalpha(metad, (double) m->nsteps, m->gammethod,
        m->alpha0, &m->qT, m->qprec, m->alpha_nint, m->fnalpha);

  /* multiple production runs */
  {
    int i, itr, n = metad->n;
    ave_t ef[1];
    double erri, errf, *v0;

    xnew(v0, n);
    for ( i = 0; i < n; i++ ) v0[i] = metad->v[i];
    ave_clear(ef);
    erri = metad_geterror(metad);
    metad_save(metad, "vi.dat");
    fprintf(stderr, "starting production metadynamics run of %ld steps..., a %g\n", m->gam_nsteps, metad->a);
    for ( itr = 0; itr < m->ntrials; itr++ ) {
      for ( i = 0; i < n; i++ ) metad->v[i] = v0[i];
      errf = prodrun(m, metad, lj);
      metad_save(metad, "vf.dat");
      ave_add(ef, errf);
      printf("%4ld: %14g %14g %14g\n", itr, errf, ef->ave, erri);
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
  m->nsteps = 100000000L;
  m->alpha0 = 1e-5;
  m->pbc = 0;
  m->pregamma = 1;
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  work(m);
  invtpar_finish(m);
  return 0;
}
