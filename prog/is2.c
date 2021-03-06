#include "invt.h"
#ifndef IS2_LB
#define IS2_LB 6
#endif
#include "is2.h"
#include "metad.h"
#include "cmvar.h"
#include "ave.h"

int l = 0; /* explicit value */

static void is2par_help(void)
{
  fprintf(stderr, "  --l=:          side length of Ising model, default %d\n", l);
}

static int is2par_keymatch(invtpar_t *m, const char *key, const char *val)
{
  if ( strcmp(key, "l") == 0 ) {
    l = invtpar_getint(m, key, val);
  } else {
    return -1;
  }
  return 0;
}

/* run with decreasing magnitude */
static void decmagrun(invtpar_t *m, metad_t *metad, is2_t *is)
{
  int id, h;
  int icur, inew, enew, acc;
  long t;

  icur = metad_getindex(metad, is->E);
  for ( t = 1; ; t++ ) {
    //IS2_PICK(is, id, h);
    id = is2_pick(is, &h);
    enew = is->E + h * 2;
    acc = metad_acc(metad, icur, enew, &inew);
    if ( acc ) {
      icur = inew;
      //IS2_FLIP(is, id, h);
      is2_flip(is, id, h);
    }
    metad_updatev(metad, icur);
    if ( t % 1000 == 0 ) {
      int sacc = metad_wlcheck(metad, m->fluc, m->magred);
      if ( sacc && metad->a < m->alpha0 )
        break;
    }
    if ( t % 1000000 == 0 ) {
      printf("t %ld, fl %g\n", t, metad->hfl);
      metad_save(metad, "vbias.dat");
    }
  }
}

/* constant updating magnitude run
 * computes the gamma values */
static int gammrun(invtpar_t *m, metad_t *metad, is2_t *is)
{
  long t;
  int id, h, acc;
  int icur, inew, enew, iold;

  metad_varv_clear(metad);
  metad->a = m->alpha0;
  iold = icur = metad_getindex(metad, is->E);

  fprintf(stderr, "starting constant magnitude %g metadynamics run of %ld steps...\n",
      metad->a, m->gam_nsteps);
  for ( t = 1; t <= m->gam_nsteps; t++ ) {
    //IS2_PICK(is, id, h);
    id = is2_pick(is, &h);
    enew = is->E + h * 2;
    acc = metad_acc(metad, icur, enew, &inew);
    if ( acc ) {
      icur = inew;
      //IS2_FLIP(is, id, h);
      is2_flip(is, id, h);
    }
    metad_updatev(metad, icur);
    metad->tmat[icur*metad->n + iold] += 1;
    iold = icur;
    if ( t % m->gam_nstave == 0 )
      metad_varv_add(metad, metad->v);
  }

  /* estimate gamma */
  metad_getgamma_varv(metad, m->alpha0, metad->gamma, "gamma.dat");
  metad_getgamma_tmat(metad, 1, metad->tgamma, "tgamma.dat");

  return 0;
}

/* production run */
static double prodrun(invtpar_t *m, metad_t *metad, is2_t *is, int prod, long nsteps)
{
  long t;
  int id, h, acc;
  int icur, inew, enew;
  double t0;

  metad->winn = 1; // make sure it is single-bin
  metad->a = m->alpha0;
  t0 = 2/metad->a;
  icur = metad_getindex(metad, is->E);
  for ( t = 1; t <= nsteps; t++ ) {
    //IS2_PICK(is, id, h);
    id = is2_pick(is, &h);
    enew = is->E + h * 2;
    acc = metad_acc(metad, icur, enew, &inew);
    if ( acc ) {
      icur = inew;
      //IS2_FLIP(is, id, h);
      is2_flip(is, id, h);
    }
    if ( prod ) {
      if ( m->opta ) {
        metad->a = intq_interpa(metad->intq, (double) t);
      } else {
        metad->a = 1.0/(t + t0);
      }
    }
    metad_updatev(metad, icur);
  }
  return metad_geterror(metad);
}

/* compute or load the exact density of states */
static int addvref(metad_t *metad)
{
  char fn[128], s[128];
  FILE *fp;
  double *lndos, v0;
  int i, id, n = l * l;

  sprintf(fn, "is2dos/is2lndos%dx%d.dat", l, l);
  if ((fp = fopen(fn, "r")) != NULL) {
    xnew(lndos, n + 1);
    for ( i = 0; i <= n; i++ ) {
      fgets(s, sizeof s, fp);
      sscanf(s, "%d%lf", &id, &lndos[i]);
    }
    fclose(fp);
  } else {
    fprintf(stderr, "cannot read %s, computing the exact density of states... ", fn);
    lndos = is2dos(l, l);
    is2dos_save(lndos, l, l, "is2dos/");
    fprintf(stderr, "done.\n");
  }
  fprintf(stderr, "lndos %g, %g, %g, ... %g\n", lndos[0], lndos[1], lndos[2], lndos[n]);

  /* map the exact density of states to the grid */
  for ( i = 0; i < metad->n; i++ ) {
    id = (metad->imin + 2*n)/4 + i;
    if ( i == 0 ) v0 = lndos[id];
    metad->vref[i] = lndos[id] - v0;
  }
  metad_shiftv(metad, metad->vref);
  free(lndos);

  return 0;
}

static int invt_is2_run(invtpar_t *m)
{
  is2_t *is;
  int id, h;
  long t, itr;
  metad_t *metad;

  int emax = 0; // -IS2_N;
  int emin = emax - (m->n - 1) * 4;
  int kcerr;

  //mtscramble( clock() );

  is = is2_open(l);

  /* equilibration at the infinite temperature */
  for ( t = 1; ; t++ ) {
    id = is2_pick(is, &h);
    is2_flip(is, id, h);
    if ( is->E >= emin && is->E < emax ) break;
  }

  metad = metad_open(emin, emax, 4,
      m->pbc, METAD_SHIFT_TAIL, m->gaussig, m->kc, m->win, m->winn);
  addvref(metad);
  metad_errtrunc(metad, metad->vref, &kcerr, "vtrunc.dat");

  /* gradually reduce the updating magnitude */
  decmagrun(m, metad, is);

  /* run with constant updating magnitude */
  if ( m->gammethod != GAMMETHOD_NONE
    && m->gammethod != GAMMETHOD_LOAD )
    gammrun(m, metad, is);

  /* compute the optimal schedule */
  metad_getalphaerr(metad, m->opta, (double) m->nsteps,
      m->gammethod, m->fngamma, m->sampmethod,
      metad->vref, m->alpha0, (double) m->nequil, &m->qT,
      m->qprec, m->alpha_nint, m->fnalpha);

  /* multiple production runs */
  if ( m->ntrials > 0 ) {
    int i;
    ave_t ei[1], ef[1];
    double erri, errf, *v0;

    xnew(v0, metad->n);
    for ( i = 0; i < metad->n; i++ ) v0[i] = metad->v[i];
    ave_clear(ei);
    ave_clear(ef);
    fprintf(stderr, "starting production metadynamics run of %ld/%ld steps..., a %g\n", m->nequil, m->nsteps, metad->a);
    for ( itr = 0; itr < m->ntrials; itr++ ) {
      for ( i = 0; i < metad->n; i++ ) metad->v[i] = v0[i];
      erri = prodrun(m, metad, is, 0, m->nequil);
      ave_add(ei, erri);
      metad_save(metad, "vi.dat");
      errf = prodrun(m, metad, is, 1, m->nsteps);
      ave_add(ef, errf);
      metad_save(metad, "vf.dat");
      printf("%4ld: %14g %14g %14g %14g %10s\n", itr, errf, ef->ave, erri, ei->ave, " ");
    }
    free(v0);
  }

  is2_close(is);
  metad_close(metad);

  return 0;
}

int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  m->nequil = 10000;
  m->nsteps = 100000000;
  m->alpha0 = 1e-6;
  m->gam_nsteps = 100000000;
  m->gam_nstave = 100;
  m->gammethod = GAMMETHOD_TMAT;
  m->pbc = 0;
  m->userhelp = is2par_help;
  m->usermatch = is2par_keymatch;
  invtpar_doargs(m, argc, argv);
  if ( l <= 0 ) l = IS2_L;
  invtpar_dump(m);
  printf("l %d\n", l);
  invt_is2_run(m);
  invtpar_finish(m);
  return 0;
}
