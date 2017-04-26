#include "invt.h"
#define IS2_LB 4
#include "is2.h"
#include "metad.h"

const double IS2_TC = 2.3;
const int IS2_EMIN = -2*IS2_N + 8;
const int IS2_EMAX = 0;

typedef invtpar_t par_t;

/* constant updating magnitude run */
static int invt_is2_prerun(par_t *m, metad_t *metad, is2_t *is)
{
  long t;
  int id, h, acc;
  int icur, inew, enew, iold;
  cmvar_t *cm;

  cm = cmvar_open(metad->n, 0);
  metad->a = m->alpha0;
  m->gam_nstave = 100;
  //strcpy(m->fncorr, "corr.dat");
  //metad->corr = corr_open(metad->n - 1, m->nsteps / m->nstcorr);
  iold = icur = metad_getindex(metad, is->E);

  fprintf(stderr, "starting constant magnitude metadynamics run of %ld steps...\n", m->gam_nsteps);
  for ( t = 1; t <= m->gam_nsteps; t++ ) {
    IS2_PICK(is, id, h);
    enew = is->E + h * 2;
    acc = metad_acc(metad, icur, enew, &inew);
    if ( acc ) {
      icur = inew;
      IS2_FLIP(is, id, h);
    }
    metad_updatev(metad, icur);
    if ( t % m->gam_nstave == 0 ) {
      metad->tmat[icur*metad->n + iold] += 1;
      iold = icur;
      //getcosmodes(icur, metad->n, metad->u, metad->costab);
      ///* the first mode should always be zero,
      // * so we start from the second mode, u + 1 */
      //corr_add(metad->corr, metad->u + 1);
      cmvar_add(cm, metad->v);
    }
  }

  /* compute gamma */
  {
    int i;
    double *tgamma;
    xnew(tgamma, cm->n);
    metad_getgamma_tmat(metad, tgamma, m->gam_nstave);
    cmvar_get(cm);
    for ( i = 1; i < metad->n; i++ ) {
      printf("%4d: %g %g\n", i, tgamma[i], cm->uvar[i]*2/m->alpha0);
    }
    free(tgamma);
  }
  //if ( metad->corr != NULL ) {
  //  /* print out the thermodynamic fluctuations */
  //  corr_printfluc(metad->corr, 0, NULL);
  //  /* compute the correlation functions and save them to file
  //   * the maximal span, 100/alpha, should be large enough */
  //  corr_save(metad->corr, m->nstcorr, 200000,
  //      m->corrtol, 0, m->fncorr);
  //}
  cmvar_close(cm);

  return 0;
}

static int invt_is2_run(par_t *m)
{
  is2_t *is;
  int id, h;
  long t;
  int enew, icur, inew, acc;
  metad_t *metad;

  mtscramble( clock() );

  is = is2_open(IS2_L);

  /* equilibration at the critical temperature */
  is2_setuproba(1.0/IS2_TC, is->uproba);
  for ( t = 1; t <= m->nequil; t++ ) {
    IS2_PICK(is, id, h);
    if ( h <= 0 || mtrand() <= is->uproba[h] ) {
      IS2_FLIP(is, id, h);
    }
  }

  /* typical WL run */
  metad = metad_open(IS2_EMIN, IS2_EMAX, 4, 0);
  icur = metad_getindex(metad, is->E);
  for ( t = 1; ; t++ ) {
    IS2_PICK(is, id, h);
    enew = is->E + h * 2;
    acc = metad_acc(metad, icur, enew, &inew);
    if ( acc ) {
      icur = inew;
      IS2_FLIP(is, id, h);
    }
    metad_updatev(metad, icur);
    if ( t % 1000 == 0 ) {
      int sacc = metad_wlcheck(metad);
      if ( sacc && metad->a < m->alpha0 )
        break;
    }
  }

  /* constant magnitude run */
  invt_is2_prerun(m, metad, is);

  metad_save(metad, "vbias.dat");
  is2_close(is);
  metad_close(metad);

  return 0;
}

int main(int argc, char **argv)
{
  par_t m[1];

  invtpar_init(m);
  m->nequil = 10000;
  m->nsteps = 100000000;
  m->alpha0 = 1e-6;
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  invt_is2_run(m);
  invtpar_finish(m);
  return 0;
}
