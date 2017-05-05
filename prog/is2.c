#include "invt.h"
#define IS2_LB 6
#include "is2.h"
#include "metad.h"
#include "cmvar.h"
#include "ave.h"



/* constant updating magnitude run
 * computes the gamma values */
static int invt_is2_prerun(invtpar_t *m, metad_t *metad, is2_t *is)
{
  long t;
  int id, h, acc;
  int icur, inew, enew, iold;
  cmvar_t *cm;

  cm = cmvar_open(metad->n, metad->pbc);
  metad->a = m->alpha0;
  m->gam_nstave = 100;
  //strcpy(m->fncorr, "corr.dat");
  //metad->corr = corr_open(metad->n - 1, m->nsteps / m->nstcorr);
  iold = icur = metad_getindex(metad, is->E);

  fprintf(stderr, "starting constant magnitude %g metadynamics run of %ld steps...\n",
      metad->a, m->gam_nsteps);
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

  /* estimate gamma */
  {
    int i, n = metad->n;

    cmvar_get(cm);
    for ( i = 1; i < n; i++ ) {
      double lam = metad->lambda[i];
      metad->gamma[i] = (lam > 0.1) ? cm->uvar[i]*2/m->alpha0/lam : 1;
    }

    metad_getgamma_tmat(metad, metad->gamma_tmat, m->gam_nstave);
    for ( i = 1; i < n; i++ ) {
      printf("%4d: %14g %14g %14g\n", i,
          metad->gamma[i], metad->gamma_tmat[i], metad->lambda[i]);
    }
    savegamma(n, metad->gamma, "gamma.dat");
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

/* production run */
static double invt_is2_simul(invtpar_t *m, metad_t *metad, is2_t *is)
{
  long t;
  int id, h, acc;
  int icur, inew, enew;
  double t0;

  metad->a = m->alpha0;
  t0 = 2/metad->a;
  icur = metad_getindex(metad, is->E);
  for ( t = 1; t <= m->nsteps; t++ ) {
    IS2_PICK(is, id, h);
    enew = is->E + h * 2;
    acc = metad_acc(metad, icur, enew, &inew);
    if ( acc ) {
      icur = inew;
      IS2_FLIP(is, id, h);
    }
    if ( m->opta ) {
      metad->a = intq_interpa(metad->intq, (double) t);
    } else {
      metad->a = 1.0/(t + t0);
    }
    metad_updatev_wl(metad, icur);
  }
  return metad_geterror(metad);
}

static int invt_is2_addexact(metad_t *metad, int l)
{
  char fn[128], s[128];
  FILE *fp;
  double *lndos, v0;
  int i, id, n;

  sprintf(fn, "is2dos/is2lndos%dx%d.txt", l, l);
  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fn);
    return -1;
  }
  n = l*l;
  xnew(lndos, n + 1);
  for ( i = 0; i <= n; i++ ) {
    fgets(s, sizeof s, fp);
    sscanf(s, "%d%lf", &id, &lndos[i]);
  }
  fclose(fp);

  /* map the exact dos to the grid */
  for ( i = 0; i < metad->n; i++ ) {
    id = (metad->xmin + 2*n)/4 + i;
    if ( i == 0 ) v0 = lndos[id];
    metad->vref[i] = lndos[id] - v0;
  }
  metad_trimv(metad, metad->vref);
  free(lndos);

  return 0;
}

static int invt_is2_run(invtpar_t *m)
{
  is2_t *is;
  int id, h;
  long t, itr;
  int enew, icur, inew, acc;
  metad_t *metad;

  int emax = 0; // -IS2_N;
  int emin = emax - (m->n - 1) * 4;

  //mtscramble( clock() );

  is = is2_open(IS2_L);

  /* equilibration at the infinite temperature */
  for ( t = 1; ; t++ ) {
    IS2_PICK(is, id, h);
    IS2_FLIP(is, id, h);
    if ( is->E >= emin && is->E < emax ) break;
  }

  /* typical WL run */
  metad = metad_open(emin, emax, 4,
      m->pbc, m->gaussig, m->okmax, m->win, m->winn);
  invt_is2_addexact(metad, IS2_L);
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
      int sacc = metad_wlcheck(metad, m->flatness, m->magred);
      if ( sacc && metad->a < m->alpha0 )
        break;
    }
    if ( t % 1000000 == 0 ) {
      printf("t %ld, fl %g\n", t, metad->hflatness);
      metad_save(metad, "vbias.dat");
    }
  }

  /* constant magnitude run */
  invt_is2_prerun(m, metad, is);

  {
    int i, n = metad->n;
    ave_t ef[1];
    double erri, errf, *v0;
    xnew(v0, n);
    for ( i = 0; i < n; i++ ) v0[i] = metad->v[i];
    ave_clear(ef);
    /* compute the optimal schedule */
    if ( m->opta ) metad_getalpha(metad, m, (double) m->nsteps);
    erri = metad_geterror(metad);
    metad_save(metad, "vi.dat");
    fprintf(stderr, "starting production metadynamics run of %ld steps..., a %g\n", m->gam_nsteps, metad->a);
    for ( itr = 0; itr < m->ntrials; itr++ ) {
      for ( i = 0; i < n; i++ ) metad->v[i] = v0[i];
      errf = invt_is2_simul(m, metad, is);
      metad_save(metad, "vf.dat");
      ave_add(ef, errf);
      printf("%4ld: %14g %14g %14g\n", itr, errf, ef->ave, erri);
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
  m->gam_nstave = 1000;
  m->pbc = 0;
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  invt_is2_run(m);
  invtpar_finish(m);
  return 0;
}
