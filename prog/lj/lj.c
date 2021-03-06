#include "lj.h"
#include "../ave.h"
#include <time.h>

int np = 100;
double rho = 0.8;
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
int loadxerr = 0;

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
  fprintf(stderr, "  --loadxerr: load xerr from file, default %d\n", loadxerr);
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
  } else if ( strcmpfuzzy(key, "loadxerr") == 0 ) {
    loadxerr = 1;
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
static int gammrun(invtpar_t *m, metad_t *metad, lj_t *lj, long nsteps)
{
  long t;
  int ir0, ir;
  double dr, sacc = 0;

  metad->a = m->alpha0;

  if ( hblksz <= 0 ) {
    hblksz = (int)(0.1/metad->a);
    if ( hblksz <= 0 ) hblksz = 0;
  }

  ir0 = dist01(metad, lj, &dr);

  fprintf(stderr, "starting gamma run of %ld steps (updating magnitude %g)...%20s\n",
      nsteps, metad->a, "");
  metad_varv_clear(metad);
  metad_hblk_clear(metad);
  metad_clearav(metad); /* clear the average correction data */
  for ( t = 1; t <= nsteps; t++ ) {
    sacc += lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    metad_updatev(metad, ir);
    metad_vav_wupdate(metad, ir, metad->a);
    metad->tmat[ir*metad->n + ir0] += 1;
    metad_hblk_add(metad, ir);
    ir0 = ir;
    if ( t % m->gam_nstave == 0 ) {
      double *v = ( metad->winn == 1 ) ? metad->v : metad->vav;
      metad_varv_add(metad, v);
    }
    if ( t % hblksz == 0 ) {
      metad_hblk_dump(metad);
    }
    if ( t % 10000 == 0 ) fprintf(stderr, "t %ld/%ld = %5.2f%%, acc %.2f%% \r", t, nsteps, 100.*t/nsteps, 100*sacc/t);
    if ( t % nstsave == 0 || t == nsteps ) {
      metad_getgamma_varv(metad, m->alpha0, metad->gamma, "gamma.dat");
      metad_getgamma_tmat(metad, 1, metad->tgamma, "tgamma.dat");
      metad_getgamma_hblk(metad, metad->hgamma, "hgamma.dat");
      metad_save(metad, fnvbias);
    }
  }
  fprintf(stderr, "\n");
  return 0;
}


/* production run
 * `fixa` use a fixed updating magnitude
 * `*hfl` is the histogram fluctuation
 * `*erc` is the error of the average histogram-corrected bias potential */
static double prodrun(invtpar_t *m, metad_t *metad, lj_t *lj,
    int fixa, long nsteps, const char *fn, double *hfl, double *erc)
{
  int ir;
  long t;
  double t0, dr, err = 0;

  metad->a = m->alpha0;
  t0 = 2/metad->a;
  metad_clearh(metad);
  metad_clearav(metad); /* clear the average correction data */
  //if ( fixa )
  //  metad_varv_clear(metad);
  *erc = 0;
  for ( t = 1; t <= nsteps; t++ ) {
    lj_metroblk(lj, metad);
    ir = dist01(metad, lj, &dr);
    if ( !fixa ) { /* use the optimal schedule */
      if ( m->opta ) {
        //double a0 = intq_interpa(metad->intq, (double) t);
        metad->a = intq_evala(metad->intq, (double) t);
      } else {
        metad->a = 1.0/(t + t0);
      }
    }
    metad_updatev(metad, ir);
    metad_updateav(metad, ir); /* for the average histogram-corrected bias potential */
    //if ( fixa && t % 100 == 0 )
    //  metad_add_varv(metad);
    if ( t % nstsave == 0 || t == nsteps )
      metad_save(metad, fn);
    if ( t % 10000 == 0 || t == nsteps ) {
      *hfl = metad_hfl(metad, -1);
      err = metad_geterrav(metad, erc);
      fprintf(stderr, "%s t %8ld/%8ld = %5.2f%%, a %.3e, err %9.3e/%9.3e, "
                      "hist. fl %9.3e/%5.3f%% %g %20s\r",
          (fixa ? "fix-a" : "var-a"), t, nsteps, 100.*t/nsteps,
          metad->a, err, *erc, (*hfl) * (*hfl), 100 * (*hfl), metad->avcnt," ");
    }
  }
  //if ( fixa )
  //  metad_getgamma_varv(metad, m->alpha0, metad->tgamma, "gam.dat");
  *hfl *= *hfl;
  return err;
}

/* save the error components */
static int metad_save_cmvar(metad_t *metad,
    cmvar_t *cmi, cmvar_t *cmf, cmvar_t *cci, cmvar_t *ccf,
    const char *fn)
{
  int k, n = metad->n;
  double uk;
  FILE *fp;

  cmvar_get(cmi);
  cmvar_get(cmf);
  cmvar_get(cci);
  cmvar_get(ccf);
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  getcosmodes(metad->vref, n, metad->vft, metad->costab);
  metad_saveheader(metad, fp);
  fprintf(fp, "# %ld\n", cmi->cnt);
  for ( k = 1; k < n; k++ ) {
    uk = metad->vft[k];
    fprintf(fp, "%d %g %g %g %g %g %g %g %g %g %g %g\n", k,
        cmf->uvar[k], cmf->uave[k] - uk,
        cmi->uvar[k], cmi->uave[k] - uk,
        ccf->uvar[k], ccf->uave[k] - uk,
        cci->uvar[k], cci->uave[k] - uk,
        uk, metad->xefref[k], metad->xeiref[k]);
  }
  fclose(fp);
  return 0;
}

/* load the error components */
static int metad_load_cmvar(metad_t *metad,
    cmvar_t *cmi, cmvar_t *cmf, cmvar_t *cci, cmvar_t *ccf,
    const char *fn)
{
  int k, ki, cc, n = metad->n;
  long cnt;
  double uk;
  double *mfave, *miave, *cfave, *ciave;
  double *mfvar, *mivar, *cfvar, *civar;
  FILE *fp;
  char s[1024];

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  /* skip the information lines */
  fgets(s, sizeof s, fp);
  if ( s[0] != '#' ) goto ERR;
  k = atoi(s + 1);
  if ( k != n ) {
    fprintf(stderr, "Wrong number of bins %d vs %d(file)\n",
      n, k);
    goto ERR;
  }

  fgets(s, sizeof s, fp);
  if ( s[0] != '#' ) goto ERR;
  cnt = atol(s + 1);

  xnew(mfave, n); xnew(mfvar, n);
  xnew(miave, n); xnew(mivar, n);
  xnew(cfave, n); xnew(cfvar, n);
  xnew(ciave, n); xnew(civar, n);
  getcosmodes(metad->vref, n, metad->vft, metad->costab);
  for ( k = 1; k < n; k++ ) {
    uk = metad->vft[k];
    fgets(s, sizeof s, fp);
    cc = sscanf(s, "%d%lf%lf%lf%lf%lf%lf%lf%lf", &ki,
        &mfvar[k], &mfave[k], &mivar[k], &miave[k],
        &cfvar[k], &cfave[k], &civar[k], &ciave[k]);
    if ( cc != 9 || ki != k ) {
      fprintf(stderr, "error for mode %d/%d, cc %d, file %s\n",
          k, ki, cc, fn);
      goto ERR;
    }
    mfave[k] += uk;
    miave[k] += uk;
    cfave[k] += uk;
    ciave[k] += uk;
    //printf("var %g\n", mfvar[k]); getchar();
  }
  cmvar_set(cmi, cnt, miave, mivar);
  cmvar_set(cmf, cnt, mfave, mfvar);
  cmvar_set(cci, cnt, ciave, civar);
  cmvar_set(ccf, cnt, cfave, cfvar);
  free(mfave); free(mfvar);
  free(miave); free(mivar);
  free(cfave); free(cfvar);
  free(ciave); free(civar);
  fclose(fp);
  return 0;
ERR:
  fclose(fp);
  return -1;
}

static int work(invtpar_t *m)
{
  lj_t *lj;
  int ir, kcerr, n;
  double dr, errtrunc, hfl, errc;
  metad_t *metad;

  mtscramble(clock());

  lj = lj_open(np, rho, rcdef);
  lj_energy(lj);

  metad = metad_openf(0, lj->l*0.5, delr, m->pbc,
      METAD_SHIFT_TAIL, m->gaussig, m->kc, m->win, m->winn);
  n = metad->n;
  metad->a = m->alpha0;
  /* load the reference bias potential */
  metad_load(metad, metad->vref, fnvref);
  /* estimate the ``truncation error'' for the Gaussian updating scheme
   * does not alter vref */
  errtrunc = metad_errtrunc(metad, metad->vref, &kcerr, "vtrunc.dat");
  ir = dist01(metad, lj, &dr);
  fprintf(stderr, "n %d, rmax %g, r %d/%g, err trunc %g, mode %d\n", n, metad->xmax, ir, dr, errtrunc, kcerr);

  /* constant magnitude run, to activate, use the options
   * --gam=varv --gamnsteps=10000000 */
  if ( m->gammethod != GAMMETHOD_NONE
    && m->gammethod != GAMMETHOD_LOAD ) {
    prodrun(m, metad, lj, 1, m->nequil, "vg.dat", &hfl, &errc); /* equilibration */
    gammrun(m, metad, lj, m->gam_nsteps);
    ///* compute the optimal schedule based the histogram-corrected
    // * bias potential obtained from the gamma run */
    //metad_getalphaerr(metad, m->opta, (double) m->nsteps,
    //    m->gammethod, m->fngamma, m->sampmethod,
    //    metad->vcorr, m->alpha0, (double) m->nequil, &m->qT,
    //    m->qprec, m->alpha_nint, "alpha_vc.dat");
  }

  /* compute the optimal schedule and the error */
  metad_getalphaerr(metad, m->opta, (double) m->nsteps,
      m->gammethod, m->fngamma, m->sampmethod,
      metad->vref, m->alpha0, (double) m->nequil, &m->qT,
      m->qprec, m->alpha_nint, m->fnalpha);

  /* multiple production runs */
  if ( m->ntrials > 0 ) {
    int i, itr;
    ave_t ei[1], ef[1], eci[1], ecf[1], fi[1], ff[1];
    double erri, errf, erci, ercf, hfli, hflf, alpha0;
    cmvar_t *cmi, *cmf, *cci, *ccf;
    FILE *fplog;

    ave_clear(ei); ave_clear(eci);
    ave_clear(ef); ave_clear(ecf);
    ave_clear(fi);
    ave_clear(ff);
    cmi = cmvar_open(n, 0); /* initial */
    cmf = cmvar_open(n, 0); /* final */
    cci = cmvar_open(n, 0); /* initial corrected */
    ccf = cmvar_open(n, 0); /* final correct */
    if ( loadxerr ) {
      i = metad_load_cmvar(metad, cmi, cmf, cci, ccf, fnxerr);
      if ( i != 0 ) goto TRIAL_END;
      metad_save_cmvar(metad, cmi, cmf, cci, ccf, "new_xerr.dat");
      //printf("please check new_xerr.dat...\n"); getchar();
    } else {
      metad_save_cmvar(metad, cmi, cmf, cci, ccf, "xerr0.dat");
    }
    alpha0 = m->opta ? metad->intq->aarr[0] : metad->a/2;
    fprintf(stderr, "starting %ld metadynamics run of %ld/%ld steps..., a %g (%g), err %g -> %g, loadxerr %d\n",
        m->ntrials, m->nequil, m->nsteps, metad->a, alpha0, metad->eiref, metad->efref, loadxerr);
    /* write the header information for the log file */
    fplog = fopen(fnlog, "a");
    metad_saveheader(metad, fplog);
    fprintf(fplog, "# %ld %ld %ld %g %g %g\n",
        m->ntrials, m->nequil, m->nsteps, metad->a, metad->eiref, metad->efref);
    fclose(fplog);
    for ( itr = 0; itr < m->ntrials; itr++ ) {
      for ( i = 0; i < n; i++ ) metad->v[i] = 0; // v0[i];
      erri = prodrun(m, metad, lj, 1, m->nequil, "vi.dat", &hfli, &erci);
      ave_add(ei, erri);
      ave_add(fi, hfli);
      ave_add(eci, erci);
      cmvar_add(cmi, metad->v);
      cmvar_add(cci, metad->vcorr);
      errf = prodrun(m, metad, lj, 0, m->nsteps, "vf.dat", &hflf, &ercf);
      ave_add(ef, errf);
      ave_add(ff, hflf);
      ave_add(ecf, ercf);
      cmvar_add(cmf, metad->v);
      cmvar_add(ccf, metad->vcorr);
      printf("%4d: %9.7f %9.7f/%9.7f(%5.2f/%5.2f) %9.7f %9.7f/%9.7f(%5.2f/%5.2f) | %9.7f %9.7f %9.7f %9.7f\n",
          itr, errf, ef->ave, ecf->ave, ef->ave/metad->efref, ecf->ave/metad->efref,
               erri, ei->ave, eci->ave, ei->ave/metad->eiref, eci->ave/metad->eiref,
               hflf, ff->ave, hfli, fi->ave);
      /* save the log file */
      fplog = fopen(fnlog, "a");
      fprintf(fplog, "%d %.8f %.8f %.8f %.8f %.8f %.8f\n",
          itr, errf, erri, hflf, hfli, ercf, erci);
      fclose(fplog);
      /* save the components */
      metad_save_cmvar(metad, cmi, cmf, cci, ccf, fnxerr);
    }
TRIAL_END:
    cmvar_close(cmi);
    cmvar_close(cmf);
    cmvar_close(cci);
    cmvar_close(ccf);
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
