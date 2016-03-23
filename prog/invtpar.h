#ifndef INVTPAR_H__
#define INVTPAR_H__



/* parameters for the invt.c */




#include "util.h"



/* maximal number of neighbors
 * in the multiple-bin updating scheme */
#define NBMAX 1025



typedef struct {
  double c; /* constant for the updating magnitude */
  double t0; /* offset of equation c/(t + t0) */
  int n; /* number of bins */
  double *p; /* target distribution */
  double alpha0; /* initial updating magnitude */
  int fixa; /* fix alpha during the entire process */

  int optc; /* use the optimal c */

  int opta; /* use the analytically optimal alpha(t) */
  int alpha_nint; /* number of integration points for the exactly optimal alpha(t) */
  char fnalpha[FILENAME_MAX]; /* output file for the exactly optimal alpha */

  int pbc; /* periodic boundary condition */
  int winn; /* width of the updating window function */
  double win[NBMAX + 1]; /* shape of the window function */
  double wingaus; /* standard deviation of the Gaussian window */
  int winmax; /* explicit width for the Gaussian window */
  double initrand; /* magnitude of the initial Gaussian noise */
  int kcutoff; /* cutoff wave number of the initial noise */
  int sampmethod; /* sampling method */
  double localg; /* hopping probability for the local sampling process */
  double tcorr; /* correlation time for the sampling method */

  int pregamma; /* use simulation to estimate the gamma values */
  long gam_nsteps; /* number of steps */
  int gam_nstave; /* interval of accumulating averages */
  char fngamma[FILENAME_MAX]; /* file name for the gamma values */

  /* molecular dynamics parameters */
  double mddt; /* MD time step */
  double tp; /* temperature */
  double thermdt; /* thermostat time step */
#if 0
  double dwa; /* parameter a in sin(x) * (a - b * sin(x)) */
  double dwb; /* parameter b in sin(x) * (a - b * sin(x)) */
#endif

  int docorr; /* compute correlation functions */
  int nstcorr; /* time interval of computing the autocorrelation functions */
  double corrtol; /* error tolerance of the autocorrelation functions */
  char fncorr[FILENAME_MAX]; /* file name for the autocorrelation functions */

  long nequil; /* equilibration time */
  long nsteps; /* number of steps */
  long ntrials; /* number of trials */

#ifdef SCAN /* for predict.c */
  int cscan;
  double cmin;
  double cdel;
  double cmax;

  int nbscan;
  double nbmin;
  double nbdel;
  double nbmax;

  int sigscan;
  double sigmin;
  double sigdel;
  double sigmax;
#endif /* SCAN */

  int verbose; /* verbose level */
  const char *prog; /* name of the program */
} invtpar_t;



enum {
  SAMPMETHOD_METROGLOBAL = 0,
  SAMPMETHOD_METROLOCAL,
  SAMPMETHOD_HEATBATH,
  SAMPMETHOD_OU, /* Ornstein-Uhlenbeck process */
  SAMPMETHOD_MD,
  SAMPMETHOD_COUNT
};


#define MAX_OPT_ALIASES 8

const char *sampmethod_names[][MAX_OPT_ALIASES] = {
  {"global Metropolis", "global", "g"},
  {"local Metropolis", "local", "l"},
  {"heat-bath", "h"},
  {"Ornstein-Uhlenbeck", "OU", "o",
   "Harmonic oscillator", "HO"},
  {"molecular dynamics", "MD", "d"},
  {""}
};



/* initialize the default parameters */
static void invtpar_init(invtpar_t *m)
{
  int i;

  m->c = 1.0;
  m->t0 = 0;
  m->n = 100;
  m->p = NULL;
  m->alpha0 = 0.0;
  m->fixa = 0;
  m->optc = 0;

  m->opta = 0;
  m->alpha_nint = 1000;
  m->fnalpha[0] = '\0';;

  m->pbc = 0;
  m->winn = 1; /* single bin update */
  m->winmax = 0;
  for ( i = 1; i <= NBMAX; i++ ) {
    m->win[i] = 0;
  }
  m->win[0] = 1; /* single-bin case */
  m->wingaus = 0;
  m->initrand = 0;
  m->kcutoff = 0;
  m->sampmethod = 0;
  m->localg = 0.5;
  m->tcorr = 1.0; /* only used for the Ornstein-Uhlenbeck process */

  m->pregamma = 0;
  m->gam_nsteps = 100000000L;
  m->gam_nstave = 0;
  m->fngamma[0] = '\0';

  /* molecular dynamics parameters */
  m->mddt = 0.01;
  m->tp = 1.0;
  m->thermdt = m->mddt;
#if 0
  m->dwa = 1.0;
  m->dwb = 2.0;
#endif

  m->docorr = 0; /* don't compute correlation functions */
  m->nstcorr = 0; /* do correlation functions */
  m->corrtol = 1e-4;
  m->fncorr[0] = '\0';

  m->nequil = m->n * 10000L;
  m->nsteps = 100000000L;
  m->ntrials = 100;
  m->prog = "invt";
  m->verbose = 0;

#ifdef SCAN /* for predict.c */
  m->cscan = 0;
  m->cmin = 0.1;
  m->cdel = 0.01;
  m->cmax = 5.0;

  m->nbscan = 0;
  m->nbmin = 0.0;
  m->nbdel = 0.01;
  m->nbmax = 0.25;

  m->sigscan = 0;
  m->sigmin = 0.0;
  m->sigdel = 0.2;
  m->sigmax = 10.0;
#endif /* SCAN */
}



static void invtpar_mkgauswin(invtpar_t *m)
{
  int i;
  double x, s, sig = m->wingaus;

  m->winn = NBMAX;

  if ( m->winn > m->n ) {
    m->winn = m->n;
  }

  /* truncate the Gaussian at 5 * sigma */
  if ( m->winn >= sig * 5 ) {
    m->winn = (int) (sig * 5 + 0.5);
    if ( m->winn < 1 ) {
      m->winn = 1;
    }
  }

  /* further limit the window width
   * if an explicit value is given */
  if ( m->winmax > 0 && m->winn > m->winmax ) {
    m->winn = m->winmax;
  }

  m->win[0] = 1;
  s = m->win[0];
  for ( i = 1; i < m->winn; i++ ) {
    x = i / sig;
    m->win[i] = exp(-0.5 * x * x);
    s += m->win[i] * 2;
  }

  /* normalize the window function, such that
   * win[0] + 2 * (win[1] + ... + win[n - 1]) = 1 */
  for ( i = 0; i < m->winn; i++ ) {
    m->win[i] /= s;
  }
}



/* compute dependent parameters
 * call this function only once
 * the second call will miss supposedly default parameters */
static void invtpar_compute(invtpar_t *m)
{
  double x = 0;
  int i;

  if ( m->alpha0 <= 0 ) {
    /* alpha0 is to be divided by p[i] = 1/n
     * the coefficient can be as large as 1.0 for global MC moves
     * but should be much smaller for local MC moves
     * the current value is good for n up to 100 */
    m->alpha0 = 0.01 / m->n;
  }

  if ( fabs(m->t0) <= 0 ) { /* means t0 == 0 */
    /* set the default t0 such that
     * alpha is continuous at the beginning
     * of the production run */
    m->t0 = m->c / m->alpha0;
  }

  /* turn on correlation function computation
   * if the output file or frequency is set */
  if ( m->fncorr[0] != '\0' || m->nstcorr != 0 ) {
    m->docorr = 1;
  }

  if ( m->docorr ) {
    /* set default parameters for computing
     * autocorrelation function */
    if ( m->nstcorr <= 0 ) {
      m->nstcorr = (int) (1.0 / m->alpha0 + 0.5);
    }

    if ( m->fncorr[0] == '\0' ) {
      strcpy(m->fncorr, "corr.dat");
    }

    m->fixa = 1;
    m->ntrials = 1;
  }

  if ( m->fngamma[0] != '\0' ) {
    m->pregamma = 1;
  }

  if ( m->pregamma ) {
    if ( m->gam_nstave <= 0 ) {
      m->gam_nstave = m->n;
    }
  }

  /* initialize the window function */
  if ( m->wingaus > 0 ) {

    /* construct the Gaussian window */
    invtpar_mkgauswin(m);

  } else {

    /* normalize the window function
     * such that it sums to 1.0 */
    for ( x = 0, i = 1; i < m->winn; i++ ) {
      x += m->win[i];
    }
    m->win[0] = 1 - 2 * x;

  }

  /* initialize the target distribution */
  if ( m->p ) {
    /* in case this function has been called twice
     * release the previous memory */
    free(m->p);
  }
  /* initialize a normalized distribution */
  xnew(m->p, m->n);
  for ( i = 0; i < m->n; i++ ) {
    m->p[i] = 1.0 / m->n;
  }

  /* set the initial cutoff wave number */
  if ( m->kcutoff <= 0 ) {
    m->kcutoff = m->n;
  }

  /* currently, turn on PBC for MD automatically */
  if ( m->sampmethod == SAMPMETHOD_MD ) {
    m->pbc = 1;
  }

  if ( m->localg > 0.5 ) {
    fprintf(stderr, "local hopping probability cannot exceed 0.5\n");
    m->localg = 0.5;
  }
}



/* clean up */
static void invtpar_finish(invtpar_t *m)
{
  free(m->p);
  m->p = NULL;
}



/* print help message and die */
static void invtpar_help(const invtpar_t *m)
{
  fprintf(stderr, "Inverse-time formula for Wang-Landau and metadynamics simulations\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  %s [Options] [input.cfg]\n\n", m->prog);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  --n=:          set the number of bins, default %d\n", m->n);
  fprintf(stderr, "  --c=:          set c in alpha = c/(t + t0), default %g\n", m->c);
  fprintf(stderr, "  --a0=:         set the initial alpha during equilibration, default %g\n", m->alpha0);
  fprintf(stderr, "  --t0=:         set t0 in alpha = c/(t + t0), if unset, t0 = c/a0, default %g\n", m->t0);
  fprintf(stderr, "  --fixa:        fix the alpha during the entire process, default %d\n", m->fixa);
  fprintf(stderr, "  --optc:        use the optimal c for alpha(t), default %d\n", m->optc);
  fprintf(stderr, "  --opta:        use the exact optimal schedule alpha(t), default %d\n", m->opta);
  fprintf(stderr, "  --nint:        set the number of integration points for the exact optimal schedule alpha(t), default %d\n", m->alpha_nint);
  fprintf(stderr, "  --fnalpha=:    set the output file to output the exact optimal schedule, alpha(t), default %s\n", m->fnalpha);
  fprintf(stderr, "  --pbc:         use periodic boundary condition, default %d\n", m->pbc);
  fprintf(stderr, "  --nb=:         explicitly set the update window parameters, separated by comma, like --nb=0.2,0.4\n");
  fprintf(stderr, "  --sig=:        set the standard deviation Gaussian window, default %g\n", m->wingaus);
  fprintf(stderr, "  --wmax=:       set the width truncation of Gaussian window, default %d\n", m->winmax);
  fprintf(stderr, "  --initrand=:   magnitude of the initial random error, default %g\n", m->initrand);
  fprintf(stderr, "  --kcutoff=:    cutoff of wave number of the initial random error, default %d\n", m->kcutoff);
  fprintf(stderr, "  --samp=:       set the sampling scheme, g=global Metropolis, l=local Metropolis, h=heat-bath, d=molecular dynamics, o=Ornstein-Uhlenbeck, default %s\n", sampmethod_names[m->sampmethod][0]);
  fprintf(stderr, "  --localg=:     set the hopping probability for the local sampling process, [0, 0.5], default %g\n", m->localg);
  fprintf(stderr, "  --gam:         use a preliminary simulation to estimate the integrals of autocorrelation functions of the eigenmodes, default %d\n", m->pregamma);
  fprintf(stderr, "  --gamnsteps=:  set the number of steps in the preliminary simulation, default %ld\n", m->gam_nsteps);
  fprintf(stderr, "  --gamnstave=:  set the interval of accumulating data in the preliminary simulation, default %d\n", m->gam_nstave);
  fprintf(stderr, "  --fngamma=:    set the file for the gamma values, default %s\n", m->fngamma);
  fprintf(stderr, "  --corr:        compute correlation functions, default %d\n", m->docorr);
  fprintf(stderr, "  --nstcorr=:    set the number of steps of setting the correlation function, default %d\n", m->nstcorr);
  fprintf(stderr, "  --corrtol=:    set the tolerance level to truncate the autocorrelation function, default %g\n", m->corrtol);
  fprintf(stderr, "  --fncorr=:     set the file name for the correlation function, default %s\n", m->fncorr);
  fprintf(stderr, "  --try=:        set the number of trials, default %ld\n", m->ntrials);
  fprintf(stderr, "  --nsteps=:     set the number of simulation steps, default %ld\n", m->nsteps);
  fprintf(stderr, "  --equil=:      set the number of equilibration steps, default %ld\n", m->nequil);
#ifdef SCAN /* for predict.c */
  fprintf(stderr, "  --cmin=:       set the minimal c in c-scan, default %g\n", m->cmin);
  fprintf(stderr, "  --cmax=:       set the maximal c in c-scan, default %g\n", m->cmax);
  fprintf(stderr, "  --dc=:         set the increment of c in c-scan, default %g\n", m->cdel);
  fprintf(stderr, "  --nbmin=:      set the minimal nearest-neighbor updating magnitude nb in nb-scan, default %g\n", m->nbmin);
  fprintf(stderr, "  --nbmax=:      set the maximal nearest-neighbor updating magnitude nb in nb-scan, default %g\n", m->nbmax);
  fprintf(stderr, "  --dnb=:        set the increment of nb in nb-scan, default %g\n", m->nbdel);
  fprintf(stderr, "  --sigmin=:     set the minimal width of the Gaussian scheme in sig-scan, default %g\n", m->sigmin);
  fprintf(stderr, "  --sigmax=:     set the maximal width of the Gaussian scheme in sig-scan, default %g\n", m->sigmax);
  fprintf(stderr, "  --dsig=:       set the increment of the width in sig-scan, default %g\n", m->sigdel);
#endif /* SCAN */
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc.\n");
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
}



/* get integer */
static int invtpar_getint(invtpar_t *m,
    const char *key, const char *val)
{
  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    invtpar_help(m);
  }

  /* if 'e' exists in the string, scan the floating point
   * number, and then convert it to an integer */
  if ( strchr(val, 'e') != NULL
    || strchr(val, 'E') != NULL ) {
    return (int) ( atof(val) + 0.5 );
  } else {
    return atoi(val);
  }
}



/* get long integer */
static long invtpar_getlong(invtpar_t *m,
    const char *key, const char *val)
{
  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    invtpar_help(m);
  }

  /* if 'e' exists in the string, scan the floating point
   * number, and then convert it to a long integer */
  if ( strchr(val, 'e') ) {
    return (long) ( atof(val) + 0.5 );
  } else {
    return atol(val);
  }
}



/* get double */
static double invtpar_getdouble(invtpar_t *m,
    const char *key, const char *val)
{
  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    invtpar_help(m);
  }
  return atof(val);
}



/* parse a string `val` into an array */
static int invtpar_readarray(invtpar_t *m,
    const char *key, const char *val,
    double *arr, int n)
{
  const char *delims = ",;:";
  char *p, *s;
  int i = 0;

  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    invtpar_help(m);
  }

  /* make a copy of val */
  xnew(s, strlen(val) + 1);
  strcpy(s, val);

  p = strtok(s, delims);
  while ( p != NULL ) {
    arr[ i++ ] = atof(p);
    p = strtok(NULL, delims);
    if ( i >= n )
      break;
  }

  free(s);

  return i;
}



/* select an option */
static int invtpar_selectoption(invtpar_t *m,
    const char *key, const char *val,
    const char *names[][MAX_OPT_ALIASES], int cnt)
{
  int i, j;

  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    invtpar_help(m);
  }

  /* loop over options */
  for ( i = 0; i < cnt; i++ ) {
    /* loop over aliases */
    for ( j = 0; j < MAX_OPT_ALIASES; j++ ) {
      if ( names[i][j] == NULL
        || names[i][j][0] == '\0' ) {
        break;
      }

      if ( strcmpfuzzy(names[i][j], val) == 0 ) {
        return i;
      }
    }
  }

  /* try to treat `val` is a number */
  if ( isdigit(val[0]) ) {
    i = atoi(val);
    if ( i < 0 || i >= cnt ) i = 0;
  } else {
    i = 0;
  }

  return i;
}



/* get a boolean/integer value */
static int invtpar_getbool(invtpar_t *m,
    const char *key, const char *val)
{
  if ( val == NULL || val[0] == '\0' ) {
    return 1;
  }

  if ( isdigit(val[0]) ) {
    return atoi(val);
  }

  if ( strcmpfuzzy(val, "no") == 0
    || strcmpfuzzy(val, "false") == 0
    || strcmpfuzzy(val, "n") == 0
    || strcmpfuzzy(val, "f") == 0 ) {
    return 0;
  } else if ( strcmpfuzzy(val, "yes") == 0
           || strcmpfuzzy(val, "true") == 0
           || strcmpfuzzy(val, "y") == 0
           || strcmpfuzzy(val, "t") == 0 ) {
    return 1;
  } else {
    fprintf(stderr, "unknown value [%s] for %s\n", val, key);
    invtpar_help(m);
  }

  return 1;
}



/* match string key and value pairs */
static int invtpar_keymatch(invtpar_t *m,
    const char *key, const char *val)
{
  if ( strcmp(key, "n") == 0 )
  {
    m->n = invtpar_getint(m, key, val);
  }
  else if ( strcmpfuzzy(key, "c") == 0
         || strcmpfuzzy(key, "invt-c") == 0 )
  {
    m->c = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "t0") == 0 )
  {
    m->t0 = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "alpha0") == 0
         || strcmpfuzzy(key, "a0") == 0 )
  {
    m->alpha0 = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "fixa") == 0
         || strcmpfuzzy(key, "fixalpha") == 0 )
  {
    m->fixa = invtpar_getbool(m, key, val);
  }
  else if ( strcmpfuzzy(key, "optc") == 0 )
  {
    m->optc = invtpar_getbool(m, key, val);
  }
  else if ( strcmpfuzzy(key, "opta") == 0 )
  {
    m->opta = invtpar_getbool(m, key, val);
  }
  else if ( strcmpfuzzy(key, "nint") == 0 )
  {
    m->alpha_nint = invtpar_getint(m, key, val);
  }
  else if ( strcmpfuzzy(key, "fna") == 0
         || strcmpfuzzy(key, "fnalpha") == 0 )
  {
    strcpy(m->fnalpha, val);
  }
  else if ( strcmpfuzzy(key, "pbc") == 0 )
  {
    m->pbc = invtpar_getbool(m, key, val);
  }
  else if ( strcmpfuzzy(key, "nb") == 0
         || strstartswith(key, "neighbo")
         || strcmpfuzzy(key, "win") == 0
         || strcmpfuzzy(key, "window") == 0 )
  {
    m->winn = 1 + invtpar_readarray(m, key, val,
        m->win + 1, NBMAX);
  }
  else if ( strcmpfuzzy(key, "wgaus") == 0
         || strcmpfuzzy(key, "wingaus") == 0
         || strcmpfuzzy(key, "width-Gaussian") == 0
         || strcmpfuzzy(key, "sigma") == 0
         || strcmpfuzzy(key, "sig") == 0 )
  {
    m->wingaus = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "wmax") == 0
         || strcmpfuzzy(key, "w-cutoff") == 0 )
  {
    m->winmax = invtpar_getint(m, key, val);
  }
  else if ( strcmpfuzzy(key, "initrand") == 0 )
  {
    m->initrand = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "kcutoff") == 0 )
  {
    m->kcutoff = invtpar_getint(m, key, val);
  }
  else if ( strcmpfuzzy(key, "sampling-method") == 0
         || strcmpfuzzy(key, "sampmethod") == 0
         || strcmpfuzzy(key, "samp") == 0 )
  {
    m->sampmethod = invtpar_selectoption(m, key, val,
        sampmethod_names, SAMPMETHOD_COUNT);
  }
  else if ( strcmp(key, "g") == 0
         || strcmpfuzzy(key, "localg") == 0 )
  {
    m->localg = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "gam") == 0
         || strcmpfuzzy(key, "pre") == 0
         || strcmpfuzzy(key, "gamma") == 0 )
  {
    /* the value 1 means to compute
     * the value 2 means to load */
    if ( val != NULL ) {
      m->pregamma = invtpar_getint(m, key, val);
    } else {
      m->pregamma = 1;
    }
  }
  else if ( strcmpfuzzy(key, "gamnsteps") == 0
         || strcmpfuzzy(key, "gamsteps") == 0
         || strcmpfuzzy(key, "prensteps") == 0
         || strcmpfuzzy(key, "presteps") == 0 )
  {
    m->gam_nsteps = invtpar_getlong(m, key, val);
  }
  else if ( strcmpfuzzy(key, "gamnstave") == 0
         || strcmpfuzzy(key, "gamave") == 0
         || strcmpfuzzy(key, "prenstave") == 0
         || strcmpfuzzy(key, "preave") == 0 )
  {
    m->gam_nstave = invtpar_getint(m, key, val);
  }
  else if ( strcmpfuzzy(key, "fngamma") == 0
         || strcmpfuzzy(key, "fngam") == 0 )
  {
    strcpy(m->fngamma, val);
  }
  else if ( strcmpfuzzy(key, "tcorr") == 0
         || strcmpfuzzy(key, "corr-time") == 0 )
  {
    m->tcorr = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "mddt") == 0 )
  {
    m->mddt = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "tp") == 0
         || strstartswith(key, "temp") )
  {
    m->tp = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "thermdt") == 0 )
  {
    m->thermdt = invtpar_getdouble(m, key, val);
  }
#if 0
  else if ( strcmpfuzzy(key, "dwa") == 0 )
  {
    m->dwa = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "dwb") == 0 )
  {
    m->dwb = invtpar_getdouble(m, key, val);
  }
#endif
  else if ( strcmpfuzzy(key, "corr") == 0
         || strcmpfuzzy(key, "docorr") == 0 )
  {
    m->docorr = invtpar_getbool(m, key, val);
  }
  else if ( strcmpfuzzy(key, "nstcorr") == 0 )
  {
    m->nstcorr = invtpar_getint(m, key, val);
  }
  else if ( strcmpfuzzy(key, "corrtol") == 0 )
  {
    m->corrtol = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "fncorr") == 0 )
  {
    strcpy(m->fncorr, val);
  }
  else if ( strcmpfuzzy(key, "equil") == 0
           || strcmpfuzzy(key, "nequil") == 0 )
  {
    m->nequil = invtpar_getlong(m, key, val);
  }
  else if ( strcmpfuzzy(key, "steps") == 0
         || strcmpfuzzy(key, "nsteps") == 0 )
  {
    m->nsteps = invtpar_getlong(m, key, val);
  }
  else if ( strcmpfuzzy(key, "try") == 0
         || strstartswith(key, "repeat")
         || strstartswith(key, "trial")
         || strstartswith(key, "ntrial") )
  {
    m->ntrials = invtpar_getlong(m, key, val);
  }
  else if ( strcmpfuzzy(key, "verbose") == 0 )
  {
    m->verbose = invtpar_getbool(m, key, val);
  }
#ifdef SCAN /* for predict.c */
  /* c-scan paramerters */
  else if ( strcmpfuzzy(key, "cmin") == 0 )
  {
    m->cmin = invtpar_getdouble(m, key, val);
    m->cscan = 1;
  }
  else if ( strcmpfuzzy(key, "cdel") == 0
         || strcmpfuzzy(key, "dc") == 0 )
  {
    m->cdel = invtpar_getdouble(m, key, val);
    m->cscan = 1;
  }
  else if ( strcmpfuzzy(key, "cmax") == 0 )
  {
    m->cmax = invtpar_getdouble(m, key, val);
    m->cscan = 1;
  }
  /* nb-scan parameters */
  else if ( strcmpfuzzy(key, "nbmin") == 0 )
  {
    m->nbmin = invtpar_getdouble(m, key, val);
    m->nbscan = 1;
  }
  else if ( strcmpfuzzy(key, "nbdel") == 0
         || strcmpfuzzy(key, "dnb") == 0 )
  {
    m->nbdel = invtpar_getdouble(m, key, val);
    m->nbscan = 1;
  }
  else if ( strcmpfuzzy(key, "nbmax") == 0 )
  {
    m->nbmax = invtpar_getdouble(m, key, val);
    m->nbscan = 1;
  }
  /* sig-scan parameters */
  else if ( strcmpfuzzy(key, "sigmin") == 0 )
  {
    m->sigmin = invtpar_getdouble(m, key, val);
    m->sigscan = 1;
  }
  else if ( strcmpfuzzy(key, "sigdel") == 0
         || strcmpfuzzy(key, "dsig") == 0 )
  {
    m->sigdel = invtpar_getdouble(m, key, val);
    m->sigscan = 1;
  }
  else if ( strcmpfuzzy(key, "sigmax") == 0 )
  {
    m->sigmax = invtpar_getdouble(m, key, val);
    m->sigscan = 1;
  }
#endif /* SCAN */
  else
  {
    return -1;
  }

  return 0;
}



/* loading parameters from a configuration file
 * each line of the configuration file is
 *  key = val
 * */
static int invtpar_load(invtpar_t *m, const char *fn)
{
  FILE *fp;
  char buf[800], *p, *key, *val;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  while ( fgets(buf, sizeof buf, fp) ) {
    /* remove leading and trailing spaces */
    strstrip(buf);

    /* skip a blank or comment line */
    if ( buf[0] == '\0' || buf[0] == '#' ) {
      continue;
    }

    /* parse the line to a key and a value */
    key = buf;
    /* find the end of the key */
    if ( (p = strchr(key, '=')) != NULL ) {
      *p = '\0'; /* end the key part */
      strstrip(key);

      /* find the beginning of the value */
      for ( val = p + 1; isspace( *val ) ; )
        val++;
      strstrip(val);

      /* remove the trailing semicolon */
      p = val + strlen(val) - 1;
      if ( *p == ';' ) *p = '\0';
      strstrip(val);

    } else {
      val = NULL;
    }

    if ( invtpar_keymatch(m, key, val) != 0 ) {
      fprintf(stderr, "Warning: unknown options %s = %s in %s\n",
          key, val, fn);
    }
  }

  fclose(fp);

  /* do not call invtpar_compute(m) now
   * because it will call the function prematurally */
  /* invtpar_compute(m); */

  return 0;
}



/* scan parameters from command-line arguments */
static int invtpar_doargs(invtpar_t *m, int argc, char **argv)
{
  int i, j, ch;
  char *p, *q;

  for ( i = 1; i < argc; i++ ) {
    /* test if it is an argument */
    if ( argv[i][0] != '-' ) {
      invtpar_load(m, argv[i]);
      continue;
    }

    /* long argument, like --help */
    if ( argv[i][1] == '-' ) {
      /* try to parse the argment
         e.g., `--prog=aaa' is parsed to `--prog' and `aaa' */
      p = argv[i] + 2;
      /* let q point to the argument of the option */
      if ( (q = strchr(p, '=')) != NULL ) {
        *q++ = '\0';
      } else {
        q = NULL;
      }

      if ( invtpar_keymatch(m, p, q) != 0 ) {
        if ( strcmpfuzzy(p, "help") != 0 ) {
          fprintf(stderr, "Unknown option %s, key [%s], val [%s]\n",
              argv[i], p, (q != NULL ? q : "NULL") );
        }
        invtpar_help(m);
      }
      continue;
    }

    /* it is a short option
     * loop over characters in the options
     * in this way, `-vo' is understood as `-v -o' */
    for ( j = 1; (ch = argv[i][j]) != '\0'; j++ ) {
      if ( strchr("nc", ch) != NULL ) {
        /* handle options that require an argument */
        q = p = argv[i] + j + 1;
        if ( *p != '\0' ) {
          /* the argument follows immediately after the option
           * e.g., -oa.dat */
          q = p;
        } else if ( ++i < argc ) {
          /* the option and argument are separated by a space
           * then the argument belongs to the next argv[] element,
           * hence ++i
           * e.g., -o a.dat */
          q = argv[i];
        } else {
          fprintf(stderr, "-%c requires an argument!\n", ch);
          invtpar_help(m);
        }

        if ( ch == 'n' ) { /* the number of bins */
          m->n = atoi(q);
        } else if ( ch == 'c' ) { /* constant for inverse time */
          m->c = atof(q);
        }
        break; /* skip the rest of the characters in the option */
      } else if ( ch == 'v' ) {
        m->verbose++;
      } else if ( ch == 'h' ) {
        invtpar_help(m);
      } else {
        fprintf(stderr, "unknown option %s, j %d, ch %c\n", argv[i], j, ch);
        invtpar_help(m);
      }
    }
  }

  invtpar_compute(m);

  return 0;
}



/* print out the key parameters */
static void invtpar_dump(const invtpar_t *m)
{
  int i;
  double sum = 0;

  fprintf(stderr, "%ld trials: n %d, alpha = %g/(t + %g), alpha0 %g, "
      "pbc %d, %s, %ld steps; equil %ld steps\n",
      m->ntrials, m->n, m->c, m->t0, m->alpha0, m->pbc,
      sampmethod_names[m->sampmethod][0], m->nsteps, m->nequil);

  if ( m->pregamma ) {
    fprintf(stderr, "preliminary run: %ld steps, averaging every %d steps\n",
      m->gam_nsteps, m->gam_nstave);
  }

  fprintf(stderr, "update window function (%d bins): ", m->winn);
  for ( i = 0; i < m->winn; i++ ) {
    sum += m->win[i];
    fprintf(stderr, "%g ", m->win[i]);
  }
  sum = sum * 2 - m->win[0];
  fprintf(stderr, "| sum %g\n", sum);


#ifdef SCAN
  if ( m->cscan ) {
    fprintf(stderr, "c-scan window: %g:%g:%g\n",
        m->cmin, m->cdel, m->cmax);
  }

  if ( m->nbscan ) {
    fprintf(stderr, "nb-scan window: %g:%g:%g\n",
        m->nbmin, m->nbdel, m->nbmax);
  }

  if ( m->sigscan ) {
    fprintf(stderr, "sig-scan window: %g:%g:%g\n",
        m->sigmin, m->sigdel, m->sigmax);
  }
#endif /* SCAN */
}



#endif /* INVTPAR__ */


