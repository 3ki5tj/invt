#ifndef INVTPAR_H__
#define INVTPAR_H__



/* parameters for the invt.c */




#include "util.h"



/* maximal number of neighbors
 * in the multiple-bin updating scheme */
#define NBMAX 127



typedef struct {
  double c; /* constant for the updating magnitude */
  double t0; /* offset of equation c/(t + t0) */
  int n; /* number of bins */
  double *p; /* target distribution */
  double alpha0; /* initial updating magnitude */
  int fixa; /* fix alpha during the entire process */
  int winn; /* width of the updating window function */
  double win[NBMAX + 1]; /* shape of the window function */
  double wingaus; /* width of the Gaussian window */
  double initrand; /* magnitude of the initial Gaussian noise */
  int kcutoff; /* cutoff wave number of the initial noise */
  int sampmethod; /* sampling method */
  double tcorr; /* correlation time for the sampling method */

  int docorr; /* compute correlation functions */
  int nstcorr; /* time interval of computing autocorrelation function */
  double corrtol; /* error tolerance of the autocorrelation function */
  char fncorr[FILENAME_MAX]; /* file name for the autocorrelation function */

  long nequil; /* equilibration time */
  long nsteps; /* number of steps */
  long ntrials; /* number of trials */
  int verbose; /* verbose level */
  const char *prog; /* name of the program */
} invtpar_t;



enum {
  SAMPMETHOD_METROGLOBAL = 0,
  SAMPMETHOD_METROLOCAL,
  SAMPMETHOD_HEATBATH,
  SAMPMETHOD_COUNT
};


#define MAX_OPT_ALIASES 8

const char *sampmethod_names[][MAX_OPT_ALIASES] = {
  {"global Metropolis", "global", "g"},
  {"local Metropolis", "local", "l"},
  {"heat-bath", "h"},
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
  m->winn = 1; /* single bin update */
  for ( i = 1; i < m->winn; i++ ) {
    m->win[i] = 0;
  }
  m->win[0] = 1;
  m->wingaus = 0;
  m->initrand = 0;
  m->kcutoff = 0;
  m->sampmethod = 0;
  m->tcorr = 0.0; /* perfect sampling */

  m->docorr = 0; /* don't compute correlation functions */
  m->nstcorr = 0; /* do correlation functions */
  m->corrtol = 1e-8;
  m->fncorr[0] = '\0';

  m->nequil = m->n * 10000L;
  m->nsteps = 100000000L;
  m->ntrials = 100;
  m->prog = "invt";
  m->verbose = 0;
}



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

  if ( m->t0 <= 0 ) {
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

  /* construct the Gaussian window */
  if ( m->wingaus > 0 ) {
    m->winn = NBMAX;

    if ( m->winn >= m->n / 2 ) {
      m->winn = m->n / 2;
    }

    /* truncate the Gaussian at 5 sigma */
    if ( m->winn >= m->wingaus * 5 ) {
      m->winn = (int) (m->wingaus * 5 + 0.5);
    }

    x = 1;
    for ( i = 1; i < m->winn; i++ ) {
      m->win[i] = exp(-0.5*i*i/m->wingaus/m->wingaus);
      x += m->win[i] * 2;
    }

    /* normalize the window function, such that
     * win[0] + 2 * (win[1] + ... + win[n - 1]) = 1 */
    for ( i = 0; i < m->winn; i++ ) {
      m->win[i] /= x;
    }

  } else {

    for ( x = 0, i = 1; i < m->winn; i++ ) {
      x += m->win[i];
    }
    m->win[0] = 1 - 2 * x;
  }

  /* initialize the target distribution */
  if ( m->p ) {
    free(m->p);
  }
  xnew(m->p, m->n);
  for ( i = 0; i < m->n; i++ ) {
    m->p[i] = 1.0 / m->n;
  }

  /* set the initial cutoff wave number */
  if ( m->kcutoff <= 0 ) {
    m->kcutoff = m->n;
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
  fprintf(stderr, "  --nb=:         explicitly set the update window parameters, separated by comma, like --nb=0.2,0.4\n");
  fprintf(stderr, "  --sig=:        set the Gaussian window width, default %g\n", m->wingaus);
  fprintf(stderr, "  --initrand=:   magnitude of the initial random error, default %g\n", m->initrand);
  fprintf(stderr, "  --kcutoff=:    cutoff of wave number of the initial random error, default %d\n", m->kcutoff);
  fprintf(stderr, "  --samp=:       set the sampling scheme, g=global Metropolis, l=local Metropolis, h=heat-bath, default %s\n", sampmethod_names[m->sampmethod][0]);
  fprintf(stderr, "  --corr:        compute correlation functions, default %d\n", m->docorr);
  fprintf(stderr, "  --nstcorr=:    set the number of steps of setting the correlation function, default %d\n", m->nstcorr);
  fprintf(stderr, "  --corrtol=:    set the tolerance level to truncate the autocorrelation function, default %g\n", m->corrtol);
  fprintf(stderr, "  --fncorr=:     set the file name for the correlation function, default %s\n", m->fncorr);
  fprintf(stderr, "  --try=:        set the number of trials, default %ld\n", m->ntrials);
  fprintf(stderr, "  --step=:       set the number of simulation steps, default %ld\n", m->nsteps);
  fprintf(stderr, "  --equil=:      set the number of equilibration steps, default %ld\n", m->nequil);
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
  return atoi(val);
}



/* get long integer */
static long invtpar_getlong(invtpar_t *m,
    const char *key, const char *val)
{
  if ( val == NULL ) {
    fprintf(stderr, "no value for %s\n", key);
    invtpar_help(m);
  }
  return atol(val);
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
  else if ( strstartswith(key, "fixa") )
  {
    m->fixa = 1;
  }
  else if ( strstartswith(key, "nb")
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
  else if ( strcmpfuzzy(key, "tcorr") == 0
         || strcmpfuzzy(key, "corr-time") == 0 )
  {
    m->tcorr = invtpar_getdouble(m, key, val);
  }
  else if ( strcmpfuzzy(key, "corr") == 0
         || strcmpfuzzy(key, "docorr") == 0 )
  {
    m->docorr = 1;
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
    m->verbose = val ? atoi(val) : 1;
  }
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
      fprintf(stderr, "Unknown options %s = %s in %s\n",
          key, val, fn);
      getchar();
    }
  }

  fclose(fp);

  invtpar_compute(m);

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

  fprintf(stderr, "%ld trials: n %d, alpha = %g/(t + %g), alpha0 %g, "
      "%s, %ld steps;\n",
      m->ntrials, m->n, m->c, m->t0, m->alpha0,
      sampmethod_names[m->sampmethod][0], m->nsteps);
  fprintf(stderr, "update window function (%d bins): ", m->winn);
  for ( i = 0; i < m->winn; i++ ) {
    fprintf(stderr, "%g ", m->win[i]);
  }
  fprintf(stderr, "\n");
}



#endif /* INVTPAR__ */


