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
  int nbn; /* width of the updating window function */
  double nbs[NBMAX + 1]; /* shape of the window function */
  double wgaus; /* width of the Gaussian window */
  double initrand; /* magnitude of the initial Gaussian noise */
  int sampmethod; /* sampling method */
  double tcorr; /* correlation time */
  long nequil; /* equilibration time */
  long nsteps; /* number of steps */
  int ntrials; /* number of trials */
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
  {"Metropolis-global", "global", "g"},
  {"Metropolis-local", "local", "l"},
  {"heat-bath", "h"},
  {""}
};



/* parse a string into an array */
static int readarray(double *arr, int n, char *s)
{
  const char *delims = ",;:";
  char *p;
  int i = 0;

  p = strtok(s, delims);
  while ( p != NULL ) {
    arr[ i++ ] = atof(p);
    p = strtok(NULL, delims);
    if ( i >= n )
      break;
  }
  return i;
}



/* select an option */
static int selectoption(const char *s,
    const char *names[][MAX_OPT_ALIASES], int cnt)
{
  int i, j;

  /* loop over options */
  for ( i = 0; i < cnt; i++ ) {
    /* loop over aliases */
    for ( j = 0; j < MAX_OPT_ALIASES; j++ ) {
      if ( names[i][j] == NULL
        || names[i][j][0] == '\0' ) {
        break;
      }

      if ( strcmpfuzzy(names[i][j], s) == 0 ) {
        return i;
      }
    }
  }

  /* try to treat `s` is a number */
  if ( isdigit(s[0]) ) {
    i = atoi(s);
    if ( i < 0 || i >= cnt ) i = 0;
  } else {
    i = 0;
  }

  return i;
}



/* initialize the default parameters */
static void invtpar_init(invtpar_t *m)
{
  int i;

  m->c = 1.0;
  m->t0 = 0;
  m->n = 10;
  m->p = NULL;
  m->alpha0 = 0.0;
  m->nequil = m->n * 10000L;
  m->nsteps = 100000000L;
  m->tcorr = 0.0; /* perfect sampling */
  m->ntrials = 100;
  m->nbn = 1; /* single bin update */
  for ( i = 1; i < m->nbn; i++ ) {
    m->nbs[i] = 0;
  }
  m->nbs[0] = 1;
  m->wgaus = 0;
  m->initrand = 0;
  m->sampmethod = 0;
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

  /* construct the Gaussian window */
  if ( m->wgaus > 0 ) {
    m->nbn = NBMAX;

    if ( m->nbn >= m->n / 2 ) {
      m->nbn = m->n / 2;
    }

    /* truncate the Gaussian at 5 sigma */
    if ( m->nbn >= m->wgaus * 5 ) {
      m->nbn = (int) (m->wgaus * 5 + 0.5);
    }

    x = 1;
    for ( i = 1; i < m->nbn; i++ ) {
      m->nbs[i] = exp(-0.5*i*i/m->wgaus/m->wgaus);
      x += m->nbs[i] * 2;
    }

    /* normalize */
    for ( i = 0; i < m->nbn; i++ ) {
      m->nbs[i] /= x;
    }

  } else {

    for ( x = 0, i = 1; i < m->nbn; i++ ) {
      x += m->nbs[i];
    }
    m->nbs[0] = 1 - 2 * x;
  }

  /* initialize the target distribution */
  if ( m->p ) {
    free(m->p);
  }
  xnew(m->p, m->n);

  for ( i = 0; i < m->n; i++ ) {
    m->p[i] = 1.0 / m->n;
  }
}



/* clean up */
static void invtpar_finish(invtpar_t *m)
{
  free(m->p);
  m->p = NULL;
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
    strstrip(buf); /* remove trailing spaces */
    if ( buf[0] == '\0' || buf[0] == '#' ) continue;

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

    if ( strcmp(key, "n") == 0 ) {
      m->n = atoi(val);
    } else if ( strcmpfuzzy(key, "c") == 0
             || strcmpfuzzy(key, "invt-c") ) {
      m->c = atof(val);
    } else if ( strcmpfuzzy(key, "t0") == 0 ) {
      m->t0 = atof(val);
    } else if ( strcmpfuzzy(key, "alpha0") == 0
             || strcmpfuzzy(key, "a0") == 0 ) {
      m->alpha0 = atof(val);
    } else if ( strstartswith(key, "nb")
             || strstartswith(key, "neighbo")
             || strcmpfuzzy(key, "window") == 0 ) {
      m->nbn = 1 + readarray(m->nbs + 1, NBMAX, val);
    } else if ( strcmpfuzzy(key, "wgaus") == 0
             || strcmpfuzzy(key, "width-Gaussian") == 0
             || strcmpfuzzy(key, "sigma") == 0
             || strcmpfuzzy(key, "sig") == 0 ) {
      m->wgaus = atof(val);
    } else if ( strcmpfuzzy(key, "initrand") == 0 ) {
      m->initrand = atof(val);
    } else if ( strcmpfuzzy(key, "sampling-method") == 0
             || strcmpfuzzy(key, "sampmethod") == 0 ) {
      m->sampmethod = selectoption(val,
          sampmethod_names, SAMPMETHOD_COUNT);
    } else if ( strcmpfuzzy(key, "tcorr") == 0
             || strcmpfuzzy(key, "corr-time") == 0 ) {
      m->tcorr = atof(val);
    } else if ( strcmpfuzzy(key, "equil") == 0
             || strcmpfuzzy(key, "nequil") == 0 ) {
      m->nequil = atoi(val);
    } else if ( strcmpfuzzy(key, "steps") == 0
             || strcmpfuzzy(key, "nsteps") == 0 ) {
      m->nsteps = atoi(val);
    } else if ( strcmpfuzzy(key, "try") == 0
             || strcmpfuzzy(key, "repeat") == 0
             || strstartswith(key, "trial")
             || strstartswith(key, "ntrial") ) {
      m->ntrials = atoi(val);
    } else if ( strcmpfuzzy(key, "verbose") == 0 ) {
      m->verbose = val ? atoi(val) : 1;
    } else {
      fprintf(stderr, "Unknown options %s = %s in %s\n",
          key, val, fn);
      getchar();
    }
  }

  fclose(fp);

  invtpar_compute(m);

  return 0;
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
  fprintf(stderr, "  --t0=:         set t0 in alpha = c/(t + t0), default %g\n", m->t0);
  fprintf(stderr, "  --nb=:         explicitly set the update window parameters, separated by comma, like --nb=0.2,0.4\n");
  fprintf(stderr, "  --sig=:        set the Gaussian window width, default %g\n", m->wgaus);
  fprintf(stderr, "  --initrand=:   randomize the error initially, default %g\n", m->initrand);
  fprintf(stderr, "  --try=:        set the number of trials, default %d\n", m->ntrials);
  fprintf(stderr, "  --step=:       set the number of simulation steps, default %ld\n", m->nsteps);
  fprintf(stderr, "  --equil=:      set the number of equilibration steps, default %ld\n", m->nequil);
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc.\n");
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
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

      if ( strcmp(p, "n") == 0 ) {
        m->n = atoi(q);
      } else if ( strcmpfuzzy(p, "c") == 0
               || strcmpfuzzy(p, "invtc") == 0 ) {
        m->c = atof(q);
      } else if ( strcmpfuzzy(p, "alpha0") == 0
               || strcmpfuzzy(p, "a0") == 0 ) {
        m->alpha0 = atof(q);
      } else if ( strcmp(p, "t0") == 0 ) {
        m->t0 = atof(q);
      } else if ( strstartswith(p, "nb")
               || strstartswith(p, "neighbor")
               || strcmpfuzzy(p, "window")  == 0 ) {
        m->nbn = 1 + readarray(m->nbs + 1, NBMAX, q);
      } else if ( strcmpfuzzy(p, "wgaus") == 0
               || strcmpfuzzy(p, "width-Gaussian") == 0
               || strcmpfuzzy(p, "sigma") == 0
               || strcmpfuzzy(p, "sig") == 0 ) {
        m->wgaus = atof(q);
      } else if ( strcmpfuzzy(p, "initrand") == 0 ) {
        m->initrand = atof(q);
      } else if ( strstartswith(p, "samp") ) {
        m->sampmethod = selectoption(q,
            sampmethod_names, SAMPMETHOD_COUNT);
      } else if ( strcmpfuzzy(p, "tcorr") == 0
               || strcmpfuzzy(p, "corr-time") == 0 ) {
        m->tcorr = atof(q);
      } else if ( strcmp(p, "equil") == 0
               || strcmp(p, "nequil") == 0 ) {
        m->nequil = atoi(q);
      } else if ( strstartswith(p, "steps")
               || strstartswith(p, "nsteps") ) {
        m->nsteps = atoi(q);
      } else if ( strcmp(p, "try") == 0
               || strcmp(p, "repeat") == 0
               || strstartswith(p, "trial")
               || strstartswith(p, "ntrial") ) {
        m->ntrials = atoi(q);
      } else if ( strcmp(p, "verbose") == 0 ) {
        m->verbose = atoi(q);
      } else if ( strcmp(p, "help") == 0 ) {
        invtpar_help(m);
      } else {
        fprintf(stderr, "Unknown option %s\n", argv[i]);
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

  fprintf(stderr, "ntrials %d, n %d, alpha = %g/t, alpha0 %g, "
      "%s, tcorr %g\n",
      m->ntrials, m->n, m->c, m->alpha0,
      sampmethod_names[m->sampmethod][0], m->tcorr);
  fprintf(stderr, "update window function (%d bins): ", m->nbn);
  for ( i = 0; i < m->nbn; i++ ) {
    fprintf(stderr, "%g ", m->nbs[i]);
  }
  fprintf(stderr, "\n");
}



#endif /* INVTPAR__ */


