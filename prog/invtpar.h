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
  double alpha0; /* initial updating magnitude */
  long nequil; /* equilibration time */
  long nsteps; /* number of steps */
  int nbn; /* number of neighboring bins */
  double nbs[NBMAX + 1]; /* magnitudes of neighbors */
  int sampmethod; /* sampling method */
  double tcorr; /* correlation time */
  int ntrials; /* number of trials */
  int verbose; /* verbose level */
  const char *prog; /* name of the program */
} invtpar_t;



enum {
  SAMPMETHOD_MCLOCAL = 0,
  SAMPMETHOD_MCGLOBAL,
  SAMPMETHOD_HEATBATH,
  SAMPMETHOD_COUNT
};

const char *sampmethod_names[][8] = {
  {"MC-local", "local", "l"},
  {"MC-global", "global", "g"},
  {"heat-bath", "h"},
  {""}
};



static void invtpar_init(invtpar_t *p)
{
  int i;

  p->c = 1.0;
  p->t0 = 0;
  p->n = 10;
  p->alpha0 = 0.0;
  p->nequil = p->n * 10000L;
  p->nsteps = 100000000L;
  p->tcorr = 0.0; /* perfect sampling */
  p->ntrials = 100;
  p->nbn = 1; /* single bin update */
  for ( i = 1; i < p->nbn; i++ ) {
    p->nbs[i] = 0;
  }
  p->nbs[0] = 1;
  p->prog = "invt";
  p->verbose = 0;
}



static void invtpar_compute(invtpar_t *m)
{
  if ( m->alpha0 <= 0 ) {
    m->alpha0 = 1.0 / m->n;
  }
  if ( m->t0 <= 0 ) {
    m->t0 = m->c / m->alpha0;
  }
}



static int invtpar_load(invtpar_t *m, const char *fn)
{
  FILE *fp;
  char buf[800], *p, *key, *val;
  int inpar;

  if ( (fp = fopen(fn, "r")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  while ( fgets(buf, sizeof buf, fp) ) {
    strstrip(buf); /* remove trailing spaces */
    if ( buf[0] == '\0' || buf[0] == '#' ) continue;

    /* parse the line to a key and a value */
    /* find the end of the key */
    inpar = 0; /* within (...) */
    for ( p = buf; *p; p++ ) {
      if ( (!inpar && isspace(*p)) || *p == '=' ) {
        *p = '\0'; /* end the key part */
        break;
      }
      if ( !inpar && (*p == '(' || *p == '[') )
        inpar = 1; /* enter a parentheses block */
      else if ( inpar && (*p == ')' || *p == ']') )
        inpar = 0; /* leave a parentheses block */
      *p = (char) tolower(*p);
    }
    key = buf;

    /* find the beginning of the value */
    for ( p++; isspace(*p) || *p == '=' ; )
      p++;
    val = p;
    for ( ; *p; p++ )
      *p = (char) tolower(*p);

    if ( strcmp(key, "n") == 0 ) {
      m->n = atoi(val);
    } else if ( strcmp(key, "c") == 0 ) {
      m->c = atof(val);
    } else if ( strcmp(key, "t0") == 0 ) {
      m->t0 = atof(val);
    } else if ( strcmpfuzzy(key, "equil") == 0
             || strcmpfuzzy(key, "nequil") == 0 ) {
      m->nequil = atoi(val);
    } else if ( strcmpfuzzy(key, "steps") == 0
             || strcmpfuzzy(key, "nsteps") == 0 ) {
      m->nsteps = atoi(val);
    } else if ( strcmpfuzzy(key, "steps") == 0
             || strcmpfuzzy(key, "nsteps") == 0 ) {
      m->nsteps = atoi(val);
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
  fprintf(stderr, "  --try=:        set the number of trials, default %d\n", m->ntrials);
  fprintf(stderr, "  --step=:       set the number of simulation steps, default %ld\n", m->nsteps);
  fprintf(stderr, "  --equil=:      set the number of equilibration steps, default %ld\n", m->nequil);
  fprintf(stderr, "  -v:            be verbose, -vv to be more verbose, etc.\n");
  fprintf(stderr, "  -h, --help:    display this message\n");
  exit(1);
}



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
      } else if ( strcmp(p, "c") == 0 ) {
        m->c = atof(q);
      } else if ( strcmp(p, "t0") == 0 ) {
        m->t0 = atof(q);
      } else if ( strcmp(p, "nb") == 0
          || strcmp(p, "neighbor") == 0 ) {
        double sum = 0;
        m->nbn = 1 + readarray(m->nbs + 1, NBMAX, q);
        for ( j = 1; j < m->nbn; j++ ) {
          sum += m->nbs[j];
        }
        m->nbs[0] = 1 - 2 * sum;
      } else if ( strcmp(p, "equil") == 0
               || strcmp(p, "nequil") == 0 ) {
        m->nequil = atoi(q);
      } else if ( strstartswith(p, "ntrial")
          || strstartswith(p, "trials")
          || strcmp(p, "ntry") == 0
          || strcmp(p, "try") == 0 ) {
        m->ntrials = atoi(q);
      } else if ( strstartswith(p, "steps")
               || strstartswith(p, "nsteps") ) {
        m->nsteps = atoi(q);
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



static void invtpar_dump(const invtpar_t *m)
{
  int i;

  fprintf(stderr, "ntrials %d, n %d, c %g\n",
      m->ntrials, m->n, m->c);
  fprintf(stderr, "update matrix (%d): ", m->nbn);
  for ( i = 0; i < m->nbn; i++ ) {
    fprintf(stderr, "%g ", m->nbs[i]);
  }
  fprintf(stderr, "\n");
}

#endif /* INVTPAR__ */


