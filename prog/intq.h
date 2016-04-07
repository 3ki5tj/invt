#ifndef INTQ_H__
#define INTQ_H__



/* integrating the exact optimal alpha(t) */



typedef struct {
  int m; /* number of integration point */
  double T; /* total simulation time */
  double *tarr; /* time grid */
  double *qarr;
  double *aarr;
  double *dinva; /* instantaneous eigenvalue */
  int n;
  const double *lambda; /* eigenvalues of the updating magnitude */
  const double *gamma; /* autocorrelation integrals */
  double Ea, Er, E; /* square root errors */
} intq_t;



static intq_t *intq_open(double T, int m,
    int n, const double *lambda, const double *gamma)
{
  intq_t *intq;

  xnew(intq, 1);
  intq->m = m;
  xnew(intq->tarr, m + 1);
  xnew(intq->aarr, m + 1);
  xnew(intq->qarr, m + 1);
  xnew(intq->dinva, m + 1);
  intq->T = T;
  intq->n = n;
  intq->lambda = lambda;
  intq->gamma = gamma;
  intq->Ea = intq->Er = intq->E = 0;
  return intq;
}



static void intq_close(intq_t *intq)
{
  if ( intq != NULL ) {
    free(intq->tarr);
    free(intq->aarr);
    free(intq->qarr);
    free(intq->dinva);
    free(intq);
  }
}



/* compute the integrand for `intq` */
static double intq_getvel(intq_t *intq, double dq)
{
  int i, n = intq->n;
  double lam, y = 0;

  for ( i = 1; i < n; i++ ) {
    lam = intq->lambda[i];
    y += intq->gamma[i] * lam * lam * exp( 2 * lam * dq );
  }

  return sqrt( y );
}



/* compute the optimal schedule, in terms of q(t)
 * the results are saved to tarr[0..m], qarr[0..m]
 * */
static void intq_getq(intq_t *intq, double qt)
{
  double c, dq, y;
  int j, m = intq->m;

  /* integrate over q over uniform grid */
  dq = qt / m;
  intq->tarr[0] = 0;
  intq->qarr[0] = 0;
  y = intq_getvel(intq, intq->qarr[0] - qt);
  for ( j = 1; j <= m; j++ ) {
    intq->qarr[j] = j * dq;
    intq->tarr[j] = intq->tarr[j - 1] + y * 0.5 * dq;
    y = intq_getvel(intq, intq->qarr[j] - qt);
    intq->tarr[j] += y * 0.5 * dq;
  }

  /* normalize the time array */
  c = intq->T / intq->tarr[m];
  for ( j = 1; j <= m; j++ ) {
    intq->tarr[j] *= c;
  }
}



/* differentiate q(t) to get alpha(t) */
static void intq_diffq(intq_t *intq)
{
  double y;
  int j, m = intq->m;

  /* differentiate q(t) */
  y = ( intq->qarr[1] - intq->qarr[0] )
    / ( intq->tarr[1] - intq->tarr[0] );
  intq->aarr[0] = y;
  for ( j = 1; j < m; j++ ) {
    y = ( intq->qarr[j + 1] - intq->qarr[j - 1] )
      / ( intq->tarr[j + 1] - intq->tarr[j - 1] );
    intq->aarr[j] = y;
  }
  y = ( intq->qarr[m] - intq->qarr[m - 1] )
    / ( intq->tarr[m] - intq->tarr[m - 1] );
  intq->aarr[m] = y;
}



/* differentiate 1/a(t) to get instantaneous eigenvalue */
static void intq_diffinva(intq_t *intq)
{
  double y;
  int j, m = intq->m;

  y = ( 1.0 / intq->aarr[1] - 1.0 / intq->aarr[0] )
    / (       intq->tarr[1] -       intq->tarr[0] );
  intq->dinva[0] = y;
  for ( j = 1; j < m; j++ ) {
    y = ( 1.0 / intq->aarr[j + 1] - 1.0 / intq->aarr[j - 1] )
      / (       intq->tarr[j + 1] -       intq->tarr[j - 1] );
    intq->dinva[j] = y;
  }
  y = ( 1.0 / intq->aarr[m] - 1.0 / intq->aarr[m - 1] )
    / (       intq->tarr[m] -       intq->tarr[m - 1] );
  intq->dinva[m] = y;
}



/* compute the optimal schedule, alpha(t)
 * for a period t, and fixed
 * qt = Integral {from 0 to t} alpha(t') dt'
 * the results are saved to tarr[0..m], qarr[0..m], and aarr[0..m]
 * */
static void intq_geta(intq_t *intq, double qt)
{
  intq_getq(intq, qt);

  /* differentiate q(t) to get alpha(t) */
  intq_diffq(intq);
}



/* compute the square-root asymptotic error
 * assuming intq_geta() has been called
 * */
static double intq_asymerr(intq_t *intq, double qt)
{
  double err, y, dt;
  int j, m = intq->m;

  /* compute the error */
  err = 0;
  for ( j = 0; j <= m; j++ ) {
    y = intq_getvel(intq, intq->qarr[j] - qt) * intq->aarr[j];
    if ( j == 0 ) {
      dt = intq->tarr[1] - intq->tarr[0];
    } else if ( j == m ) {
      dt = intq->tarr[m] - intq->tarr[m - 1];
    } else {
      dt = intq->tarr[j+1] - intq->tarr[j-1];
    }
    err += y * y * 0.5 * dt;
  }

  intq->Ea = err;

  return sqrt( err );
}



/* compute the square-root residue error */
static double intq_reserr(intq_t *intq, double a0, double qt)
{
  int i;
  double err = 0;

  for ( i = 1; i < intq->n; i++ ) {
    err += 0.5 * a0 * intq->gamma[i] * intq->lambda[i]
         * exp(-2.0 * intq->lambda[i] * qt);
  }
  intq->Er = err;
  return sqrt( err );
}



/* compute the error components */
static double intq_errcomp(intq_t *intq, double a0,
    double qt, double *xerr)
{
  int i, j, m = intq->m;
  double y, dq, dt, err, errtot = 0, lam, gam;

  /* loop over modes, starting from mode 1
   * since mode 0 is always zero */
  for ( i = 1; i < intq->n; i++ ) {
    lam = intq->lambda[i];
    gam = intq->gamma[i];

    /* residual error */
    err = 0.5 * a0 * gam * lam * exp(-2.0 * lam * qt);

    /* asymptotic error */
    for ( j = 0; j <= m; j++ ) {
      dq = intq->qarr[j] - qt;
      y = gam * lam * lam * exp( 2 * lam * dq ) * intq->aarr[j];
      if ( j == 0 ) {
        dt = intq->tarr[1] - intq->tarr[0];
      } else if ( j == m ) {
        dt = intq->tarr[m] - intq->tarr[m - 1];
      } else {
        dt = intq->tarr[j+1] - intq->tarr[j-1];
      }
      err += y * y * 0.5 * dt;
    }

    xerr[i] = err;

    errtot += err;
  }

  return sqrt( errtot );
}



/* compute the error from the optimal alpha(t) */
static double intq_geterr(intq_t *intq, double a0,
    double qt)
{
  /* compute the optimal schedule alpha(t) */
  intq_geta(intq, qt);

  /* compute the asymptotic error */
  intq_asymerr(intq, qt);

  /* compute the residual error */
  intq_reserr(intq, a0, qt);

  intq->E = intq->Ea + intq->Er;
  //fprintf(stderr, "a0 %g, qt %g, Er %g Ea %g, E %g\n", a0, qt, intq->Er, intq->Ea, intq->E); getchar();

  return sqrt( intq->E );
}



/* variate the value of qt to compute
 * the optimal schedule alpha(t) and the error
 * this function is adapted from estbestc_invt()
 * in `invt.h` */
static double intq_minerr(intq_t *intq, double a0,
    double *qt, double prec, int verbose)
{
  double ql, qm, qr, qn;
  double el, em, er, en;
  int it;

  /* specify the initial bracket */
  qn = *qt;
  if ( qn <= 0 ) {
    qn = log(1 + intq->T * a0);
  }
  ql = 0.5 * qn;
  qm = 1.0 * qn;
  qr = 2.0 * qn;

  /* compute the errors at the initial backet */
  el = intq_geterr(intq, a0, ql);
  //intq_save(intq, 1, 1, "ql.dat");
  em = intq_geterr(intq, a0, qm);
  //intq_save(intq, 1, 1, "qm.dat");
  er = intq_geterr(intq, a0, qr);
  //intq_save(intq, 1, 1, "qr.dat");

  for ( it = 1; ; it++ ) {
    if ( verbose >= 1 ) {
      fprintf(stderr, "%d: %g (%g) - %g (%g) - %g (%g)\n",
          it, ql, el, qm, em, qr, er);
      if ( verbose >= 3 ) {
        fprintf(stderr, "Press enter to continue...\n");
        getchar();
      }
    }

    /* update the bracket */
    if ( el < em && el < er ) {
      /* el is the least of the three, extend to the left
       * <-- L - M - R */
      qr = qm;
      qm = ql;
      ql = ql * 0.5;
      er = em;
      em = el;
      el = intq_geterr(intq, a0, ql);
    } else if ( er < em && er < el ) {
      /* er is the least of the three, extend to the right
       * L - M - R --> */
      ql = qm;
      qm = qr;
      qr = qr * 2.0;
      el = em;
      em = er;
      er = intq_geterr(intq, a0, qr);
    } else {
      /* break the loop */
      if ( qr - ql < prec ) {
        break;
      }

      /* em is the least of the three */
      if ( qm - ql > qr - qm ) {
        /* refine the left side */
        qn = (ql + qm) * 0.5;
        en = intq_geterr(intq, a0, qn);
        if ( en > em ) {
          /* L - (N - M - R) */
          ql = qn;
          el = en;
        } else {
          /* (L - N - M) - R */
          qr = qm;
          qm = qn;
          er = em;
          em = en;
        }
      } else {
        /* refine the right side */
        qn = (qm + qr) * 0.5;
        en = intq_geterr(intq, a0, qn);
        if ( en > em ) {
          /* (L - M - N) - R */
          qr = qn;
          er = en;
        } else {
          /* L - (M - N - R) */
          ql = qm;
          qm = qn;
          el = em;
          em = en;
        }
      }
    }
  }

  *qt = qm;

  return em;
}



/* compute the updating magnitude alpha at a given t
 * by interpolation
 * assuming intq_geta() has been called */
__inline static double intq_interpa(intq_t *intq, double t, int *id)
{
  int i, m = intq->m;
  double r, a;

  /* find the proper bin */
  i = *id;
  /* seek backward */
  if ( i < m ) {
    for ( ; i > 0; i-- ) {
      if ( t >= intq->tarr[i] )
        break;
    }
  }
  /* seek forward */
  if ( i >= 0 ) {
    for ( ; i < m; i++ ) {
      if ( t < intq->tarr[i + 1] )
        break;
    }
  }
  //if ( i != *id ) { printf("i %d -> %d\n", *id, i); getchar(); }
  *id = i;

  /* compute alpha(t) by interpolation or extrapolation */
  if ( i < 0 ) { /* extrapolate */
    i = 0;
  } else if ( i >= m ) { /* extrapolate */
    i = m - 1;
  }

  /* interpolation formula */
  r = ( t                 - intq->tarr[i] )
    / ( intq->tarr[i + 1] - intq->tarr[i] );

  a = (1 - r) * intq->aarr[i] + r * intq->aarr[i + 1];
  if ( a < 0 ) a = 0;
  return a;
}



/* save optimal protocol to file */
__inline static int intq_save(intq_t *intq,
    double c, double t0, const char *fn)
{
  FILE *fp;
  int i, m = intq->m;
  double a1, q1, T, qT;

  if ( fn == NULL || fn[0] == '\0' ) {
    return -1;
  }

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  T = intq->T;
  qT = intq->qarr[m];
  if ( t0 <= 0 )
    t0 = T / (exp(qT/c) - 1);

  /* differentiate 1/a(t) */
  intq_diffinva(intq);

  fprintf(fp, "# %d %d %g %g %g %g %g %g\n",
      m, intq->n, intq->E, intq->Ea, intq->Er,
      T, qT, t0);
  for ( i = 0; i < m; i++ ) {
    a1 = c / (intq->tarr[i] + t0);
    q1 = c * log( 1 + intq->tarr[i] / t0 );
    fprintf(fp, "%12.3f\t%20.8e\t%12.6f\t%20.8e\t"
        "%20.8e\t%12.6f\t%20.8e\n",
        intq->tarr[i], intq->aarr[i], intq->qarr[i], intq->dinva[i],
        a1, q1, 1.0 / c);
  }

  fclose(fp);

  fprintf(stderr, "saved optimal schedule to %s\n",
     fn);

  return 0;
}



/* return the square-root error from the optimal schedule  */
static double esterror_opt(double T, double a0, double *qt,
    double qprec, int m, intq_t **intq_ptr,
    int n, double *xerr, const double *lambda, const double *gamma,
    int verbose)
{
  intq_t *intq;
  double err;

  intq = intq_open(T, m, n, lambda, gamma);

  /* compute the optimal schedule and error */
  err = intq_minerr(intq, a0, qt, qprec, verbose);

  /* compute the error components */
  if ( xerr != NULL ) {
    intq_errcomp(intq, a0, *qt, xerr);
  }

  if ( intq_ptr != NULL ) {
    *intq_ptr = intq;
  } else {
    /* free the intq object if it no longer needed */
    intq_close(intq);
  }

  return err;
}



/* estimate t0 for computing the normalized error
 * (T + t0) * E */
__inline static double intq_estt0(double T, double qT)
{
  return T / (exp(qT) - 1);
}



#endif /* INTQ_H__ */
