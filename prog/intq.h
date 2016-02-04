#ifndef INTQ_H__
#define INTQ_H__



/* integrating the exact optimal alpha(t) */



typedef struct {
  int m; /* number of integration point */
  double t; /* number of times */
  double qt;
  double *tarr;
  double *qarr;
  double *aarr;
  int n;
  const double *lamarr; /* eigenvalues of the updating magnitude */
  const double *gamma; /* autocorrelation integrals */
  double Ea, Er, E; /* square root errors */
} intq_t;



static intq_t *intq_open(double t, double qt, int m,
    int n, const double *lamarr, const double *gamma)
{
  intq_t *intq;

  xnew(intq, 1);
  intq->m = m;
  xnew(intq->tarr, m + 1);
  xnew(intq->aarr, m + 1);
  xnew(intq->qarr, m + 1);
  intq->t = t;
  intq->qt = qt;
  intq->n = n;
  intq->lamarr = lamarr;
  intq->gamma = gamma;
  intq->Ea = intq->Er = intq->E = 0;
  return intq;
}



static void intq_close(intq_t *intq)
{
  free(intq->tarr);
  free(intq->aarr);
  free(intq->qarr);
  free(intq);
}



/* compute the integrand for `intq` */
static double intq_getvel(intq_t *intq, double q)
{
  int i, n = intq->n;
  double lam, y = 0, dq = q - intq->qt;

  for ( i = 1; i < n; i++ ) {
    lam = intq->lamarr[i];
    y += intq->gamma[i] * lam * lam * exp( 2 * lam * dq );
  }

  return sqrt( y );
}



/* compute the optimal protocol
 * for a period t, and fixed
 * qt = Integral {from 0 to t} alpha(t') dt'
 * the results are saved to tarr[0..m], qarr[0..m], and aarr[0..m]
 * return the square-root asymptotic error
 * */
static void intq_geta(intq_t *intq)
{
  double c, dq, y;
  int j, m = intq->m;

  /* integrate over q over uniform grid */
  dq = intq->qt / m;
  intq->tarr[0] = 0;
  intq->qarr[0] = 0;
  y = intq_getvel(intq, intq->qarr[0]);
  for ( j = 1; j <= m; j++ ) {
    intq->qarr[j] = j * dq;
    intq->tarr[j] = intq->tarr[j - 1] + y * 0.5 * dq;
    y = intq_getvel(intq, intq->qarr[j]);
    intq->tarr[j] += y * 0.5 * dq;
  }

  /* normalize the time array */
  c = intq->t / intq->tarr[m];
  for ( j = 1; j <= m; j++ ) {
    intq->tarr[j] *= c;
  }

  /* differentiate q(t) */
  y = ( intq->qarr[1] - intq->qarr[0] )
    / ( intq->tarr[1] - intq->tarr[0] );
  intq->aarr[0] = y;
  for ( j = 1; j < m; j++ ) {
    intq->aarr[j] += 0.5 * y;
    y = ( intq->qarr[j + 1] - intq->qarr[j] )
      / ( intq->tarr[j + 1] - intq->tarr[j] );
    intq->aarr[j] += 0.5 * y;
  }
  intq->aarr[m] = y;
}



/* compute the square-root asymptotic error
 * assuming intq_geta() has been called
 * */
static double intq_asymerr(intq_t *intq)
{
  double err, y, dt;
  int j, m = intq->m;

  /* compute the error */
  err = 0;
  for ( j = 0; j <= m; j++ ) {
    y = intq_getvel(intq, intq->qarr[j]) * intq->aarr[j];
    y = y * y;
    if ( j == 0 ) {
      dt = intq->tarr[1] - intq->tarr[0];
    } else if ( j == m ) {
      dt = intq->tarr[m] - intq->tarr[m - 1];
    } else {
      dt = intq->tarr[j+1] - intq->tarr[j-1];
    }
    err += y * 0.5 * dt;
  }

  intq->Ea = err;

  return sqrt( err );
}



/* compute the square-root residue error */
static double intq_reserr(intq_t *intq, double a0)
{
  int i;
  double err = 0;

  for ( i = 1; i < intq->n; i++ ) {
    err += 0.5 * a0 * intq->gamma[i] * intq->lamarr[i]
         * exp(-2.0 * intq->lamarr[i] * intq->qt);
  }
  intq->Er = err;
  return sqrt( err );
}



/* from the error from the optimal alpha(t) */
static double intq_geterr(intq_t *intq, double a0)
{
  /* compute the optimal schedule alpha(t) */
  intq_geta(intq);

  /* compute the asymptotic error */
  intq_asymerr(intq);

  /* compute the residual error */
  intq_reserr(intq, a0);

  intq->E = intq->Ea + intq->Er;

  return sqrt( intq->E );
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
static int intq_save(intq_t *intq,
    double c, double t0, const char *fn)
{
  FILE *fp;
  int i, m = intq->m;
  double a1, q1;

  if ( fn == NULL || fn[0] == '\0' ) {
    return -1;
  }

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }

  fprintf(fp, "# %d %d %g %g %g\n",
      m, intq->n, intq->E, intq->Ea, intq->Er);
  for ( i = 0; i < m; i++ ) {
    a1 = c / (intq->tarr[i] + t0);
    q1 = c * log( 1 + intq->tarr[i] / t0 );
    fprintf(fp, "%g\t%20.8e\t%12.6f\t%20.8e\t%12.6f\n",
        intq->tarr[i], intq->aarr[i], intq->qarr[i], a1, q1);
  }

  fclose(fp);

  fprintf(stderr, "saved optimal schedule to %s\n",
     fn);

  return 0;
}



/* return the square-root error from the optimal schedule  */
__inline static double opterror_ez(double c, double t, double a0,
    int m, const char *fnarr,
    int n, int winn, double *win, int sampmethod,
    intq_t **intq_ptr, int verbose)
{
  double t0, qt, err;
  double *lamarr, *gamma;
  intq_t *intq;

  /* estimate the eigenvalues of the w matrix,
   * for the updating scheme */
  lamarr = esteigvals(n, winn, win);

  /* estimate the autocorrelation integrals
   * of the eigenmodes of the w matrix,
   * for the updating scheme */
  gamma = estgamma(n, sampmethod);

  qt = c * log(1 + a0 * t / c);
  t0 = c / a0;

  intq = intq_open(t, qt, m, n, lamarr, gamma);

  /* compute the optimal schedule and error */
  err = intq_geterr(intq, a0);

  /* save the schedule to file */
  intq_save(intq, c, t0, fnarr);

  if ( verbose >= 1 ) {
    fprintf(stderr, "estimated final error %g, sqr: %e\n",
        sqrt(err), err);
  }

  if ( intq_ptr != NULL ) {
    *intq_ptr = intq;
  } else {
    free(lamarr);
    free(gamma);
    intq_close(intq);
  }
  return err;
}




#endif /* INTQ_H__ */
