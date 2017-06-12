#ifndef INTQ_H__
#define INTQ_H__



/* integrating the exact optimal alpha(t) */



typedef struct {
  int m; /* number of integration point */
  double T; /* total simulation time */
  double *tarr; /* time grid */
  double *qarr; /* q(t) */
  double *aarr; /* alpha(t) */
  double *dinva; /* instantaneous eigenvalue d(1/alpha(t))/dt */
  int n; /* number of bins or modes */
  int K; /* cutoff mode */
  int pbc; /* periodic boundary condition */
  const double *lambda; /* eigenvalues of the updating magnitude */
  const double *gamma; /* autocorrelation integrals */
  double Ea, Er, E; /* square errors */
  double EKa, EKr, EK; /* mode-limited errors */
  int curr_id; /* bin index for interpolation */
} intq_t;



static intq_t *intq_open(double T, int m,
    int n, int K, int pbc,
    const double *lambda, const double *gamma)
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
  intq->K = K;
  intq->pbc = pbc;
  intq->lambda = lambda;
  intq->gamma = gamma;
  intq->Ea = intq->Er = intq->E = 0;
  intq->EKa = intq->EKr = intq->EK = 0;
  intq->curr_id = 0;
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


/* set the q grid points */
static void intq_setqgrid(intq_t *intq, double qT)
{
  int i, m = intq->m;
  double dlnq = log(qT+1)/m;
  for ( i = 0; i <= m; i++ ) {
    intq->qarr[i] = (qT+1) * (1 - exp(-i*dlnq));
  }
}



/* compute the integrand for `intq` */
static double intq_getmassx(intq_t *intq, double dq,
    double *massK)
{
  int k, n = intq->n, K = intq->K, pbc = intq->pbc;
  double lam, y, mass = 0;

  *massK = 0;
  for ( k = 1; k < n; k++ ) {
    lam = intq->lambda[k];
    y = intq->gamma[k] * lam * lam * exp( 2 * lam * dq );
    mass += y;
    if ( K < 0 || k <= K || (pbc && k >= n - K) ) {
      *massK += y; /* for mode limited optimal schedule */
    }
  }

  if ( mass < 0 || mass > 1e10 ) {
    fprintf(stderr, "mass %g, dq %g\n", mass, dq);
  }
  mass = sqrt( mass );
  *massK = sqrt( *massK );
  return mass;
}



/* compute the optimal schedule, in terms of q(t)
 * the results are saved to tarr[0..m], qarr[0..m]
 * */
static double intq_getq(intq_t *intq, double qT)
{
  int j, m = intq->m;
  double c, dq, yK = 0, dt;

  intq_setqgrid(intq, qT);
  intq->tarr[0] = 0;
  intq->qarr[0] = 0;
  intq_getmassx(intq, intq->qarr[0] - qT, &yK);
  for ( j = 1; j <= m; j++ ) {
    dq = intq->qarr[j] - intq->qarr[j-1];
    dt = yK * 0.5 * dq;
    intq_getmassx(intq, intq->qarr[j] - qT, &yK);
    dt += yK * 0.5 * dq;
    intq->tarr[j] = intq->tarr[j-1] + dt;
  }

  /* normalize the time array */
  c = intq->T / intq->tarr[m];
  for ( j = 1; j <= m; j++ ) {
    intq->tarr[j] *= c;
  }
  return c;
}



/* differentiate q(t) to get alpha(t) */
static void intq_diffq(intq_t *intq)
{
  double y, lam = intq->lambda[1];
  int j, j1, j0, m = intq->m;

  /* differentiate q(t) */
  for ( j = 0; j <= m; j++ ) {
    j1 = (j < m ? j + 1 : m);
    j0 = (j > 0 ? j - 1 : 0);
    y = ( intq->qarr[j1] - intq->qarr[j0] )
      / ( intq->tarr[j1] - intq->tarr[j0] );
    intq->aarr[j] = (1 - exp(-y*lam)) / lam;
  }
}



/* differentiate 1/a(t) to get instantaneous eigenvalue */
static void intq_diffinva(intq_t *intq)
{
  double y;
  int j, j0, j1, m = intq->m;

  for ( j = 0; j <= m; j++ ) {
    j1 = (j < m ? j + 1 : m);
    j0 = (j > 0 ? j - 1 : 0);
    y = ( 1.0 / intq->aarr[j1] - 1.0 / intq->aarr[j0] )
      / (       intq->tarr[j1] -       intq->tarr[j0] );
    intq->dinva[j] = y;
  }
}



/* compute the optimal schedule, alpha(t)
 * for a period t, and fixed
 * qT = Integral {from 0 to t} alpha(t') dt'
 * the results are saved to tarr[0..m], qarr[0..m], and aarr[0..m]
 * */
static void intq_geta(intq_t *intq, double qT)
{
  intq_getq(intq, qT);

  /* differentiate q(t) to get alpha(t) */
  intq_diffq(intq);
}



/* retrieve the initial updating magnitude */
__inline static double intq_getinita(intq_t *intq)
{
  return intq->aarr[0];
}


/* compute the asymptotic error
 * assuming intq_geta() has been called
 * */
static double intq_asymerr(intq_t *intq, double qT)
{
  double y, yK, dt;
  int j, m = intq->m;

  /* compute the error */
  intq->Ea = intq->EKa = 0;
  for ( j = 0; j <= m; j++ ) {
    y = intq_getmassx(intq, intq->qarr[j] - qT, &yK);
    y *= intq->aarr[j];
    yK *= intq->aarr[j];
    if ( j == 0 ) {
      dt = 0.5 * (intq->tarr[1] - intq->tarr[0]);
    } else if ( j == m ) {
      dt = 0.5 * (intq->tarr[m] - intq->tarr[m - 1]);
    } else {
      dt = 0.5 * (intq->tarr[j+1] - intq->tarr[j-1]);
    }
    intq->Ea += y * y * dt;
    intq->EKa += yK * yK * dt;
  }

  return sqrt( intq->Ea );
}



/* compute the square-root residue error */
static double intq_reserr(intq_t *intq, double a0, const double *xerr, double qT)
{
  int i, n = intq->n, K = intq->K, pbc = intq->pbc;
  double y;

  intq->Er = intq->EKr = 0;
  for ( i = 1; i < n; i++ ) {
    y = 0.5 * a0 * intq->gamma[i] * intq->lambda[i];
    if ( xerr != NULL ) y += xerr[i];
    y *= exp(-2.0 * intq->lambda[i] * qT);
    intq->Er += y;
    if ( K < 0 || i <= K || (!pbc && i >= n - K) )
      intq->EKr += y;
  }
  return sqrt( intq->Er );
}



/* compute the error components */
__inline static double intq_errcomp(intq_t *intq, double a0,
    double qT, double *xerr, double *xerr_r, double *xerr_a)
{
  int i, j, m = intq->m, n = intq->n;
  // K = intq->K, pbc = intq->pbc;
  double err_r, err_a, err;
  double a, y, dq, dt, errtot = 0, lam, gam;

  /* loop over modes, starting from mode 1
   * since mode 0 is always zero */
  for ( i = 1; i < n; i++ ) {
    lam = intq->lambda[i];
    gam = intq->gamma[i];

    /* residual error */
    err_r = 0.5 * a0 * gam * lam * exp(-2.0 * lam * qT);

    /* asymptotic error */
    err_a = 0;
    for ( j = 0; j <= m; j++ ) {
      dq = intq->qarr[j] - qT;
      a = intq->aarr[j];
      y = gam * lam * lam * exp( 2 * lam * dq ) * a * a;
      if ( j == 0 ) {
        dt = 0.5 * (intq->tarr[1] - intq->tarr[0]);
      } else if ( j == m ) {
        dt = 0.5 * (intq->tarr[m] - intq->tarr[m - 1]);
      } else {
        dt = 0.5 * (intq->tarr[j+1] - intq->tarr[j-1]);
      }
      err_a += y * dt;
    }

    if ( xerr_r != NULL ) xerr_r[i] = err_r;
    if ( xerr_a != NULL ) xerr_a[i] = err_a;

    err = err_r + err_a;
    xerr[i] = err;
    errtot += err;
  }

  return sqrt( errtot );
}



/* compute the error from the optimal alpha(t) */
static double intq_geterrx(intq_t *intq, double a0,
    const double *xerr, double qT)
{
  /* compute the optimal schedule alpha(t) */
  intq_geta(intq, qT);

  /* compute the asymptotic error */
  intq_asymerr(intq, qT);

  /* compute the residual error */
  intq_reserr(intq, a0, xerr, qT);

  intq->E = intq->Ea + intq->Er;
  intq->EK = intq->EKa + intq->EKr;

  return sqrt( intq->E );
}



/* evaulate
 *   mint = Integrate M(Q) dQ,
 * where
 *   M(Q) = sqrt{ Sum_k Gamma_k lambda_k^2 e^{-2 lambda_k Q} }.
 * Ea should be equal to mint^2 / T.
 * */
__inline static double intq_getmint(intq_t *intq, double qT,
    double *mint_ubound, const char *fnmass)
{
  int i, k, m = intq->m, n = intq->n;
  double q, dq, lambda, gamma, xp, y, mass, mint;
  double mass_ub, mint_ub;
  FILE *fpmass = NULL;

  intq_setqgrid(intq, qT);
  mint = 0;
  mint_ub = 0;

  if ( fnmass != NULL &&
      (fpmass = fopen(fnmass, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fnmass);
  }

  /* compute Int M(Q) dQ */
  for ( i = 0; i <= m; i++ ) {
    q = qT - intq->qarr[m-i];
    dq = (i < m ? intq->qarr[m-i] - intq->qarr[m-i-1] : 0)
       + (i > 0 ? intq->qarr[m-i+1] - intq->qarr[m-i] : 0);
    y = 0;
    mass_ub = 0;
    for ( k = 1; k < n; k++ ) {
      gamma = intq->gamma[k];
      lambda = intq->lambda[k];
      xp = exp(-2 * lambda * q);
      y += gamma * lambda * lambda * xp;
      mass_ub += lambda * sqrt(gamma * xp);
    }
    mass = sqrt( y );
    mint += mass * dq * 0.5;
    mint_ub += mass_ub * dq * 0.5;
    if ( fpmass != NULL ) {
      fprintf(fpmass, "%g\t%g\t%g\t%g\t%g\n", q, mass, mint, mass_ub, mint_ub);
    }
  }

  if ( mint_ubound != NULL ) {
    *mint_ubound = mint_ub;
  }
  fprintf(stderr, "qT %g, mint %g <= %g\n", qT, mint, mint_ub);

  if ( fpmass != NULL ) {
    fclose(fpmass);
  }
  return mint;
}



/* evaulation the function
 *   F(qT) = Int {0 to qT} M(Q) dQ - M(qT) T initalpha
 * where
 *   M^2(Q) =  Sum_k Gamma_k lambda_k^2 e^{-2 lambda_k Q}.
 * and the derivative
 *  F'(qT) = M(qT) + (T a0/M(qT)) Sum_k Gamma_k lambda_k^3 e^{-2 lambda_k Q}
 * */
static double intq_optqTfunc(intq_t *intq,
    double initalpha, double qT, double *df)
{
  int i, k, m = intq->m;
  int n = intq->n, K = intq->K, pbc = intq->pbc;
  double f, dq, q, lambda, gamma, xp, y, mint, mass;

  intq_setqgrid(intq, qT);
  mint = 0;

  /* compute Int M(Q) dQ */
  for ( i = 0; i <= m; i++ ) {
    q = qT - intq->qarr[m-i];
    dq = (i < m ? intq->qarr[m-i] - intq->qarr[m-i-1] : 0)
       + (i > 0 ? intq->qarr[m-i+1] - intq->qarr[m-i] : 0);
    y = 0;
    for ( k = 1; k < n; k++ ) {
      if ( K >= 0 && k > K && (!pbc || k < n - K) )
        continue;
      gamma = intq->gamma[k];
      lambda = intq->lambda[k];
      xp = exp(-2 * lambda * q);
      y += gamma * lambda * lambda * xp;
    }
    mass = sqrt( y );
    mint += mass * dq * 0.5;
  }
  f = mint - mass * intq->T * initalpha;

  for ( y = 0, k = 1; k < n; k++ ) {
    if ( K >= 0 && k > K && (!pbc || k < n - K) )
      continue;
    gamma = intq->gamma[k];
    lambda = intq->lambda[k];
    xp = exp(-2 * lambda * q);
    y += gamma * lambda * lambda * lambda * xp;
    //printf("k %d, gamma %g, lambda %g, xp %g\n", k, gamma, lambda, xp);
  }

  *df = mass + y / mass * intq->T * initalpha;
  //printf("qT %g, f %g, mass %g, df %g, y %g, K %d\n", qT, f, mass, *df, y, K);
  return f;
}



/* compute the optimal q(T) such that the initial
 * updating magnitude is initalpha
 * Solving the equation
 *   Int {0 to q(T)} M(Q) dQ - M(q(T)) T initalpha == 0
 * by the Newton-Raphson method
 * */
static double intq_optqT(intq_t *intq, double initalpha,
    double tol, int verbose)
{
  double qT = log(1 + intq->T * initalpha / 2), dlnq;
  double f = DBL_MAX, df;
  double qTmin = 1e-10, qTmax = 1e10, ql = 1e-10, qr = 1e10;
  int i;

  for ( i = 0; fabs(f) > 1e-10; i++ ) {
    f = intq_optqTfunc(intq, initalpha, qT, &df);

    /* update the bracket */
    if ( f < 0 && qT > ql ) ql = qT;
    if ( f > 0 && qT < qr ) qr = qT;

    if ( df < fabs(f/qT) ) df = fabs(f/qT);
    dlnq = -f/df/qT;

    if ( qT > qTmax && f < 0 || qT < qTmin && f > 0 || fabs(dlnq) < tol )
      break; /* break if qT is too large and the error is still decreasing */

    /* limit qT within the bracket */
    if ( qT > qTmax ) {
      qT = qTmax - 0.1 * (qTmax - qTmin);
    } else if ( qT < qTmin ) {
      qT = qTmin + 0.1 * (qTmax - qTmin);
    }

    if ( verbose ) {
      fprintf(stderr, "%d: qT %g -> %g (%g, %g), dlnq %g, f %g, df %g\n",
          i, qT, qT*exp(dlnq), qTmin, qTmax, dlnq, f, df);
    }
    qT *= exp(dlnq);
    if ( qT < ql || qT > qr ) qT = sqrt(ql * qr);
  }
  return qT;
}



/* evaulation the function
 *   F(qT) = M(qT) Int {0 to qT} M(Q) dQ - M^2(qT) a0 T / 2
 *         - T Sum_k err_k lambda_k e^{-2 lambda_k Q}
 * where
 *   M^2(Q) =  Sum_k Gamma_k lambda_k^2 e^{-2 lambda_k Q}.
 * and the derivative
 *  F'(qT) = M^2(qT) + M'(qT) IntM - a0 T M(qT) M'(qT))
 *         + 2 T Sum_k xerr_k lambda_k^2 e^{-2 lambda_k Q}.
 * where M'(qT) = Sum_k Gamma_k lambda_k^3 e^{-2 lambda_k Q} / M^2(Q)
 * */
static double intq_optqTfuncx(intq_t *intq, double a0,
    double *xerr, double qT, double *df)
{
  int i, k, m = intq->m;
  int n = intq->n, K = intq->K, pbc = intq->pbc;
  double f, q, dq, lambda, gamma, xp, y, z, w, mint, mass, dmass, T;

  intq_setqgrid(intq, qT);
  mint = 0;

  /* compute Int M(Q) dQ */
  for ( i = 0; i <= m; i++ ) {
    q = qT - intq->qarr[m-i];
    dq = (i < m ? intq->qarr[m-i] - intq->qarr[m-i-1] : 0)
       + (i > 0 ? intq->qarr[m-i+1] - intq->qarr[m-i] : 0);
    y = 0;
    for ( k = 1; k < n; k++ ) {
      if ( K >= 0 && k > K && (!pbc || k < n - K) )
        continue;
      gamma = intq->gamma[k];
      lambda = intq->lambda[k];
      xp = exp(-2 * lambda * q);
      y += gamma * lambda * lambda * xp;
    }
    mass = sqrt( y );
    mint += mass * dq * 0.5;
  }

  y = z = w = 0;
  for ( k = 1; k < n; k++ ) {
    if ( K >= 0 && k > K && (!pbc || k < n - K) )
      continue;
    gamma = intq->gamma[k];
    lambda = intq->lambda[k];
    xp = exp(-2 * lambda * q);
    y += gamma * lambda * lambda * lambda * xp;
    if ( k > 0 ) {
      z += xerr[k] * lambda * xp;
      w += xerr[k] * lambda * lambda * xp;
    }
  }

  dmass = -y/mass;
  T = intq->T;
  f = mass*mint - mass*mass*a0*T/2 - z*T;
  *df = mass*mass + dmass*mint - a0*T*mass*dmass + 2*w*T;
  return f;
}




/* compute the optimal q(T) such that the initial
 * updating magnitude is initalpha
 * Solving the equation (derivative of the total error)
 *   M(q(T))/T Int {0 to q(T)} M(Q) dQ - M^2(q(T)) a0/2
 *    - Sum_k xerr_k lambda_k e^{-2 lambda_k qT} == 0
 * by the Newton-Raphson method
 * */
static double intq_optqTx(intq_t *intq, double a0,
    double *xerr, double tol, int verbose)
{
  double qT = log(1 + intq->T * a0), dlnq;
  double f = DBL_MAX, df;
  double qTmin = 1e-10, qTmax = 1e10, ql = 1e-10, qr = 1e10;
  int i;

  for ( i = 0; fabs(f) > 1e-10; i++ ) {
    f = intq_optqTfuncx(intq, a0, xerr, qT, &df);
    //double f1, df1; f1 = intq_optqTfuncx(intq, a0, xerr, qT+0.01, &df1);

    /* update the bracket */
    if ( f < 0 && qT > ql ) ql = qT;
    if ( f > 0 && qT < qr ) qr = qT;

    if ( df < fabs(f/qT) ) df = fabs(f/qT);
    dlnq = -f/df/qT;
    //printf("f %g, df %g,%g q %g, dlnq %g\n", f, df, (f1-f)/0.01, qT, dlnq); getchar();

    if ( qT > qTmax*0.1 && f < 0 || qT < qTmin*10 && f > 0 || fabs(dlnq) < tol )
      break; /* break if qT is too large and the error is still decreasing */

    if ( verbose ) {
      fprintf(stderr, "%d: qT %g -> %g (%g, %g), dlnq %g, f %g, df %g\n",
          i, qT, qT*exp(dlnq), qTmin, qTmax, dlnq, f, df);
    }
    qT *= exp(dlnq);
    if ( qT < ql || qT > qr ) qT = sqrt(ql * qr);
    //printf("qT %g, qTmin %g\n", qT, qTmin); getchar();
  }
  return qT;
}




#if 0 /* backup routine, replaced by intq_optqT() */
/* variate the value of qT to compute
 * the optimal schedule alpha(t) and the error
 * this function is adapted from estbestc_invt()
 * in `invt.h` */
static double intq_minerr(intq_t *intq, double a0,
    double *qT, double prec, int verbose)
{
  double ql, qm, qr, qn;
  double el, em, er, en;
  int it;

  /* specify the initial bracket */
  qn = *qT;
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

  *qT = qm;

  return em;
}
#endif



/* low level interpolation code */
__inline static double intq_interp(double t, int *id,
    int m, const double *tarr, const double *yarr)
{
  int i;
  double r, a;

  /* find the proper bin */
  i = *id;
  if ( i < 0 ) i = 0;
  if ( i >= m ) i = m - 1;
  /* seek backward */
  if ( i < m ) {
    for ( ; i > 0; i-- ) {
      if ( t >= tarr[i] )
        break;
    }
  }
  /* seek forward */
  if ( i >= 0 ) {
    for ( ; i < m; i++ ) {
      if ( t < tarr[i + 1] )
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
  r = ( t - tarr[i] ) / ( tarr[i + 1] - tarr[i] );

  a = (1 - r) * yarr[i] + r * yarr[i + 1];
  if ( a < 0 ) a = 0;
  return a;
}



/* compute the updating magnitude alpha at a given t
 * newer discretization */
__inline static double intq_evala(intq_t *intq, double t)
{
  double t0 = (t >= 1 ? t - 1 : 0), q0, q1, lam = intq->lambda[1];
  q0 = intq_interp(t0, &intq->curr_id, intq->m, intq->tarr, intq->qarr);
  q1 = intq_interp(t,  &intq->curr_id, intq->m, intq->tarr, intq->qarr);
  return (1 - exp(lam*(q0 - q1)))/lam;
}



/* compute the updating magnitude alpha at a given t
 * by interpolation
 * assuming intq_geta() has been called */
__inline static double intq_interpa(intq_t *intq, double t)
{
  return intq_interp(t, &intq->curr_id, intq->m, intq->tarr, intq->aarr);
}



/* save optimal protocol to file */
__inline static int intq_save(intq_t *intq,
    double c, double t0, int resample, const char *fn)
{
  FILE *fp;
  int i, id, m = intq->m;
  double t, a, q, dinva, T, qT;

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

  fprintf(fp, "# %d %d %g %g %g %g %.10f %g %g\n",
      m, intq->n, intq->E, intq->Ea, intq->Er,
      T, qT, t0, c);
  for ( i = 0; i <= m; i++ ) {
    if ( resample ) {
      t = T * i / m;
      id = i;
      q = intq_interp(t, &id, m, intq->tarr, intq->qarr);
      a = intq_interp(t, &id, m, intq->tarr, intq->aarr);
      dinva = intq_interp(t, &id, m, intq->tarr, intq->dinva);
    } else {
      t = intq->tarr[i];
      a = intq->aarr[i];
      q = intq->qarr[i];
      dinva = intq->dinva[i];
    }
    fprintf(fp, "%12.3f\t%20.8e\t%12.6f\t%20.8e\n",
        t, a, q, dinva);
  }

  fclose(fp);

  fprintf(stderr, "saved optimal schedule to %s\n", fn);

  return 0;
}



/* return the square-root error from the optimal schedule */
static double esterror_opt(double T, double a0,
    double initalpha, double *qT, double qprec,
    int m, intq_t **intq_ptr,
    int n, int K, int pbc,
    const double *lambda, const double *gamma,
    int verbose)
{
  intq_t *intq;
  double err;

  intq = intq_open(T, m, n, K, pbc, lambda, gamma);

  /* compute the optimal schedule and error */
  if ( *qT <= 0 ) {
    /* by default, the initial updating magnitude
     * takes the optimal value, which is
     * half of the equilibration value */
    if ( initalpha <= 0 ) initalpha = a0 / 2;
    *qT = intq_optqT(intq, initalpha, qprec, verbose);
    //printf("qT %g\n", *qT); getchar();
  }
  err = intq_geterrx(intq, a0, NULL, *qT);

  if ( intq_ptr != NULL ) {
    *intq_ptr = intq;
  } else {
    /* free the intq object if it no longer needed */
    intq_close(intq);
  }

  return err;
}



/* return the square-root error from the optimal schedule */
static double esterror_optx(double T, double a0,
    double *xerr, double *qT, double qprec,
    int m, intq_t **intq_ptr,
    int n, int K, int pbc,
    const double *lambda, const double *gamma,
    int verbose)
{
  intq_t *intq;
  double err;

  intq = intq_open(T, m, n, K, pbc, lambda, gamma);

  /* compute the optimal schedule and error */
  if ( *qT <= 0 ) {
    *qT = intq_optqTx(intq, a0, xerr, qprec, verbose);
    //printf("qT %g\n", *qT); getchar();
  }
  err = intq_geterrx(intq, a0, xerr, *qT);

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
