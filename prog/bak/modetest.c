/* test the eigenmode decomposition code
 * that is the pair of functions
 * getcosmodes() and fromcosmodes()
 * */
#include "invt.h"



static void printarr(const double *a, int n, const char *delim)
{
  int i;

  for ( i = 0; i < n; i++ ) {
    printf("%10.4f", a[i]);
    if ( i < n - 1 ) printf(delim);
  }
  printf("\n");
}



static int output(const double *v1, const double *v2,
    const double *u, int n, const char *fn)
{
  int i;
  FILE *fp;


  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "%d\t%g\t%g\t%g\n", i, v1[i], v2[i], u[i]);
  }

  fclose(fp);

  return 0;
}



int main(int argc, char **argv)
{
  int n = 2, kcutoff = 10, i;
  double *u, *u1, *v, *v1, *v2, *costab;
  double errv, erru, errv2;

  if ( argc > 1 ) {
    n = atoi(argv[1]);
  }

  if ( argc > 2 ) {
    kcutoff = atoi(argv[2]);
  }

  xnew(u, n);
  xnew(u1, n);
  xnew(v, n);
  xnew(v1, n);
  xnew(v2, n);
  costab = makecostab(n);

  /* random noise */
  for ( i = 0; i < n; i++ ) {
    v[i] = randgaus();
  }

  /* normalize the error */
  normalize(v, n, 1.0, NULL);

  /* print out the field for a small n */
  if ( n <= 10 ) {
    printf("   i:      vi        vi^2\n");
    for ( i = 0; i < n; i++ ) {
      printf("%4d: %10.6f %10.6f\n", i, v[i], v[i] * v[i]);
    }
    printf("\n");
  }

  /* errv^2 = sum_i v_i^2 / n */
  errv = geterror(v, n, NULL);

  /* convert them to modes */
  getcosmodes(v, n, u, costab);

  /* print out the field for a small n */
  if ( n <= 10 ) {
    printf("   i:      ui        ui^2\n");
    for ( i = 0; i < n; i++ ) {
      printf("%4d: %10.6f %10.6f\n", i, u[i], u[i] * u[i]);
    }
    printf("\n");
  }

  /* erru^2 = sum_i u_i^2 */
  erru = getuerror(u, n);

  /* convert back */
  fromcosmodes(v1, n, u, costab);

  /* truncate high frequency modes */
  for ( i = 0; i < n; i++ ) {
    u1[i] = ( i < kcutoff ) ? u[i] : 0;
  }

  /* convert back */
  fromcosmodes(v2, n, u1, costab);

  errv2 = geterror(v2, n, NULL);

  printf("n %d, errv %g, erru %g, errv2 %g\n",
      n, errv, erru, errv2);

  /* print out u and v */
  if ( n <= 12 ) {
    printarr(v, n, "  ");
    printarr(u, n, "  ");
    printarr(v1, n, "  ");
    printarr(v2, n, "  ");
  }

  output(v1, v2, u, n, "out.dat");

  free(u);
  free(u1);
  free(v);
  free(v1);
  free(v2);
  free(costab);

  return 0;
}
