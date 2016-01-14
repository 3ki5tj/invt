/* test the eigenmode decomposition code */
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


static int output(const double *v1, const double *v2, int n,
    const char *fn)
{
  int i;
  FILE *fp;


  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }

  for ( i = 0; i < n; i++ ) {
    fprintf(fp, "%d\t%g\t%g\n", i, v1[i], v2[i]);
  }

  fclose(fp);

  return 0;
}


int main(void)
{
  int n = 100, i, j;
  double *u, *u1, *v, *v1, *v2, *costab;

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

  /* convert them to modes */
  getcosmodes(v, n, u, costab);

  /* convert back */
  fromcosmodes(v1, n, u, costab);

  /* truncate high frequency modes */
  for ( i = 0; i < n; i++ ) {
    u1[i] = ( i <= 10 ) ? u[i] : 0;
  }

  /* convert back */
  fromcosmodes(v2, n, u1, costab);

  /* print out u and v */
  if ( n <= 12 ) {
    printarr(v, n, "  ");
    printarr(u, n, "  ");
    printarr(v1, n, "  ");
    printarr(v2, n, "  ");
  }

  output(v1, v2, n, "out.dat");

  free(u);
  free(u1);
  free(v);
  free(v1);
  free(v2);
  free(costab);
}
