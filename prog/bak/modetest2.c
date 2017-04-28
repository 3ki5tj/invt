/* test the residue error due to mode truncation */
#include "../util.h"
#include "../cosmodes.h"

int main(int argc, char **argv)
{
  int n = 1000, kcutoff = 50, pbc = 0, i;
  double *u, *v, *v1, *costab;
  double x, y, err, err1, maxdev;
  FILE *fp;

  if ( argc > 1 ) n = atoi(argv[1]);
  if ( argc > 2 ) kcutoff = atoi(argv[2]);
  if ( argc > 3 ) pbc = atoi(argv[3]);

  xnew(u, n);
  xnew(v, n);
  xnew(v1, n);
  costab = mkcostab(n, pbc);

  /* random noise */
  for ( i = 0; i < n; i++ ) {
    x = i * M_PI / n;
    if ( pbc ) x *= 2;
    v[i] = x;
  }

  /* convert them to modes */
  getcosmodes(v, n, u, costab);

  /* eliminate high frequency modes */
  for ( err1 = 0, i = kcutoff + 1; i < (pbc ? n - kcutoff : n); i++ ) {
    err1 += u[i] * u[i];
    u[i] = 0;
  }

  /* convert back */
  fromcosmodes(v1, n, u, costab);

  fp = fopen("modecmp.dat", "w");
  /* compute the error */
  maxdev = 0;
  for ( err = 0, i = 0; i < n; i++ ) {
    y = fabs(v1[i] - v[i]);
    x = i*M_PI/n;
    if ( pbc ) x *= 2;
    fprintf(fp, "%g %g %g %g\n", x, v[i], v1[i], u[i]);
    err += y * y;
    if ( y > maxdev ) maxdev = y;
  }
  err /= n;
  fclose(fp);

  printf("err %g, %g, sqrt %g, maxdev %g\n", err, err1, sqrt(err), maxdev);

  free(u);
  free(v);
  free(v1);
  free(costab);

  return 0;
}
