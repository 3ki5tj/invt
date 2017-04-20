#include "util.h"


/* metadynamics for integer */
typedef struct {
  int xmin, xmax, xdel;
  int n;
  double *v; /* bias potential */
  double *h; /* histogram */
  double a; /* updating magnitude */
} metad_t;



static metad_t *metad_open(int xmin, int xmax, int xdel)
{
  metad_t *metad;
  int i;

  xnew(metad, 1);
  metad->xmin = xmin;
  metad->xmax = xmax;
  metad->xdel = xdel;
  metad->n = (xmax - xmin) / xdel + 1;
  metad->xmax = metad->xmin + metad->n * metad->xdel;
  xnew(metad->v, metad->n);
  xnew(metad->h, metad->n);
  for ( i = 0; i < metad->n; i++ ) {
    metad->v[i] = 0;
    metad->h[i] = 0;
  }
  metad->a = 1;
  return metad;
}



static void metad_close(metad_t *metad)
{
  free(metad->v);
  free(metad->h);
  free(metad);
}



/* return the bin index of value x */
__inline static int metad_getindex(metad_t *metad, int x)
{
  if ( x < metad->xmin || x >= metad->xmax )
    return -1;
  return (x - metad->xmin) / metad->xdel;
}



/* decide if a transition from xold to xnew is to be accepted */
static int metad_acc(metad_t *metad, int iold, int xnew,
    int *inew)
{
  double dv, r;

  *inew = metad_getindex(metad, xnew);
  /* reject out-of-boundary moves */
  if ( *inew < 0 ) return 0;

  if ( iold == *inew ) return 1;
  dv = metad->v[*inew] - metad->v[iold];
  if ( dv <= 0 ) return 1;
  r = rand01();
  return r < exp(-dv);
}



static void metad_updatev(metad_t *metad, int i)
{
  metad->v[i] += metad->a;
  metad->h[i] += 1;
}



static void metad_trimv(metad_t *metad)
{
  int i;
  double v0 = metad->v[0];

  for ( i = 0; i < metad->n; i++ )
    metad->v[i] -= v0;
}


static int metad_save(metad_t *metad, const char *fn)
{
  FILE *fp;
  int i;

  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  metad_trimv(metad);
  fprintf(fp, "# %d %d %d %d\n",
      metad->n, metad->xmin, metad->xmax, metad->xdel);
  for ( i = 0; i < metad->n; i++ ) {
    fprintf(fp, "%d %g %g\n",
        metad->xmin + i * metad->xdel, metad->v[i], metad->h[i]);
  }
  fclose(fp);
  return 0;
}
