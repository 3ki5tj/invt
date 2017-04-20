#include "invt.h"
#define IS2_LB 4
#include "is2.h"
#include "metad.h"

const double IS2_TC = 2.3;

static int invt_is2_run(invtpar_t *m)
{
  is2_t *is;
  int id, h;
  long t;
  int enew, iold, inew, acc;
  metad_t *metad;

  mtscramble( clock() );

  is = is2_open(IS2_L);

  /* equilibration at the critical temperature */
  is2_setuproba(1.0/IS2_TC, is->uproba);
  for ( t = 1; t <= m->nequil; t++ ) {
    IS2_PICK(is, id, h);
    if ( h <= 0 || mtrand() <= is->uproba[h] ) {
      IS2_FLIP(is, id, h);
    }
  }

  /* production run */
  metad = metad_open(-IS2_N*2+8, 0, 4);
  iold = metad_getindex(metad, is->E);
  for ( t = 1; t <= m->nsteps; t++ ) {
    IS2_PICK(is, id, h);
    enew = is->E + h * 2;
    acc = metad_acc(metad, iold, enew, &inew);
    if ( acc ) {
      iold = inew;
      IS2_FLIP(is, id, h);
    }
    metad_updatev(metad, iold);
  }

  metad_save(metad, "vbias.dat");
  is2_close(is);
  metad_close(metad);

  return 0;
}

int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  m->nequil = 10000;
  m->nsteps = 1000000;
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  invt_is2_run(m);
  invtpar_finish(m);
  return 0;
}
