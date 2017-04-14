#include "invt.h"
#include "is2.h"

static int invt_is2_run(invtpar_t *m)
{
  mtscramble( clock() );

  return 0;
}

int main(int argc, char **argv)
{
  invtpar_t m[1];

  invtpar_init(m);
  invtpar_doargs(m, argc, argv);
  invtpar_dump(m);
  invt_is2_run(m);
  invtpar_finish(m);
  return 0;
}
