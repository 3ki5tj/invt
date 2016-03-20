#include <stdio.h>
#include "erf.h"

int main(void)
{
  printf("%g %g\n", erfc(1.0), erf(1.0));
  return 0;
}
