#include "invt.h"

/* test window trimming function */

static void test2(void)
{
  int winn = 2;
  double win[2] = {-0.2, 0.6};

  printf("\n\nTest n = 2 case\n");
  trimwindow(2, &winn, win, 0);
  printf("win %g, %g\n", win[0], win[1]);
}


static void test3(void)
{
  int winn = 3;
  double win[3] = {0.2, 0.4, 0.0};

  printf("\n\nTest n = 3 case\n");
  trimwindow(3, &winn, win, 0);
  printf("win %g, %g, %g, (%d bins)\n", win[0], win[1], win[2], winn);
}


int main(void)
{
  test2();
  test3();
  return 0;
}
