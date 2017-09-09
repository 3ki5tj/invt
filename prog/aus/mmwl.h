#ifndef MMWL_H__
#define MMWL_H__



/* Wang-Landau scheme for moments */



#include <math.h>



typedef struct {
  double mm[3]; /* statistical moments */
  double fl[3]; /* fluctuation */
  double alphawl; /* Wang-Landau updating magnitude */
  int invt; /* using the 1/t schedule */
  double t0;
} mmwl_t;



__inline static void mmwl_init(mmwl_t *mm, double a)
{
  mm->mm[0] = 0;
  mm->mm[1] = 0;
  mm->mm[2] = 0;
  mm->fl[0] = 0;
  mm->fl[1] = 0;
  mm->fl[2] = 0;
  mm->alphawl = a;
  mm->invt = 0;
  mm->t0 = 1;
}



__inline static void mmwl_add(mmwl_t *mm, double y1, double y2)
{
  mm->mm[0] += 1;
  mm->mm[1] += y1;
  mm->mm[2] += y2;
}



__inline static double mmwl_getalpha(mmwl_t *mm)
{
  return mm->invt ? 1./(mm->mm[0] + mm->t0) : mm->alphawl;
}



/* compute the fluctuation of the statistical moments
 * return the maximum fluctuation */
__inline static double mmwl_calcfl(mmwl_t *mm, int dovar)
{
  double t = mm->mm[0], fl, fl2;
  if ( t <= 0 ) return 99.0;
  mm->fl[1] = mm->mm[1] / t;
  mm->fl[2] = mm->mm[2] / t;
  fl = fabs(mm->fl[1]);
  if ( dovar && (fl2 = fabs(mm->fl[2])) > fl )
    fl = fl2;
  return fl;
}



/* switch to a stage with a reduced updating magnitude */
__inline static void mmwl_switch(mmwl_t *mm, double magred)
{
  double t = mm->mm[0];
  mm->alphawl *= magred;
  if ( mm->alphawl < 1./t ) {
    mm->t0 = t;
    mm->invt = 1;
  }
  mm->mm[0] = 0;
  mm->mm[1] = 0;
  mm->mm[2] = 0;
}


__inline static int mmwl_check(mmwl_t *mm, int dovar, double fl, double magred)
{
  double flmax = mmwl_calcfl(mm, dovar);
  if ( !mm->invt && flmax < fl ) {
    mmwl_switch(mm, magred);
    return 1;
  } else {
    return 0;
  }
}



#endif /* MMWL_H__ */
