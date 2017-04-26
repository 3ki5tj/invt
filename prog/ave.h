#ifndef AVE_H__
#define AVE_H__

typedef struct {
  long cnt;
  double sum, sqr, ave, var;
} ave_t;

__inline static void ave_clear(ave_t *a)
{
  a->cnt = 0;
  a->sum = 0;
  a->sqr = 0;
  a->ave = 0;
  a->var = 0;
}

__inline static void ave_add(ave_t *a, double x)
{
  double dx;
  a->sum += x;
  if ( a->cnt > 0 ) { /* update the variance accumulator */
    dx = x - a->ave;
    a->sqr += dx * dx * a->cnt / (a->cnt + 1);
  }
  a->cnt += 1;
  a->ave = a->sum / a->cnt;
  a->var = a->sqr / a->cnt;
}


#endif
