#include "declarations.h"

void gauss_rand(double *rnd1,double *rnd2);

double randv(long *idum)
{
  double rnd1, rnd2;

  rnd1 = ran_nrc(idum);
  rnd2 = ran_nrc(idum);

  gauss_rand(&rnd1,&rnd2);
  return rnd1;
}

/****************************************************************************/

void gauss_rand(double *rnd1,double *rnd2)
{
  double tmp1 = *rnd1, tmp2 = *rnd2;

  *rnd1 = sqrt(-2.0 * log(tmp1)) * cos(TWOPI * tmp2);
  *rnd2 = sqrt(-2.0 * log(tmp1)) * sin(TWOPI * tmp2);
  return;
}
