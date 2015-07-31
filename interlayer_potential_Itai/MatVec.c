#include "declarations.h"

extern "C" void dgemv_(char &TRANSA, int &M, int &N, double &ALPHA,  double *A,  int &LDA, double *X, int &INCX, double &BETA, double *Y, int &INCY);
//double* ctof(double **in, int rows, int cols);

//void MatVec(double **A, double *VecIn, int N, double *VecOut)
void MatVec(double *A, double *VecIn, int N, double *VecOut)
{
  char TRANSA;
  int INCX, INCY;
  double ALPHA, BETA;
  //double *ATMP;
  
  TRANSA = 'N';
  ALPHA  = 1.0;
  BETA   = 0.0;
  INCX   = 1;
  INCY   = 1;
  
  //ATMP = ctof(A,N,N);
  
  //dgemv_(TRANSA,N,N,ALPHA,ATMP,N,VecIn,INCX,BETA,VecOut,INCY);
  dgemv_(TRANSA,N,N,ALPHA,A,N,VecIn,INCX,BETA,VecOut,INCY);

  //delete [] ATMP;
}
