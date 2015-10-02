#include "declarations.h"

extern "C" double dgetrf_(int &M, int &N, double *A, int &LDA, int *IPIV, int &INFO);
extern "C" double dgetri_(int &N, double *A, int &LDA, int *IPIV, double *WORK, int &LWORK, int &INFO);
//double* ctof(double **in, int rows, int cols);
//void ftoc(double *in, double **out, int rows, int cols);

//void Invert(int N, double **A, double **InvA, FILE *logfile)
void Invert(int N, double *A, double *InvA, FILE *logfile)
{
  int i;
  int LWORK, INFO;
  int *IPIV;
  double *Atmp, *WORK;
  
  Atmp = new double[N*N];
  IPIV = new int[N];
  WORK = new double[N*N];
  LWORK = N*N;
  
  //Atmp = ctof(A,N,N);
  for(i=0 ; i < N*N ; i++) Atmp[i] = A[i];
  
  dgetrf_(N,N,Atmp,N,IPIV,INFO);
      
  if(INFO < 0){
    fprintf(logfile,"Error in matrix inversion in dgetrf. info=%d! Ending session.\n",INFO);
    exit(0);
  }
  if(INFO > 0){
    fprintf(logfile,"Singular matrix inversion in dgetrf. info=%d! Ending session.\n",INFO);
    exit(0);
  }
  
  dgetri_(N,Atmp,N,IPIV,WORK,LWORK,INFO);
  
  if(INFO < 0){
    fprintf(logfile,"Error in matrix inversion in dgetri. info=%d! Ending session.\n",INFO);
    exit(0);
  }
  if(INFO > 0){
    fprintf(logfile,"Singular matrix inversion in dgetri. info=%d! Ending session.\n",INFO);
    exit(0);
  }
  
  //ftoc(Atmp,InvA,N,N);
  for(i=0 ; i < N*N ; i++) InvA[i] = Atmp[i];
  
  delete [] IPIV;
  delete [] Atmp;
  delete [] WORK;
}

/*
double* ctof(double **in, int rows, int cols)
{
  int i, j;
  double *out = new double[rows*cols];
  
  for(i=0 ; i < cols ; i++) for(j=0 ; j < rows ; j++) out[j+i*rows] = in[i][j];

  return(out);
}

void ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0 ; i < cols ; i++) for (j=0 ; j < rows ; j++) out[i][j] = in[j+i*rows];
}
*/
