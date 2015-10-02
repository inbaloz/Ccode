#include "declarations.h"

extern "C" void dsysv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* WORK, int* LWORK, int* INFO);
double* DSYSV_ctof(double **in, int rows, int cols);
void DSYSV_ftoc(double *in, double **out, int rows, int cols);

void DSYSV(int N, int NRHS, double **A, double **B, double **Vout)
{
  char UPLO;
  int LDA, LDB, LWORK, INFO;
  int *IPIV;
  double *WORK, WorkOpt, *Atmp, *Btmp;
  
  UPLO = 'U';
  LDA = LDB = ((N > 1) ? N : 1);
  IPIV = new int[N];
  INFO = -1;
  
  Atmp = DSYSV_ctof(A,N,LDA);
  Btmp = DSYSV_ctof(B,N,NRHS);
  
  // Determine the optimal LWORK
  
  LWORK = -1;
  dsysv_(&UPLO,&N,&NRHS,Atmp,&LDA,IPIV,Btmp,&LDB,&WorkOpt,&LWORK,&INFO);
  
  if(INFO == 0) LWORK = int(WorkOpt);
  else{
    cerr<<"Error: could not determine the optimal LWORK in DSYSV! Ending session.\n";
    exit(0);
  }
  
  // Perform calculation
  
  WORK = new double[LWORK];
  
  dsysv_(&UPLO,&N,&NRHS,Atmp,&LDA,IPIV,Btmp,&LDB,WORK,&LWORK,&INFO);
  
  if(INFO != 0) cerr<<"Error in DSYSV! INFO="<<INFO<<". Ending session.\n";
  
  DSYSV_ftoc(Btmp,Vout,N,NRHS);
  
  delete [] IPIV;
  delete [] WORK;
  delete [] Atmp;
  delete [] Btmp;
}

double* DSYSV_ctof(double **in, int rows, int cols)
{
  int i, j;
  double *out = new double[rows*cols];
  
  for(i=0 ; i < cols ; i++) for(j=0 ; j < rows ; j++) out[j+i*rows] = in[i][j];

  return(out);
}

void DSYSV_ftoc(double *in, double **out, int rows, int cols)
{
  int i, j;

  for (i=0 ; i < cols ; i++) for (j=0 ; j < rows ; j++) out[i][j] = in[j+i*rows];
}
