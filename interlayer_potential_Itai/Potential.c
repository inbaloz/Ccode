

#include "declarations.h"


//*****************************************************************************************************************************

//******************************************************************************************************************************
double calc_Potential(int Natoms, int *NBond, int *Bonded, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, atom *particle, double *Bmat, double *qvec, double *Avec, double *rijShiled_1, FILE* logfile, double *rij,double *angle, int *BondNum,Normal_struct *Normal,double *Fc,int *Neighb,int *NNeighb,point L,int Coulomb_flag)
{
  
  int atomi;
  double ETersoff,Epot;
  
  ETersoff = E_Tersoff(Natoms,Bonded,particle,rij,Fc,BondParams,angle,NBond,L);
  //ETersoff = 0.0;
  
  Epot   = ETersoff;
  
  
  //cerr<<"EvdW ="<<EvdW<<endl;
  //cerr<<"Ecoul ="<<Ecoul<<endl;
  //cerr<<"ETersoff ="<<ETersoff<<endl;
  //cerr<<"Epot ="<<Epot<<endl;
  
  return(Epot);
}
