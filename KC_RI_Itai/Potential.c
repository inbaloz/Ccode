

#include "declarations.h"


//*****************************************************************************************************************************
double E_Tersoff(int Natoms,int *Bonded,atom *particle,double *rij,double *Fc, BondParamsStruct *BondParams,double *angle,int *NBond,point L){
  
  double ETersoff=0.0;
  
  //#pragma omp parallel 
  //{
    
    int atomi,atomj,atomk,i,k,typej,typei,typek,index,index2,NBond_j;
    double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
    double Xsi,Tersoff_R,Tersoff_D;
    double Fc_jk,Fc_ij;
    dr_ij rji;
    dr_ij rjk;
    double inside_acos;
    
    //#pragma omp for reduction(+:ETersoff)
    for(atomj=0 ; atomj < Natoms ; atomj++){ // Go over all central atoms j.
      typej  = particle[atomj].type;
      NBond_j=NBond[atomj];
      for(i=0 ; i < NBond_j ; i++){ // Go over all j neighbors i.
	
	atomi = Bonded[atomj*MaxNBond + i];
	//if(particle[atomj].layer != particle[atomi].layer)continue;
	typei = particle[atomi].type;
	index = typei*NAtomParams + typej;
	r_ij = R_ij(particle,atomi,atomj, L);
	//r_ij = rij[atomj*MaxNBond + i];
	Tersoff_n   = BondParams[index].Tersoff_n;
	A     = BondParams[index].Tersoff_A;
	B     = BondParams[index].Tersoff_B;
	beta  = BondParams[index].Tersoff_beta;
	Tersoff_lambda1 = BondParams[index].Tersoff_lambda1;
	Tersoff_lambda2 = BondParams[index].Tersoff_lambda2;
	Tersoff_R = BondParams[index].Tersoff_R;
	Tersoff_D = BondParams[index].Tersoff_D;
	//Fc_ij = Fc[atomj*MaxNBond + i];
	Fc_ij = Fc_(r_ij,Tersoff_R,Tersoff_D);
	Xsi=0.0;
	
	for(k=0 ; k < NBond_j ; k++){ // Go over all j neighbors k != i.
	  if (k!=i){
	    atomk     = Bonded[atomj*MaxNBond + k];
	    typek     = particle[atomk].type;
	    index2    = typek*NAtomParams + typej;
	    Tersoff_c = BondParams[index2].Tersoff_c;
	    Tersoff_d = BondParams[index2].Tersoff_d;
	    Tersoff_h = BondParams[index2].Tersoff_h;
	    Tersoff_R = BondParams[index2].Tersoff_R;
	    Tersoff_D = BondParams[index2].Tersoff_D;
	    r_jk = R_ij(particle,atomk,atomj, L);
	    //r_jk      = rij[atomj*MaxNBond + k];

	    Theta_ijk = calc_val_angle(atomi,atomj,atomk,Natoms,NBond,Bonded,r_jk,r_ij,particle,L,rji,rjk,inside_acos);
	    //Theta_ijk = 0.0;
	    
	    g     = 1+sqr(Tersoff_c)/sqr(Tersoff_d)-sqr(Tersoff_c)/(sqr(Tersoff_d)+sqr(Tersoff_h-cos(Theta_ijk)));
	    
	    //Fc_jk = Fc[atomj*MaxNBond + k];
	    Fc_jk = Fc_(r_jk,Tersoff_R,Tersoff_D);
	    Xsi  += Fc_jk*g;
	  }
	}
	
	//**
	//Xsi=1;
	//**
	//if()
	b_ij      = pow(1+pow(beta*Xsi,Tersoff_n),-1/(2*Tersoff_n));
	
	ETersoff  += Fc_ij*(A*exp(-Tersoff_lambda1*r_ij) - b_ij*B*exp(-Tersoff_lambda2*r_ij));
	//ETersoff  += (A*exp(-Tersoff_lambda1*r_ij) - B*exp(-Tersoff_lambda2*r_ij));
	
      }
      //printf("atomj= %i ETersoff= %.16f,nThread= %i \n",atomj,ETersoff,omp_get_thread_num());
      
    }
    
    //}
 
  ETersoff*=0.5;
  return(ETersoff);
}
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
