#include "declarations.h"

/********************************************************************/
/* Routines required for the calculation of the potential energy.   */
/* Written by Itai Leven and Oded Hod June 2010.                    */
/* (1)ReaxFF - J. Phys. Chem. A 105, 9396-9409 (2001)- hydrocarbons */
/* (2)ReaxFF - J. Phys. Chem. A 107, 3803-3811 (2003)- SiO          */
/* (3)ReaxFF - J. Phys. Chem. A 123, 114703    (2005)- B-N-H        */
/* (4)Suppoting information for the article "Development of the     */
/*  ReaxFF BaZro3 Proton Conductor"                                 */
/********************************************************************/

/**************************************************/
/*                    E_{bond}                    */
/**************************************************/

/********************************************************************/
/* Calculate BO'(_Sigma,_Pi,and _Pi_Pi){ij} from equation (2) of the reference (4) above */

void calc_BO_tag_ij(int typei, int typej, double r_ij, BondParamsStruct *BondParams, double &Sigma, double &Pi, double &Pi_Pi)
{
  // typei      - type of atom i.
  // typej      - type of atom j.
  // particle   - atoms array.
  // L          - Box size for periodic boundary conditions calculation.
  // rij        - ij distance - passed to above function to determine cutoff compliance.
  // BondParams - Array holding the relevant bond parameters.
  
  int PiFlag,PiPiFlag;
  double r0_Sigma,r0_Pi,r0_Pi_Pi;
  double rij_r0_Sigma, rij_r0_Pi,rij_r0_Pi_Pi;
  double pbo1, pbo2, pbo3, pbo4, pbo5, pbo6;
  
  r0_Sigma = BondParams[typei*NAtomParams + typej].r0_Sigma;
  r0_Pi    = BondParams[typei*NAtomParams + typej].r0_Pi;
  r0_Pi_Pi = BondParams[typei*NAtomParams + typej].r0_Pi_Pi;
  PiFlag   = BondParams[typei*NAtomParams + typej].PiFlag;
  PiPiFlag = BondParams[typei*NAtomParams + typej].PiPiFlag;
  pbo1     = BondParams[typei*NAtomParams + typej].pbo1;
  pbo2     = BondParams[typei*NAtomParams + typej].pbo2;
  pbo3     = BondParams[typei*NAtomParams + typej].pbo3;
  pbo4     = BondParams[typei*NAtomParams + typej].pbo4;
  pbo5     = BondParams[typei*NAtomParams + typej].pbo5;
  pbo6     = BondParams[typei*NAtomParams + typej].pbo6;
  
  rij_r0_Sigma = r_ij / r0_Sigma;
  rij_r0_Pi    = r_ij / r0_Pi;
  rij_r0_Pi_Pi = r_ij / r0_Pi_Pi;
  Sigma = (1.0 + BO_cut)*exp(pbo1 * pow(rij_r0_Sigma,pbo2));
  
  Pi    = ((PiFlag == 1)   ? exp(pbo3 * pow(rij_r0_Pi,pbo4))    : 0.0);
  Pi_Pi = ((PiPiFlag == 1)   ? exp(pbo5 * pow(rij_r0_Pi_Pi,pbo6))    : 0.0);
}

/*******************************************************************/
/* Fill the BO'(Sigma,Pi,Pi_Pi) arrays using the calc_BO_Sigma_tag_ij, calc_BO_Pi_tag_ij and the calc_BO_Pi_Pi_tag_ij functions from above */

void calc_BO_tag(Bond_order *BO_tag, int Natoms,atom *particle, BondParamsStruct *BondParams, int *NBond, int *Bonded, double *rij)
{
  int atomi, j;
  int atomj,typei,typej;
  double Sigma, Pi, Pi_Pi,BOp_tot;
  double r_ij;
  
  for (atomi=0 ; atomi < Natoms ; atomi++){
    typei  = particle[atomi].type;
    for (j=0 ; j < NBond[atomi] ; j++){
      atomj  = Bonded[atomi*MaxNBond + j];
      typej  = particle[atomj].type;
      r_ij   = rij[atomi*Natoms + atomj];
      
      calc_BO_tag_ij(typei,typej,r_ij,BondParams,Sigma,Pi,Pi_Pi);
      
      BOp_tot = Sigma + Pi + Pi_Pi;
      
      if (BOp_tot >= BO_cut){     //this term appears in LAAMPS code , when its on there are disscontinuities in the energy
	BO_tag[atomi*MaxNBond + j].Sigma = Sigma - BO_cut;
	BO_tag[atomi*MaxNBond + j].Pi    = Pi;
	BO_tag[atomi*MaxNBond + j].Pi_Pi = Pi_Pi;
	BO_tag[atomi*MaxNBond + j].tot   = Sigma + Pi + Pi_Pi;
      }
      else{
	BO_tag[atomi*MaxNBond + j].Sigma = 0.0;
	BO_tag[atomi*MaxNBond + j].Pi    = 0.0;
	BO_tag[atomi*MaxNBond + j].Pi_Pi = 0.0;
	BO_tag[atomi*MaxNBond + j].tot   = 0.0;
      }
    }
  }
}

/**************************************************************/
/* Calculate Delta' from equation (3a) of the reference (4) above  */

void calc_Delta_tag(double *delta_tag, int Natoms, atom *particle, Bond_order *BO_tag, int *NBond)
{
  int atomi, j;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    delta_tag[atomi] = -particle[atomi].val;
    
    for(j=0 ; j < NBond[atomi] ; j++) delta_tag[atomi] += BO_tag[atomi*MaxNBond + j].tot;
  }
}
/***********************************************************************************/
/* Calculate Delta_boc from equation (3b) of the reference (4) above without the primes */

void calc_Delta_tag_boc(double *delta_tag_boc, int Natoms, atom *particle, Bond_order *BO_tag, int *NBond, AtomParamsStruct *AtomParams)
{
  int atomi, j;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    delta_tag_boc[atomi] = -AtomParams[particle[atomi].type].Val_boc;
    
    for(j=0 ; j < NBond[atomi] ; j++) delta_tag_boc[atomi] += BO_tag[atomi*MaxNBond + j].tot;
  }
}

/****************************************************************/
/* Calculate F1_{ij} from equation (3b) of the reference (1) above  */

double calc_F1_ij(int atomi, int j,int atomj ,double Delta_tag_i, double Delta_tag_j, atom *particle)
{
  double F1, F2_ij, F3_ij;
  
  F2_ij = calc_F2_ij(Delta_tag_i,Delta_tag_j);
  F3_ij = calc_F3_ij(Delta_tag_i,Delta_tag_j);
  
  F1 = 0.5 * ( ((particle[atomi].val + F2_ij) / (particle[atomi].val + F2_ij + F3_ij)) + ((particle[atomj].val + F2_ij) / (particle[atomj].val + F2_ij + F3_ij)) );
  
  return(F1);
}

/***********************************************************/
/* Calculate F1 and save it to an array.                   */

void calc_F1(atom *particle, double *f1, int *NBond, int *Bonded, double *delta_tag, int Natoms, BondParamsStruct *BondParams)
{
  int atomi, j, atomj,F1_Flag,typei,typej;
  double Delta_tag_i, Delta_tag_j;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    typei = particle[atomi].type;
    Delta_tag_i = delta_tag[atomi];
    for(j=0 ; j < NBond[atomi] ; j++){
      atomj = Bonded[atomi*MaxNBond + j];
      typej = particle[atomj].type;
      Delta_tag_j = delta_tag[atomj];
      F1_Flag = BondParams[typei*NAtomParams + typej].F1_Flag;
      
      if(F1_Flag == 1) f1[atomi*MaxNBond + j] = calc_F1_ij(atomi,j,atomj,Delta_tag_i,Delta_tag_j,particle);
      else{
	f1[atomi*MaxNBond + j] = 1.0;
      }
    }
  }
}
/*****************************************************************/
  /* Calculate F4/F5 from equation (3e/3f) of the reference (1) above  */

inline double calc_F45_ij(double Delta_tag_i, double BO_tag, double Pboc3, double Pboc4, double Pboc5)
{
  double F45;
  
  F45 = 1.0 / (1.0 + exp(-Pboc3 * (Pboc4 * sqr(BO_tag) - Delta_tag_i) + Pboc5));
  
  return(F45);
}

/*****************************************************************************/
/* Calculate F4/F5 and save it to an array.                */

void calc_F45(int Natoms, double *delta_tag_boc,  Bond_order *BO_tag, double *f45, int *NBond, int *Bonded, BondParamsStruct *BondParams, atom *particle, int *BondNum)
{
  int atomi, j, atomj,atomtype_i,atomtype_j;
  double BO_tag_ij, Delta_tag_i,Delta_tag_j;
  double Pboc3,Pboc4,Pboc5;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    atomtype_i  = particle[atomi].type;
    Delta_tag_i = delta_tag_boc[atomi];
    for(j=0 ; j < NBond[atomi] ; j++){
      atomj = Bonded[atomi*MaxNBond + j];
      BO_tag_ij   = BO_tag[atomi*MaxNBond + j].tot;
      Delta_tag_j = delta_tag_boc[atomj];
      atomtype_j  = particle[atomj].type;
      Pboc3   = BondParams[atomtype_i*NAtomParams + atomtype_j].Pboc3;
      Pboc4   = BondParams[atomtype_i*NAtomParams + atomtype_j].Pboc4;
      Pboc5   = BondParams[atomtype_i*NAtomParams + atomtype_j].Pboc5;
      
      f45[atomj*MaxNBond + BondNum[atomj*Natoms + atomi]] = calc_F45_ij(Delta_tag_j,BO_tag_ij,Pboc3,Pboc4,Pboc5);
      
      f45[atomi*MaxNBond + j] = calc_F45_ij(Delta_tag_i,BO_tag_ij,Pboc3,Pboc4,Pboc5);
    }
  }
}

/******************************************************************************/
/* Calculate BO_{ij} seperatley for Sigma, Pi and Pi_Pi from equation (3a) of the reference (1) above  */

void calc_BO(Bond_order *BO, Bond_order *BO_tag, double *f1, double *f45, int *NBond, int *Bonded, int *BondNum, int Natoms, atom *particle)
{
  int atomi, j, atomj;
  int index;
  double F1, F4, F5;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    for(j=0 ; j < NBond[atomi] ; j++){
      
      atomj = Bonded[atomi*MaxNBond + j];
      index = atomi*MaxNBond + j;
      F1 = f1[atomi*MaxNBond + j];
      F4 = f45[atomi*MaxNBond + j];
      F5 = f45[atomj*MaxNBond + BondNum[atomj*Natoms + atomi]];
      
      //********************new BO**********************
      BO[index].Pi    = BO_tag[index].Pi    * sqr(F1) * F4 * F5;
      BO[index].Pi_Pi = BO_tag[index].Pi_Pi * sqr(F1) * F4 * F5;
      BO[index].tot   = BO_tag[index].tot   * F1 * F4 * F5;
      BO[index].Sigma = BO[index].tot - (BO[index].Pi + BO[index].Pi_Pi);
      //***********************************************
      
      if( BO[index].tot < 1e-10 )
	BO[index].tot = 0.0;
      if( BO[index].Sigma < 1e-10 )
	BO[index].Sigma = 0.0;
      if( BO[index].Pi < 1e-10 )
	BO[index].Pi = 0.0;
      if( BO[index].Pi_Pi < 1e-10 )
	BO[index].Pi_Pi = 0.0;
      
    }
  }
}

/**********************************************************************************/
/*sums all of the Ebond_ij terms to get the total Ebond                           */

double calc_Ebond(int Natoms, int *NBond, int *Bonded,atom *particle, BondParamsStruct *BondParams, Bond_order *BO)
{
  int atomi, j, atomj, typei, typej;
  double pbe1, pbe2, D_Sigma, D_Pi, D_Pi_Pi, BO_Sigma, BO_Pi, BO_Pi_Pi;
  double Ebond;
  
  Ebond = 0.0;
  
  for (atomi=0 ; atomi < Natoms ; atomi++){
    typei    = particle[atomi].type;
    for (j=0 ; j < NBond[atomi] ; j++){
      atomj    = Bonded[atomi*MaxNBond + j];
      typej    = particle[atomj].type;
      if(atomj < atomi){
	D_Sigma  = BondParams[typei*NAtomParams + typej].D_Sigma;
	D_Pi     = BondParams[typei*NAtomParams + typej].D_Pi;
	D_Pi_Pi  = BondParams[typei*NAtomParams + typej].D_Pi_Pi;
	pbe1     = BondParams[typei*NAtomParams + typej].pbe1;
	pbe2     = BondParams[typei*NAtomParams + typej].pbe2;
	BO_Sigma = BO[atomi*MaxNBond + j].Sigma;
	BO_Pi    = BO[atomi*MaxNBond + j].Pi;
	BO_Pi_Pi = BO[atomi*MaxNBond + j].Pi_Pi;
	Ebond   -= (D_Sigma * BO_Sigma * exp(pbe1 * (1.0 - pow(BO_Sigma, pbe2))) + D_Pi * BO_Pi  + D_Pi_Pi * BO_Pi_Pi);
      }
    }
  }
  return(Ebond);
}

/**************************************************/
/*                    E_{lp}                    */
/**************************************************/

/**********************************************************************************/
// Calculates nlp_i,delta_e_i and delta_lp as in equation (4) and (5) in reference (2)

void calc_Delta_lp(int Natoms, atom *particle, int *NBond, Bond_order *BO, double *delta_lp,double *delta_e, AtomParamsStruct *AtomParams,double *n_lp)
{
  int j,atomi;
  double Delta_e_i;
  double int_2;
  
  for (atomi=0; atomi < Natoms ; atomi++){
    delta_e[atomi]  = -particle[atomi].OS_Elec;
    for (j=0 ; j < NBond[atomi] ; j++)delta_e[atomi] += BO[atomi*MaxNBond + j].tot;
    int_2       = int(delta_e[atomi]/2.0);
    n_lp[atomi] = exp(-lambda16 * sqr(2.0 + delta_e[atomi] - 2.0 * int_2)) - int_2;
    
    delta_lp[atomi] = particle[atomi].lp - n_lp[atomi];
    
  }
}

/**********************************************************************************/
/* Calculates E_lp from equation (5) in the reference (2) above                    */

double  calc_E_lp(int Natoms, AtomParamsStruct *AtomParams, int *NBond, double *delta_lp,atom *particle)
{
  int atomi, atomj, j;
  double E_lp, Plp2;
  
  E_lp = 0.0;
  
  for (atomi=0 ; atomi < Natoms ; atomi++){
    Plp2   = AtomParams[particle[atomi].type].Plp2;
    
    E_lp += Plp2 * delta_lp[atomi] / (1.0 + exp(-75.0 * delta_lp[atomi]));
    
  }
  return(E_lp);
}

/**************************************************/
/*                    E_{over}                    */
/**************************************************/

/* Calculate Delta from equation (4) of the reference (1) above without the primes */

void calc_Delta_boc(double *delta_boc, int Natoms, Bond_order *BO, int *NBond, AtomParamsStruct *AtomParams,atom *particle)
{
  int atomi, j;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    delta_boc[atomi] = -AtomParams[particle[atomi].type].Val_boc;
    
    for(j=0 ; j < NBond[atomi] ; j++) delta_boc[atomi] += BO[atomi*MaxNBond + j].tot;
  }
}

/***********************************************************************************/

/* Calculate Delta from equation (4) of the reference (1) above without the primes */

void calc_Delta(double *delta, int Natoms, atom *particle, Bond_order *BO, int *NBond)
{
  int atomi, j;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    delta[atomi] = -particle[atomi].val;
    
    for(j=0 ; j < NBond[atomi] ; j++) delta[atomi] += BO[atomi*MaxNBond + j].tot;
  }
}

/**********************************************************************************/
/* Calculate SOV from the equation (6c) from the reference (2) above              */

double SOV(int atomi, int *Bonded, double *delta, Bond_order *BO, double *delta_lp, int *NBond)
{
  int atomj, j;
  double SOV; 
  
  SOV = 0.0;
  
  for (j=0; j < NBond[atomi]; j++){
    atomj = Bonded[atomi*MaxNBond + j];
    SOV  += (delta[atomj] - delta_lp[atomj]) * (BO[atomi*MaxNBond + j].Pi + BO[atomi*MaxNBond + j].Pi_Pi);
  }
  return (SOV);
}

//********************************************************************************
/* Calculate Delta_i_tag from the equation (6b) in the reference (2) above      */

double Delta_lpcorr_i(int atomi, int *Bonded, Bond_order *BO, double *delta, double *delta_lp,int *NBond)
{
  double delta_lpcorr_i, SOV_i;
  SOV_i          = SOV(atomi,Bonded,delta,BO,delta_lp,NBond);
  delta_lpcorr_i = delta[atomi] - delta_lp[atomi] / (1.0 + Povun3 * exp(Povun4 * SOV_i));
  return (delta_lpcorr_i);
}

//********************************************************************************
/* Calculate E_over_i from the equation (6a) in the reference (2) above         */

double Eover_i(int atomi, atom *particle, BondParamsStruct *BondParams, Bond_order *BO_tag, int Natoms, int *NBond, int *Bonded, Bond_order *BO, double *delta, double *delta_lp, AtomParamsStruct *AtomParams)
{
  int j,atomj,typei,typej;
  double Sum_Povon1_BO, Eover_i;
  double delta_lpcorr_i,Povun1,D_Sigma,Povun2;
  
  Eover_i = 0.0;
  Sum_Povon1_BO = 0.0;
  typei  = particle[atomi].type;
  Povun2 = AtomParams[typei].Povun2;
  
  
  for (j=0 ; j < NBond[atomi] ; j++){
    atomj    = Bonded[atomi*MaxNBond + j];
    typej    = particle[atomj].type;
    Povun1   = BondParams[typei*NAtomParams + typej].Povun1;
    D_Sigma  = BondParams[typei*NAtomParams + typej].D_Sigma;
    Sum_Povon1_BO += Povun1 * D_Sigma * BO[atomi*MaxNBond + j].tot;
  } 
  
  delta_lpcorr_i = Delta_lpcorr_i(atomi,Bonded,BO,delta,delta_lp,NBond);
  
  Eover_i = (Sum_Povon1_BO * delta_lpcorr_i) / (delta_lpcorr_i + particle[atomi].val + 1e-8) * (1.0/(1.0 + exp(Povun2 * delta_lpcorr_i)));
  return(Eover_i);
}

//********************************************************************************
/* sums all of the E_over_i terms to get the total E_over                       */

double calc_Eover(atom *particle, BondParamsStruct *BondParams, Bond_order *BO_tag, int Natoms, int *NBond, int *Bonded, Bond_order *BO, double *delta, double *delta_lp, AtomParamsStruct *AtomParams)
{
  int atomi,typei;
  double Pbe3, Eover;
  double delta_lpcorr_i;
  
  Eover = 0.0;
  
  for (atomi=0 ; atomi < Natoms ; atomi++){
    Eover += Eover_i(atomi,particle,BondParams,BO_tag,Natoms,NBond, Bonded,BO,delta,delta_lp,AtomParams);
  }
  return(Eover);
}

/**************************************************/
/*                    E_{under}                   */
/**************************************************/

/*****************************************************************************/
/* Calculate f6 from equation (7b) of the reference (1) above                */

double calc_f6_i(int atomi, Bond_order *BO, double *delta, int *NBond, int *Bonded, double *delta_lp)
{
  double f6_i, SOV_i;
  
  SOV_i = SOV(atomi,Bonded,delta,BO,delta_lp,NBond);
  f6_i  = 1.0 / (1.0 + lambda9 * exp(lambda10 * SOV_i));
  
  return(f6_i);
}

/***************************************************************************************/
/* Calculate the total Eunder from equation (7a) of the reference (1) above            */

double calc_Eunder(int Natoms, atom *particle, AtomParamsStruct *AtomParams, double *delta, Bond_order *BO, int *NBond, int *Bonded,double *delta_lp)
{
  int atomi, typei;
  double Eunder, Povun5,Povun2, f6_i,delta_i_lp_corr;
  
  //calc_Delta(delta,Natoms,particle,BO,NBond); - This is calculated outside and passed.
  
  Eunder = 0.0;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    typei = particle[atomi].type;
    if((typei == 5 || typei == 7)){ // Eunder is relevant only for B-N,N-N bonds.
      Povun5 = AtomParams[typei].Povun5;
      Povun2 = AtomParams[typei].Povun2;
      f6_i = calc_f6_i(atomi,BO,delta,NBond,Bonded,delta_lp);
      delta_i_lp_corr  = Delta_lpcorr_i(atomi,Bonded,BO,delta,delta_lp,NBond);
      Eunder -= Povun5 * ((1.0 - exp(lambda7*delta_i_lp_corr)) / (1.0 + exp(-Povun2*delta_i_lp_corr))) * f6_i;
    }
  }
  return(Eunder);
}

/**************************************************/
/*                    E_{vdW}                     */
/**************************************************/
double calc_Tap(double r_ij)
{
  double Tap;
  
  Tap = 0.0;
  
  Tap = Tap_7 * r_ij + Tap_6;
  Tap = Tap * r_ij   + Tap_5;
  Tap = Tap * r_ij   + Tap_4;
  Tap = Tap * r_ij   + Tap_3;
  Tap = Tap * r_ij   + Tap_2;
  Tap = Tap * r_ij   + Tap_1;
  Tap = Tap * r_ij   + Tap_0;
  
  return(Tap);
}

/************************************************************************************************/
/* Calculate the tatal E_vdW (Van der Waals) from the equation (12a) in the reference (1) above */       

double calc_EvdW(int Natoms, double *rij, atom *particle, BondParamsStruct *BondParams)
{
  int atomi, atomj, typei, typej, index;
  double EvdW, lambdaW, alpha, rvdW, Epsilon, f13, term,Inner_shield;
  double ecore, acore, rcore, Tap, r_ij;
  
  EvdW = 0.0;
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    typei = particle[atomi].type;
    
    for(atomj=0 ; atomj < atomi ; atomj++){
      typej = particle[atomj].type;
      
      index = typei*NAtomParams + typej;
      
      //calculating Tap
      r_ij = rij[atomi*Natoms + atomj];
      if (r_ij < non_bond_cut){
	Tap = calc_Tap(r_ij);
	if (particle[atomj].layer == particle[atomi].layer){
	  alpha   = BondParams[index].alpha;
	  rvdW    = BondParams[index].rvdw;
	  Epsilon = BondParams[index].Epsilon;
	}
	else{
	  alpha   = BondParams[index].alpha_long;
	  rvdW    = BondParams[index].rvdw_long;
	  Epsilon = BondParams[index].Epsilon_long;
	}
	lambdaW = BondParams[index].lambdaW;
	ecore   = BondParams[index].ecore;
	acore   = BondParams[index].acore;
	rcore   = BondParams[index].rcore;
	
	f13 = pow( (pow(r_ij , lambda29) + pow((1.0 / lambdaW) , lambda29)) , (1.0 / lambda29) );
	
	term = alpha * (1.0 - (f13 / rvdW));
	Inner_shield = ecore * exp(acore * (1.0-(r_ij/rcore)));
	
	//EvdW += Tap * (Epsilon * (exp(term) - 2.0 * exp(0.5 * term)) + Inner_shield);
	
	if (particle[atomi].type == 5 && particle[atomj].type == 5)
	  {
	    //EvdW += Tap * (Epsilon * (exp(term) - 901.1466/(pow(r_ij,6)) ) + Inner_shield);
	    EvdW += (Epsilon * (exp(term)) - 901.1466/(pow(r_ij,6)) );
	  }
	else if ((particle[atomi].type == 7 && particle[atomj].type == 5) ||(particle[atomi].type == 5 && particle[atomj].type == 7))
	  {
	    //EvdW += Tap * (Epsilon * (exp(term) - 442.3059/(pow(r_ij,6)) ) + Inner_shield);
	    EvdW += (Epsilon * (exp(term)) - 442.3059/(pow(r_ij,6) ));
	  }
	else if (particle[atomi].type == 7 && particle[atomj].type == 7)
	  {
	    //EvdW += Tap * (Epsilon * (exp(term) - 257.6673/(pow(r_ij,6)) ) + Inner_shield);
	    EvdW += (Epsilon * (exp(term)) - 257.6673/(pow(r_ij,6) ));
	  }
      }
    }
  }
  return(EvdW);
}

/*************************************************************/
/* Calculate E_coulomb                                       */

double calc_Ecoulomb(int Natoms, atom *particle, double *rij, double *Bmat, double *qvec, double *Avec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rijShiled_1, FILE* logfile)
{
  int atomi, atomj;
  double Ecoulomb, Tap, r_ij;
  double rijShield,shield,gamma;
  Ecoulomb=0.0;
  
  Calc_Charge(particle,Natoms,Bmat,qvec,Avec,AtomParams,BondParams,rij,rijShiled_1,logfile);
  
  
  for (atomi = 0;atomi < Natoms;atomi++){
    if (fabs(qvec[atomi]) > 1.5){
      cerr<<"Error in Ecoulomb!!! charge on atom["<<atomi<<"]"<<" is "<<qvec[atomi]<<", bigger then 1.5 - ending session"<<endl;
      exit(0);
    }
  }
  for(atomi=0 ; atomi < Natoms ; atomi++){
    for(atomj=0 ; atomj < atomi ; atomj++){
      r_ij = rij[atomi*Natoms + atomj];
      if (r_ij < non_bond_cut){
	gamma = BondParams[particle[atomi].type*NAtomParams + particle[atomj].type].original_gamma;
	shield    = pow(rij[atomi*Natoms + atomj],3.0) + pow(1.0/gamma,3.0);
	rijShield = 1.0/pow(shield,0.333333333333);
	Tap = calc_Tap(r_ij);
	Ecoulomb += Tap*(qvec[atomi] * qvec[atomj] * rijShield);
      }
    }
  }
  
  Ecoulomb *= (Kappa * eV2kcalmole);
  
  return(Ecoulomb);
}
/**************************************************/
/*                    E_{val}                     */
/**************************************************/

/**************************************************************/
/* Calculate f7 using Eq. 8b of the above mentioned paper.    */

double calc_f7(int atomi,int j,Bond_order *BO,double Pval_4,double Pval_3)
{
  double f7,BO_ij;
  
  BO_ij = BO[atomi*MaxNBond + j].tot;
  f7 = (1.0 - exp(-Pval_3 * pow(BO_ij, Pval_4)));
  return(f7);
}
/*********************************************************************/
/* Calculate $Theta_0$ using Eq. 8d of the above mentioned paper.    */

double calc_Theta_0(int atomi, int atomj, int atomk, double *delta_boc, atom *particle, AngleParamsStruct *AngleParams, int *NBond, Bond_order *BO, double *n_lp,double *delta_e)
{
  int typei, typej, typek, index, n,int_2;
  double SBO, SBO2, Theta0, Theta00;
  double exp1, exp2, Ferm1, Ferm2,delta_boc_j;
  double Sum_BO_Pi,Sum_Exp_BO_8, BO_8,nlp;
  delta_boc_j = delta_boc[atomj];
  typei = particle[atomi].type;
  typej = particle[atomj].type;
  typek = particle[atomk].type;
  Sum_BO_Pi    = 0.0;
  Sum_Exp_BO_8 = 1.0;
  
  for (n=0; n < NBond[atomj]; n++){
    Sum_BO_Pi    += BO[atomj*MaxNBond + n].Pi + BO[atomj*MaxNBond + n].Pi_Pi;
    BO_8          = pow(BO[atomj*MaxNBond + n].tot,6.0);//maybe **6
    Sum_Exp_BO_8 *= exp(-BO_8);//maybe **6
  }
  
  int_2 = int(delta_e[atomj]/2.0);
  
  if ((delta_e[atomj] - 2.0 * int_2) >= 0.0)nlp = 0.0;
  else nlp = n_lp[atomj];
  
  SBO = Sum_BO_Pi + (1.0 - Sum_Exp_BO_8) * (-delta_boc_j - Pval_8* nlp);
  
  //**************************
  
  exp1  = exp(-Wd*(SBO-C1));
  exp2  = exp(-Wd*(SBO-C2));
  Ferm1 = 1.0 / (1.0 + exp1);
  Ferm2 = 1.0 / (1.0 + exp2);
  SBO2 = fabs(SBO * Ferm1 + (2.0 - SBO) * Ferm2);
  
  //**************************
  
  index = typei*NAtomParams*NAtomParams + typej*NAtomParams + typek;
  Theta00 = AngleParams[index].Theta00;
  Theta0 = PIE - Theta00 * (1.0 - exp(-Pval_10 * (2.0 - SBO2)));
  
  return(Theta0);
}

//******************************************************
/* Calculate E_val using Eqs. 8 of the above mentioned paper. */

double calc_Eval(int Natoms, atom *particle, int *NBond, int *Bonded, double *delta_boc, Bond_order *BO, double *angle, AngleParamsStruct *AngleParams,double *n_lp, AtomParamsStruct *AtomParams,double *delta_e)
{
  int atomi, atomj, atomk, i, k, typei, typej, typek, index;
  double f7ji, f7jk, f8ijk, Pval_1, Pval_2, Theta0, Thetaijk, EVal;
  double Pval_7,Pval_4,Pval_3;
  double Epen, BO_ij, BO_jk;
  double ExpPval;
  double BOA_jk,BOA_ij;
  EVal = 0.0;
  
  for(atomj=0 ; atomj < Natoms ; atomj++){ // Go over all central atoms j.
    typej  = particle[atomj].type;
    Pval_3 = AtomParams[typej].Pval_3;
    for(i=0 ; i < NBond[atomj] ; i++){ // Go over all j neighbors i.
      atomi = Bonded[atomj*MaxNBond + i];
      typei = particle[atomi].type;
      BO_ij = BO[atomj*MaxNBond + i].tot;
      for(k=0 ; k < i ; k++){ // Go over all j neighbors k != i.
	atomk = Bonded[atomj*MaxNBond + k];
	typek = particle[atomk].type;
	index = typei*NAtomParams*NAtomParams + typej*NAtomParams + typek;
	BO_jk = BO[atomj*MaxNBond + k].tot;
	BOA_jk = BO_jk - thb_cut;
	BOA_ij = BO_ij - thb_cut;
	
	if( (BOA_ij > 0.0) && (BOA_jk > 0.0) &&  //this term is appears in LAAMPs code turning it on makes Energy discontinuities 
	    (BO_ij > thb_cut) &&  
	    (BO_jk > thb_cut) &&  
	    (BO_ij * BO_jk > BO_cut) ) {
	  
	  Pval_1   = AngleParams[index].Pval_1;
	  Pval_2   = AngleParams[index].Pval_2;
	  Pval_7   = AngleParams[index].Pval_7;
	  Pval_4   = AngleParams[index].Pval_4;
	  f7ji  = calc_f7(atomj,i,BO,Pval_4,Pval_3);
	  f7jk  = calc_f7(atomj,k,BO,Pval_4,Pval_3);
	  f8ijk = calc_f8(atomi,atomj,atomk,Natoms,delta_boc,particle,AtomParams,Pval_7);
	  Theta0   = calc_Theta_0(atomi,atomj,atomk,delta_boc,particle,AngleParams,NBond,BO,n_lp,delta_e);
	  Thetaijk = angle[atomj * MaxNBond * MaxNBond + i * MaxNBond + k];
	  ExpPval = exp(-Pval_2 * sqr(Theta0 - Thetaijk));
	  ExpPval = Pval_1 * (1.0 - ExpPval);
	  
	  
	  EVal += (f7ji * f7jk * f8ijk * ExpPval);
	}
      }
    }
  }
  return(EVal);
}

/**************************************************/
/*                  E_{Torsion}                   */
/**************************************************/
void calc_Exp_Ptor_2_array(atom *particle, int Natoms, Bond_order *BO, int *NBond, double *Exp_Ptor_2, int *Bonded, int *BondNum)
{
  double BO_ij;
  int atomi, j, atomj;
  
  for (atomi=0; atomi < Natoms ; atomi++){
    for (j=0; j < NBond[atomi]; j++){
      atomj = Bonded[atomi*MaxNBond + j];
      BO_ij = BO[atomi*MaxNBond + j].tot;
      Exp_Ptor_2[atomi*MaxNBond + j] = Exp_Ptor_2[atomj*MaxNBond + BondNum[atomj*Natoms + atomi]] =  1.0 - exp(-Ptor_2 * BO_ij);
    }
  }
}

//**************************************************************************************************************************

void calc_F11(double *delta_boc, atom *particle, int Natoms, int *NBond, int *Bonded, double *f11)
{
  double Delta_Boc_j, Delta_Boc_k, Exp_Ptor_3;
  int atomj,atomk,k;
  for (atomj=0; atomj < Natoms ; atomj++){
    for (k=0; k < NBond[atomj]; k++){
      atomk       = Bonded[atomj*MaxNBond + k];
      Delta_Boc_j = delta_boc[atomj];
      Delta_Boc_k = delta_boc[atomk];
      Exp_Ptor_3  = exp(-Ptor_3*(Delta_Boc_j  + Delta_Boc_k));
      f11[atomj*MaxNBond + k]  = (2.0 + Exp_Ptor_3) / (1.0 + Exp_Ptor_3 + exp(Ptor_4*(Delta_Boc_j + Delta_Boc_k)));
    }
  }
}
void calc_Exp_Pcot_2_array(atom *particle, int Natoms, Bond_order *BO, int *NBond, double *Exp_Pcot_2, int *Bonded, int *BondNum){
  
  int atomi, j, atomj;
  double BO_ij;  
  for (atomi=0; atomi < Natoms ; atomi++){
    for (j=0; j < NBond[atomi]; j++){
      atomj = Bonded[atomi*MaxNBond + j];
      BO_ij = BO[atomi*MaxNBond + j].tot;
      Exp_Pcot_2[atomi*MaxNBond + j] = Exp_Pcot_2[atomj*MaxNBond + BondNum[atomj*Natoms + atomi]] =  exp(-Pcot_2 * sqr(BO_ij - 1.5));
    }
  }
}

//***************************************************************************************************************************

double calc_Etor_conj(atom *particle,double *delta_boc, Bond_order *BO, int *NBond, int *Bonded, int *BondNum, double *angle, int Natoms,dihedralParamsStruct *dihedralParams, double *Exp_Pcot_2, double *Exp_Ptor_2, double *f11,point L)
{
  int i, k, l, atomi, atomj, atomk, atoml;
  int index1,typej,typek;
  double Etor, Econj,V1, V2, V3, Ptor_1;
  double Thetaijk, Thetajkl, Omega_ijkl, sin_ijk, sin_jkl, sin_ijkl;
  double Delta_Boc_j, Delta_Boc_k, BOij, BOjk, BOkl;
  double F11,F10,F12, Exp_Ptor_2_ij, Exp_Ptor_2_jk, Exp_Ptor_2_kl;
  double Exp_Pcot_2_ij,Exp_Pcot_2_jk,Exp_Pcot_2_kl;
  double dampij, dampjk, dampkl,Pcot_1,BOjk_Pi;
  
  Etor  = 0.0;
  Econj = 0.0;
  
  for(atomj=0 ; atomj < Natoms ; atomj++){
    typej = particle[atomj].type;
    Delta_Boc_j = delta_boc[atomj];
    for(k=0 ; k < NBond[atomj] ; k++){
      atomk = Bonded[atomj*MaxNBond + k];
      typek = particle[atomk].type;
      if(atomk < atomj){ // Makes sure that each central bond jk appears only once.
	Delta_Boc_k = delta_boc[atomk];
	BOjk    = BO[atomj*MaxNBond + k].tot;
	BOjk_Pi = BO[atomj*MaxNBond + k].Pi;
	for(i=0 ; i < k ; i++){

	  atomi = Bonded[atomj*MaxNBond + i];
	  
	  BOij = BO[atomj*MaxNBond + i].tot;
	  Thetaijk = angle[atomj * MaxNBond * MaxNBond + i * MaxNBond + k];
	  sin_ijk  = sin(Thetaijk);
	  
	  for(l=0 ; l < NBond[atomk] ; l++){
	    atoml = Bonded[atomk*MaxNBond + l];
	    
	    if((atoml != atomj) && (atoml != atomi)){
	      BOkl = BO[atomk*MaxNBond + l].tot;
	      
	      if ( (BOij > thb_cut) && (BOjk > thb_cut) && (BOkl > thb_cut) ){ //this term used in lammps makes discontinuities in energy comservation!!!!
		
		Thetajkl = angle[atomk * MaxNBond * MaxNBond + BondNum[atomk*Natoms + atomj] * MaxNBond + l];
		sin_jkl  = sin(Thetajkl);
		
		Omega_ijkl = calc_dihedral(atomi,atomj,atomk,atoml,particle,L);
		sin_ijkl = sin(Omega_ijkl);
		
		index1   = typej*NAtomParams + typek;
		Ptor_1   = dihedralParams[index1].Ptor_1;
		V1       = dihedralParams[index1].V1;
		V2       = dihedralParams[index1].V2;
		V3       = dihedralParams[index1].V3;
		
		Exp_Ptor_2_ij = Exp_Ptor_2[atomj*MaxNBond + i];
		Exp_Ptor_2_jk = Exp_Ptor_2[atomj*MaxNBond + k];
		Exp_Ptor_2_kl = Exp_Ptor_2[atomk*MaxNBond + l];
		
		Exp_Pcot_2_ij = Exp_Pcot_2[atomj*MaxNBond + i];
		Exp_Pcot_2_jk = Exp_Pcot_2[atomj*MaxNBond + k];
		Exp_Pcot_2_kl = Exp_Pcot_2[atomk*MaxNBond + l];
		
		dampij = 1.0 / (1.0 + exp(-Wd3 * (BOij - C3)));
		dampjk = 1.0 / (1.0 + exp(-Wd3 * (BOjk - C3)));
		dampkl = 1.0 / (1.0 + exp(-Wd3 * (BOkl - C3)));
		
		F10      = Exp_Ptor_2_ij * Exp_Ptor_2_jk * Exp_Ptor_2_kl;
		F11      = f11[atomj*MaxNBond + k];
		F12      = Exp_Pcot_2_ij * Exp_Pcot_2_jk * Exp_Pcot_2_kl * dampij * dampjk *dampkl;
		
		Pcot_1   = dihedralParams[typej*NAtomParams + typek].Pcot_1;
		
		Etor += (F10 * sin_ijk * sin_jkl * 0.5 * (V1 * (1.0 + cos(Omega_ijkl)) +  V2 * exp(Ptor_1 * sqr(-BOjk_Pi + 2.0 - F11)) * (1.0 - cos(2.0 * Omega_ijkl)) + V3 * (1.0 + cos(3.0 * Omega_ijkl))));
		Econj += F12 * Pcot_1 * (1.0 - sqr(sin_ijkl) * sin_ijk * sin_jkl);
	      }
	    }
	  }
	}
	
	for(i=k+1 ; i < NBond[atomj] ; i++){
	  atomi = Bonded[atomj*MaxNBond + i];
	  
	  BOij = BO[atomj*MaxNBond + i].tot;
	  Thetaijk = angle[atomj * MaxNBond * MaxNBond + i * MaxNBond + k];
	  sin_ijk  = sin(Thetaijk);
	  
	  for(l=0 ; l < NBond[atomk] ; l++){
	    atoml = Bonded[atomk*MaxNBond + l];
	    
	    if((atoml != atomj) && (atoml != atomi)){
	      BOkl = BO[atomk*MaxNBond + l].tot;

	      if ( (BOij > thb_cut) && (BOjk > thb_cut) && (BOkl > thb_cut) ){ //this term used in lammps!!!!
		
		Thetajkl = angle[atomk * MaxNBond * MaxNBond + BondNum[atomk*Natoms + atomj] * MaxNBond + l];
		sin_jkl  = sin(Thetajkl);
		
		Omega_ijkl = calc_dihedral(atomi,atomj,atomk,atoml,particle,L);
		sin_ijkl = sin(Omega_ijkl);
		
		index1   = typej*NAtomParams + typek;
		Ptor_1   = dihedralParams[index1].Ptor_1;
		V1       = dihedralParams[index1].V1;
		V2       = dihedralParams[index1].V2;
		V3       = dihedralParams[index1].V3;
		
		Exp_Ptor_2_ij = Exp_Ptor_2[atomj*MaxNBond + i];
		Exp_Ptor_2_jk = Exp_Ptor_2[atomj*MaxNBond + k];
		Exp_Ptor_2_kl = Exp_Ptor_2[atomk*MaxNBond + l];
		
		Exp_Pcot_2_ij = Exp_Pcot_2[atomj*MaxNBond + i];
		Exp_Pcot_2_jk = Exp_Pcot_2[atomj*MaxNBond + k];
		Exp_Pcot_2_kl = Exp_Pcot_2[atomk*MaxNBond + l];
		
		dampij = 1.0 / (1.0 + exp(-Wd3 * (BOij - C3)));
		dampjk = 1.0 / (1.0 + exp(-Wd3 * (BOjk - C3)));
		dampkl = 1.0 / (1.0 + exp(-Wd3 * (BOkl - C3)));
		
		F10      = Exp_Ptor_2_ij * Exp_Ptor_2_jk * Exp_Ptor_2_kl;
		F11      = f11[atomj*MaxNBond + k];
		F12      = Exp_Pcot_2_ij * Exp_Pcot_2_jk * Exp_Pcot_2_kl * dampij * dampjk *dampkl;
		
		Pcot_1   = dihedralParams[typej*NAtomParams + typek].Pcot_1;
		
		Etor += (F10 * sin_ijk * sin_jkl * 0.5 * (V1 * (1.0 + cos(Omega_ijkl)) +  V2 * exp(Ptor_1 * sqr(-BOjk_Pi + 2.0 - F11)) * (1.0 - cos(2.0 * Omega_ijkl)) + V3 * (1.0 + cos(3.0 * Omega_ijkl))));
		Econj += F12 * Pcot_1 * (1.0 - sqr(sin_ijkl) * sin_ijk * sin_jkl);
	      }
	    }
	  }
	}
      }
    }
  }
  return(Etor + Econj);
}
//*****************************************************************************************************************************

double calc_Potential(int Natoms, int *NBond, int *Bonded, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, atom *particle, Bond_order *BO, point L, double *delta, double *delta_lp, Bond_order *BO_tag, double *Bmat, double *qvec, double *Avec, double *rijShiled_1, FILE* logfile, double *rij, AngleParamsStruct *AngleParams, double *delta_boc, double *n_lp, double *angle, int *Once_Bond, int *Nonce, double *f11, dihedralParamsStruct *dihedralParams, double *Exp_Ptor_2, int *BondNum,double *Exp_Pcot_2,double *delta_tag_boc,double *delta_e)
{
  
  int atomi;
  double Ebond, E_over, Elp, Ecoul, Eunder, EvdW;
  double Epot, Ekin, Eval,Etors_Econj;
  
  //Ebond  = calc_Ebond(Natoms,NBond,Bonded,particle,BondParams,BO);
  Ebond  = 0.0;
  //Elp    = calc_E_lp(Natoms,AtomParams,NBond,delta_lp,particle);
  Elp    = 0.0;;
  //E_over = calc_Eover(particle,BondParams,BO_tag,Natoms,NBond,Bonded,BO,delta,delta_lp,AtomParams);
  E_over = 0.0;
  //Eunder = calc_Eunder(Natoms,particle,AtomParams,delta,BO,NBond,Bonded,delta_lp);
  Eunder = 0.0;
  EvdW   = calc_EvdW(Natoms,rij,particle,BondParams);
  //EvdW   = 0.0;
  Ecoul  = calc_Ecoulomb(Natoms,particle,rij,Bmat,qvec,Avec,AtomParams,BondParams,rijShiled_1,logfile);
  //Ecoul  = 0.0;
  //Eval   = calc_Eval(Natoms,particle,NBond,Bonded,delta_boc,BO,angle,AngleParams,n_lp,AtomParams,delta_e);
  Eval   = 0.0;
  //Etors_Econj = calc_Etor_conj(particle,delta_boc,BO,NBond,Bonded,BondNum,angle,Natoms,dihedralParams,Exp_Pcot_2,Exp_Ptor_2,f11,L);
  Etors_Econj = 0.0;
  
  Epot   = Ebond + Elp + E_over + Eunder + EvdW + Ecoul + Eval + Etors_Econj;
  
  //**
  //FILE *Test;
  //Test = fopen("Test.dat","w");
  //fprintf(Test,"%.16f %.16f %.16f %.16f %.16f %.16f %.16f/n",Ebond,E_over,EvdW,Ecoul,Etors_Econj,Epot);
  //**
  /*
    cerr<<"Ebond ="<<Ebond<<endl;
    cerr<<"Elp   ="<<Elp<<endl;
    cerr<<"E_over ="<<E_over<<endl;
    cerr<<"Eunder ="<<Eunder<<endl;
    cerr<<"EvdW ="<<EvdW<<endl;
    cerr<<"Ecoul ="<<Ecoul<<endl;
    cerr<<"Eval ="<<Eval<<endl;
    cerr<<"Etors_Econj ="<<Etors_Econj<<endl;
    cerr<<"Epot ="<<Epot<<endl;
  */

  //return(Epot);
  return(Epot);
}
