#include "declarations.h"

/********************************************************************/
/* Routines required for the calculation of the potential energy.   */
/* Written by Itai Leven and Oded Hod June 2010.                    */
/* (1)ReaxFF - J. Phys. Chem. A 105, 9396-9409 (2001)- hydrocarbons */
/* (2)ReaxFF - J. Phys. Chem. A 107, 3803-3811 (2003)- SiO          */
/* (3)ReaxFF - J. Phys. Chem. A 123, 114703    (2005)- B-N-H        */ 
/********************************************************************/

/**************************************************/
/*                    E_{bond}                    */
/**************************************************/

/**********************************************************************/
/* Calculate dr_{ij}_dr_k                                             */

double calc_dr_ij_drk(int atomi, int atomj, int atomk, int dir, double r_ij, atom *particle, int Natoms, point L)
{
  double Length;

  //r_ij  = R_ij(particle,atomi,atomj, L);
  
  if (dir == 0)    Length = L.x;
  else if (dir ==1)Length = L.y;
  else Length = L.z;
  
  if(atomk == atomi)      return( (particle[atomi].r[dir] - particle[atomj].r[dir] - (Length * rint( (particle[atomi].r[dir] - particle[atomj].r[dir]) / (Length + TINY)))) / r_ij);
  else if(atomk == atomj) return( (particle[atomj].r[dir] - particle[atomi].r[dir] - (Length * rint( (particle[atomj].r[dir] - particle[atomi].r[dir]) / (Length + TINY)))) / r_ij);
  else                    return(0.0);
}


/**************************************************/
/*                   F_{vdW}                      */
/**************************************************/

double calc_dTap(double r_ij)
{
  double dTap;
  
  dTap = 7*Tap_7 * r_ij + 6*Tap_6;
  dTap = dTap * r_ij + 5*Tap_5;
  dTap = dTap * r_ij + 4*Tap_4;
  dTap = dTap * r_ij + 3*Tap_3;
  dTap = dTap * r_ij + 2*Tap_2;
  dTap = dTap * r_ij + Tap_1;
    
  return(dTap);
}

/****************************************************************/
double calc_dP_ij_dr(int atomi,int atomj, int atomk,int dir,atom *particle,int Natoms,double *rij,double r_ij,point *Normal,point L,int *Normal_atom,int Interlayer)
{
  
  double inside_acos,rx,ry,rz,dR_ij;
  double angle,dPij,dinside_acos,dangle;

  dr_ij *dN = new dr_ij[1];
  
  dN[0].r[0] =0.0;
  dN[0].r[1] =0.0;
  dN[0].r[2] =0.0;
  
  dR_ij = calc_dr_ij_drk(atomi,atomj,atomk,dir,r_ij,particle,Natoms,L);
  
  rx = particle[atomi].r[0] - particle[atomj].r[0];
  ry = particle[atomi].r[1] - particle[atomj].r[1];
  rz = particle[atomi].r[2] - particle[atomj].r[2];
  
  rx -= L.x * rint(rx / (L.x + TINY));
  ry -= L.y * rint(ry / (L.y + TINY));
  rz -= L.z * rint(rz / (L.z + TINY));
  
  inside_acos = (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z) / (r_ij);
  
  if (Interlayer == 1)
    calc_dNormal_k(atomi,atomk,dir,dN,particle,Normal,Normal_atom,L);
  
  if (atomk == atomi){
    if (dir == 0)dinside_acos = (Normal[atomi].x + rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/(r_ij) - (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z)*(pow(r_ij,-2.0))*dR_ij;
    if (dir == 1)dinside_acos = (Normal[atomi].y + rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/(r_ij) - (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z)*(pow(r_ij,-2.0))*dR_ij;
    if (dir == 2)dinside_acos = (Normal[atomi].z + rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/(r_ij) - (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z)*(pow(r_ij,-2.0))*dR_ij;
  }
  else if (atomk == atomj){
    if (dir == 0)dinside_acos = (-Normal[atomi].x + rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/(r_ij) - (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z)*(pow(r_ij,-2.0))*dR_ij;
    if (dir == 1)dinside_acos = (-Normal[atomi].y + rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/(r_ij) - (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z)*(pow(r_ij,-2.0))*dR_ij;
    if (dir == 2)dinside_acos = (-Normal[atomi].z + rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/(r_ij) - (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z)*(pow(r_ij,-2.0))*dR_ij;
  }
  else{
    dinside_acos = (rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/r_ij - (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z)*(pow(r_ij,-2.0))*dR_ij;
  }
  
  if (fabs(inside_acos) < 1.0){
    angle  = acos(inside_acos);
    dangle = -1.0/sqrt(1.0-sqr(inside_acos))*dinside_acos;
  }
  else if (inside_acos >= 1.0) {
    angle  = 0.0;
    dangle = 0.0;
  }
  else {
    angle  = PIE;
    dangle = 0.0;
  }
  
  dPij = dR_ij*sin(angle) + r_ij*cos(angle)*dangle;
  
  delete [] dN;
  
  return(dPij);
  //return(dN[0].r[0]+dN[0].r[1]+dN[0].r[2]);
}

//*********************************

void calc_F_vdW(point *FvdW, int Natoms, atom *particle, double *rijShiled_1, BondParamsStruct *BondParams, point L,point *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb)
{
  /*
  int atomk, atomi, typek, typei,atomj,typej, index;
  double lambdaW, alpha, rvdW, Epsilon, f13, df13_drij, term1, term2, rij_lambda29, dEvdW_drij,dInner_shield;
  double ecore, acore, rcore, dTap, r_ij, Tap,Inner_shield,C6;
  double dP_ij_drk,dP_ji_drk,R,Pij,Pji,gamma,dr_ij_drk;
  double exp_Pij_gamma, exp_Pji_gamma,Pij_gamma,Pji_gamma,reff,alpha_dr_ij_drk;
  int n,dir;
  */
  //dP_ij_drk = 0.0;
  //dP_ji_drk = 0.0;
  
#pragma omp parallel
  {
    int atomk, atomi, typek, typei,atomj,typej, index,j;
    double lambdaW, alpha, rvdW, Epsilon, f13, df13_drij, term1, term2, rij_lambda29, dEvdW_drij,dInner_shield;
    double ecore, acore, rcore, dTap, r_ij, Tap,Inner_shield,C6,Epsilon_C_vdW_exp,term3,term4,term5,term6,term7;
    double dP_ij_drk,dP_ji_drk,R,Pij,Pji,gamma,dr_ij_drk,exp_term1,_term2,C6_r_ij;
    double exp_Pij_gamma, exp_Pji_gamma,Pij_gamma,Pji_gamma,reff,alpha_dr_ij_drk,r_ik,r_jk;
    int n,dir;
    dP_ij_drk = 0.0;
    dP_ji_drk = 0.0;
#pragma omp for firstprivate(Interlayer) 
    for(atomk=0 ; atomk < Natoms ; atomk++){
      
      typek = particle[atomk].type;
      //atomj = atomk;
      //typej = typek;
      FvdW[atomk].x = FvdW[atomk].y = FvdW[atomk].z = 0.0;
      
      for(atomi=0 ; atomi < Natoms ; atomi++){
	typei = particle[atomi].type;
	r_ik  = R_ij(particle,atomi,atomk, L);
	for (j=0 ; j < NNeighb[atomi]; j++){
	  atomj = Neighb[atomi*MaxNeighb + j] ;	  
	  typej = particle[atomj].type;
	  
	  //r_jk = rijShiled_1[atomi*MaxNeighb + j];
	  //if ((r_jk < Normal_cuttoff) || (r_ik < Normal_cuttoff)){
	  
	  if (atomk==atomj || atomk == atomi || atomk ==Normal_atom[atomi*3 + 0] || atomk ==Normal_atom[atomi*3 + 1] || atomk ==Normal_atom[atomi*3 + 2] || atomk ==Normal_atom[atomj*3 + 0] || atomk ==Normal_atom[atomj*3 + 1] || atomk ==Normal_atom[atomj*3 + 2]){
	    //r_jk  = R_ij(particle,atomk,atomj, L);
	    index = typej*NAtomParams + typei;
	    //r_ij  = R_ij(particle,atomi,atomj, L);
	    r_ij = rijShiled_1[atomi*MaxNeighb + j];
	    Tap   = calc_Tap(r_ij);
	    dTap  = calc_dTap(r_ij);  
	    
	    //interlayer vdW
	    
	    reff    = BondParams[index].r_eff;	    
	    alpha   = BondParams[index].alpha_long;
	    Epsilon = BondParams[index].Epsilon_long;
	    R       = BondParams[index].rvdW_long;
	    C6      = BondParams[index].C6;
	    gamma   = BondParams[index].Trans;
	    
	    Pij   = Calc_Pij(atomi,atomj,particle,Normal,r_ij,L);
	    Pji   = Calc_Pij(atomj,atomi,particle,Normal,r_ij,L);
	    Pij_gamma = Pij/gamma; //new
	    Pji_gamma = Pji/gamma; //new
	    exp_Pij_gamma = exp(-pow((Pij_gamma),2.0)); //new
	    exp_Pji_gamma = exp(-pow((Pji_gamma),2.0)); //new
	    exp_term1 = exp(alpha*(1.0-r_ij/R));
	    //**
	    //term1 = alpha*(1.0-r_ij/R);
	    //**
	    term2   = exp(-d_TS*(r_ij/(Sr_TS*reff) - 1.0));
	    _term2  = 1.0/(1.0 + term2);
	    Epsilon_C_vdW_exp = Epsilon + C_vdW*(exp_Pij_gamma + exp_Pji_gamma);
	    C6_r_ij = C6/(pow(r_ij,6.0)+0.1);
	    term3   = exp_term1*(-alpha/R)*(Epsilon_C_vdW_exp);
	    term4   = exp_term1*C_vdW*(exp_Pij_gamma*(-2.0)*Pij_gamma/gamma);
	    term5   = exp_term1*C_vdW*(exp_Pji_gamma*(-2.0)*Pji_gamma/gamma);
	    term6   = (-(pow(_term2,2.0)*term2*(-d_TS/(Sr_TS*reff))*C6_r_ij ) + _term2*(-1)*C6*pow(pow(r_ij,6)+0.1,-2)*6*pow(r_ij,5));
	    term7   = ((exp_term1*(Epsilon_C_vdW_exp))-(_term2)*C6_r_ij)*dTap;
	    //********dirivative in the x direction******
	    
	    dir = 0;
	    dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	    dP_ij_drk = calc_dP_ij_dr(atomi,atomj,atomk,dir,particle, Natoms,rijShiled_1,r_ij,Normal,L,Normal_atom,Interlayer);
	    dP_ji_drk = calc_dP_ij_dr(atomj,atomi,atomk,dir,particle, Natoms,rijShiled_1,r_ij,Normal,L,Normal_atom,Interlayer);
	    
	    //alpha_dr_ij_drk = alpha*(dr_ij_drk/R);
	    
	    dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk - (term6*dr_ij_drk)) + dr_ij_drk*term7;
	    
	    FvdW[atomk].x -= dEvdW_drij;
	    
	    //********dirivative in the y direction******
	    
	    dir = 1;
	    dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	    dP_ij_drk = calc_dP_ij_dr(atomi,atomj,atomk,dir,particle, Natoms,rijShiled_1,r_ij,Normal,L,Normal_atom,Interlayer);
	    dP_ji_drk = calc_dP_ij_dr(atomj,atomi,atomk,dir,particle, Natoms,rijShiled_1,r_ij,Normal,L,Normal_atom,Interlayer);
	    
	    //alpha_dr_ij_drk = alpha*(dr_ij_drk/R);
	    
	    dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk - (term6*dr_ij_drk)) + dr_ij_drk*term7;
	    
	    FvdW[atomk].y -= dEvdW_drij;
	    
	    //********dirivative in the z direction******
	    
	    dir = 2;
	    dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	    dP_ij_drk = calc_dP_ij_dr(atomi,atomj,atomk,dir,particle, Natoms,rijShiled_1,r_ij,Normal,L,Normal_atom,Interlayer);
	    dP_ji_drk = calc_dP_ij_dr(atomj,atomi,atomk,dir,particle, Natoms,rijShiled_1,r_ij,Normal,L,Normal_atom,Interlayer);
	    
	    //alpha_dr_ij_drk = alpha*(dr_ij_drk/R);
	    
	    dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk - (term6*dr_ij_drk)) + dr_ij_drk*term7;
	    
	    FvdW[atomk].z -= dEvdW_drij;
	  }
	}
      }
    }
  }

}
//***********************************************
void calc_F_coul(point *Fcoul, atom *particle, int Natoms, double *Bmat, double *qvec, double *Avec, double *dqvec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rij, double *rijShiled_1, FILE* logfile, point L,int *Neighb,int *NNeighb,int Coulomb_flag)
{
  int atomi, atomj, atomk, index,j;
  double term1,term2;
  double Tap, dTap, r_ij;
  double rijShiled,shield,gamma,pow_gamma;
  double drij_dr_k_x,drij_dr_k_y,drij_dr_k_z;
  double dTerm;

  for(atomk=0 ; atomk < Natoms ; atomk++)
    Fcoul[atomk].x = Fcoul[atomk].y = Fcoul[atomk].z = 0.0;
  //**
  if(Coulomb_flag == 1){
    
  }
  //**
#pragma omp parallel for private(dqvec,Avec,atomi, atomj, atomk, index,j,term1,term2,Tap, dTap, r_ij,rijShiled,shield,gamma,drij_dr_k_x,drij_dr_k_y,drij_dr_k_z)     
  for(atomi=0 ; atomi < Natoms ; atomi++){
    for (j=0 ; j < NNeighb[atomi]; j++){
      
      atomj = Neighb[atomi*MaxNeighb + j] ;
      
      //if(atomk != atomj && atomk !=atomi && Coulomb_flag == 0)continue;
      
      //index = atomi*Natoms + atomj;
      
      //r_ij  = R_ij(particle,atomi,atomj, L);
      r_ij = rijShiled_1[atomi*MaxNeighb + j];
      Tap  = calc_Tap(r_ij);
      dTap = calc_dTap(r_ij);
      gamma = BondParams[particle[atomi].type*NAtomParams + particle[atomj].type].original_gamma;
      pow_gamma = BondParams[particle[atomi].type*NAtomParams + particle[atomj].type].original_gamma;
      shield    = pow(r_ij,3.0) + pow(1.0/gamma,3.0);
      //shield    = pow(r_ij,3.0) + pow_gamma;
      rijShiled = 1.0/pow(shield,0.333333333333);
      term1=- qvec[atomi] * qvec[atomj] * fourth(rijShiled) * sqr(r_ij);
      term2= qvec[atomi] * qvec[atomj] * rijShiled;
      
      drij_dr_k_x = calc_dr_ij_drk(atomi,atomj,atomi,0,r_ij,particle,Natoms,L);
      drij_dr_k_y = calc_dr_ij_drk(atomi,atomj,atomi,1,r_ij,particle,Natoms,L);
      drij_dr_k_z = calc_dr_ij_drk(atomi,atomj,atomi,2,r_ij,particle,Natoms,L);
      
      dTerm = (Tap * term1 + dTap * term2)*Kappa_eV_kcal;
      Fcoul[atomi].x -= drij_dr_k_x*dTerm;
      Fcoul[atomi].y -= drij_dr_k_y*dTerm;
      Fcoul[atomi].z -= drij_dr_k_z*dTerm;
      
      Fcoul[atomj].x += drij_dr_k_x*dTerm;
      Fcoul[atomj].y += drij_dr_k_y*dTerm;
      Fcoul[atomj].z += drij_dr_k_z*dTerm;
      /*
	if(Coulomb_flag == 0){
	Fcoul[atomk].x -= drij_dr_k_x*(Tap * term1 + dTap * term2);
	Fcoul[atomk].y -= drij_dr_k_y*(Tap * term1 + dTap * term2);
	Fcoul[atomk].z -= drij_dr_k_z*(Tap * term1 + dTap * term2);
	}
	else{
	
	}
      */
    }
    //Fcoul[atomi].x *= (Kappa * eV2kcalmole);
    //Fcoul[atomi].y *= (Kappa * eV2kcalmole);
    //Fcoul[atomi].z *= (Kappa * eV2kcalmole);
  }
  
  
}

//delete [] dcharge_x;
//delete [] dcharge_y;
  //delete [] dcharge_z;

/****************************************************************/
/*  dThetaijk_drn                                               */

// Calculating the derivative with respect to movments of atom-n of
// the valence angle {ijk}.

double calc_dthetaijk_drn(int atomi, int atomj, int atomk, int atomn, int dir, atom *particle, double *angle, int *BondNum, int Natoms, double *rij, point L,double Thetaijk)
{
  int indexji, indexjk;
  double ri, rj, rk, rjk_1, rji_1, cosijk;
  double Length;
  
  indexji = atomj*Natoms + atomi;
  indexjk = atomj*Natoms + atomk;
  
  //Thetaijk = angle[atomj * MaxNBond * MaxNBond + BondNum[indexji] * MaxNBond + BondNum[indexjk]];
  
  if((fabs(Thetaijk) > EPS) && (fabs(Thetaijk - PIE) > EPS)){
    if(atomn == atomi){
      ri = particle[atomi].r[dir];
      rj = particle[atomj].r[dir];
      rk = particle[atomk].r[dir];
      
      rjk_1 = 1.0 /R_ij(particle,atomk,atomj, L);
      rji_1 = 1.0 /R_ij(particle,atomi,atomj, L);
      
      if (dir == 0)    Length = L.x;
      else if (dir ==1)Length = L.y;
      else Length = L.z;
      
      cosijk = cos(Thetaijk);
      return( (1.0 / sqrt(1.0 - sqr(cosijk))) * rji_1 * ((ri - rj - (Length * rint( (ri - rj) / (Length + TINY)))) * rji_1 * cosijk - (rk - rj - (Length * rint( (rk - rj) / (Length + TINY)))) * rjk_1) );
    }
    else if(atomn == atomk){
      ri = particle[atomi].r[dir];
      rj = particle[atomj].r[dir];
      rk = particle[atomk].r[dir];
      
      rjk_1 = 1.0 / R_ij(particle,atomk,atomj, L);
      rji_1 = 1.0 / R_ij(particle,atomi,atomj, L);
      
      cosijk = cos(Thetaijk);
      
      
      if (dir == 0)    Length = L.x;
      else if (dir ==1)Length = L.y;
      else Length = L.z;
      return( (1.0 / sqrt(1.0 - sqr(cosijk))) * rjk_1 * ((rk - rj - (Length * rint( (rk - rj) / (Length + TINY)))) * rjk_1 * cosijk - (ri - rj - (Length * rint( (ri - rj) / (Length + TINY)))) * rji_1) );
    }
    else if(atomn == atomj){
      ri = particle[atomi].r[dir];
      rj = particle[atomj].r[dir];
      rk = particle[atomk].r[dir];
      
      rjk_1 = 1.0 / R_ij(particle,atomk,atomj, L);
      rji_1 = 1.0 / R_ij(particle,atomi,atomj, L);
      
      cosijk = cos(Thetaijk);
      
      if (dir == 0)    Length = L.x;
      else if (dir ==1)Length = L.y;
      else Length = L.z;
      //Length * rint( (ri - rj) / (Length + TINY))
      return( - (1.0 / sqrt(1.0 - sqr(cosijk))) * (rji_1 * ((ri - rj - (Length * rint( (ri - rj) / (Length + TINY)))) * rji_1 * cosijk - (rk - rj - (Length * rint( (rk - rj) / (Length + TINY)))) * rjk_1) + rjk_1 * ((rk - rj - (Length * rint( (rk - rj) / (Length + TINY)))) * rjk_1 * cosijk - (ri - rj - (Length * rint( (ri - rj) / (Length + TINY)))) * rji_1)) );
    }
    else return(0.0);
  }
  else return(0.0);
}


double dFc_drn(int atomi,int atomj,int atomn,double r_ij,double Tersoff_R,double Tersoff_D,double dr_ij_drn){
  
  double dFc_ij_drn;
  
  if(r_ij < (Tersoff_R-Tersoff_D) || r_ij > (Tersoff_R+Tersoff_D)) dFc_ij_drn=0.0;
  else dFc_ij_drn = -0.5*cos(PIE/2*(r_ij-Tersoff_R)/Tersoff_D)*PIE*0.5/Tersoff_D*dr_ij_drn;
  
  return (dFc_ij_drn);
}

double calc_dE_Tersoff_ij_drn(int atomi,int atomj,int atomn,int dir,int i,double *angle,int *BondNum,int Natoms,double *rij,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond){
  
  int index,atomk,typek,typej,typei,k,index2;
  double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
  double Xsi,ETersoff,dTheta_ijk_drn,dg_drn,dXsi_drn,Fc_ij,Fc_jk,Tersoff_R,Tersoff_D;
  double dr_ij_drn,dr_jk_drn,dFc_ij_drn,dFc_jk_drn,Fa,Fr,dFa_drn,dFr_drn,db_ij_drn,n,dE_Tersoff_ij_drn;
  double Tersoff_sqr_c,Tersoff_sqr_d,sqr_Tersoff_h_Cos_Theta;
  //i=BondNum[atomj*Natoms + atomi];
  r_ij       = rij[atomj*MaxNBond + i];
  //r_ij       = R_ij(particle,atomi,atomj, L);
  typei      = particle[atomi].type;
  typej      = particle[atomj].type;
  index      = typei*NAtomParams + typej;
  
  Tersoff_n  = BondParams[index].Tersoff_n;
  A          = BondParams[index].Tersoff_A;
  B          = BondParams[index].Tersoff_B;
  beta       = BondParams[index].Tersoff_beta;
  Tersoff_R  = BondParams[index].Tersoff_R;
  Tersoff_D  = BondParams[index].Tersoff_D;
  Tersoff_lambda1 = BondParams[index].Tersoff_lambda1;
  Tersoff_lambda2 = BondParams[index].Tersoff_lambda2;
  Tersoff_lambda3 = BondParams[index].Tersoff_lambda3;
  Xsi=0;
  dXsi_drn=0;
  //Fc_ij = Fc_(r_ij,Tersoff_R,Tersoff_D);
  Fc_ij = Fc[atomj*MaxNBond + i];
  dr_ij_drn = calc_dr_ij_drk(atomi,atomj,atomn,dir,r_ij,particle,Natoms,L);
  dFc_ij_drn = dFc_drn(atomi,atomj,atomn,r_ij,Tersoff_R,Tersoff_D,dr_ij_drn);
  
  for(k=0 ; k < NBond[atomj] ; k++){ // Go over all j neighbors k != i.
    atomk = Bonded[atomj*MaxNBond + k];
    if(atomk!=atomi){
      
      typek      = particle[atomk].type;
      index2     = typek*NAtomParams + typej;
      Tersoff_c  = BondParams[index2].Tersoff_c;
      Tersoff_d  = BondParams[index2].Tersoff_d;
      Tersoff_sqr_c  = BondParams[index2].Tersoff_sqr_c;
      Tersoff_sqr_d  = BondParams[index2].Tersoff_sqr_d;
      
      Tersoff_h  = BondParams[index2].Tersoff_h;
      Tersoff_R  = BondParams[index2].Tersoff_R;
      Tersoff_D  = BondParams[index2].Tersoff_D;
      //r_jk  = R_ij(particle,atomk,atomj, L);
      r_jk      = rij[atomj*MaxNBond + k];
      //Fc_jk =  Fc_(r_jk,Tersoff_R,Tersoff_D);
      Fc_jk = Fc[atomj*MaxNBond + k];
      dr_jk_drn = calc_dr_ij_drk(atomj,atomk,atomn,dir,r_jk,particle,Natoms,L);
      Theta_ijk = angle[atomj * MaxNBond * MaxNBond + i * MaxNBond + k];
      sqr_Tersoff_h_Cos_Theta = sqr(Tersoff_h-cos(Theta_ijk));
      
      dTheta_ijk_drn = calc_dthetaijk_drn(atomi,atomj,atomk,atomn,dir,particle,angle,BondNum,Natoms,rij,L,Theta_ijk);
      g      = 1+Tersoff_sqr_c/Tersoff_sqr_d-Tersoff_sqr_c/(Tersoff_sqr_d+sqr_Tersoff_h_Cos_Theta);
      //g      = 1+sqr(Tersoff_c)/sqr(Tersoff_d)-sqr(Tersoff_c)/(sqr(Tersoff_d)+sqr(Tersoff_h-cos(Theta_ijk)));
      dg_drn = Tersoff_sqr_c*pow(Tersoff_sqr_d + sqr_Tersoff_h_Cos_Theta,-2)*2*(Tersoff_h-cos(Theta_ijk))*sin(Theta_ijk)*dTheta_ijk_drn;
      //dg_drn = sqr(Tersoff_c)*pow(sqr(Tersoff_d) + sqr(Tersoff_h-cos(Theta_ijk)),-2)*2*(Tersoff_h-cos(Theta_ijk))*(sin(Theta_ijk)*dTheta_ijk_drn);
      dFc_jk_drn = dFc_drn(atomj,atomk,atomn,r_jk,Tersoff_R,Tersoff_D,dr_jk_drn);
      Xsi  += Fc_jk*g;
      dXsi_drn +=dFc_jk_drn*g + Fc_jk*dg_drn;
    }
  }
  
  //**
  //Xsi=1;
  //dXsi_drn=0;
  //**
  
  Fa      = -B*exp(-Tersoff_lambda2*r_ij);
  Fr      = A*exp(-Tersoff_lambda1*r_ij);
  b_ij    = pow(1+pow(beta*Xsi,Tersoff_n),-1/(2*Tersoff_n));
  dFa_drn = -Fa*Tersoff_lambda2*dr_ij_drn;
  dFr_drn = -Fr*Tersoff_lambda1*dr_ij_drn;
  
  if(Xsi==0)db_ij_drn=0;
  else db_ij_drn = -pow(beta,Tersoff_n)*pow(Xsi,Tersoff_n-1)/(2*pow(1+pow(beta*Xsi,Tersoff_n),1+1/(2*Tersoff_n)))*dXsi_drn;
    
  dE_Tersoff_ij_drn = dFc_ij_drn*(Fr + b_ij*Fa) + Fc_ij*(dFr_drn + db_ij_drn*Fa + b_ij*dFa_drn);
    
  return(0.5*dE_Tersoff_ij_drn);
}

void calc_F_Tersoff(double *angle,int *BondNum,int Natoms,double *rij,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,point *FTersoff,int *NBond){
  
  int i, j,k, atomk, atomi, atomj;
  int hey=0;
  
#pragma omp parallel for private(i,j,k,atomi,atomj,atomk)          
  for(atomk=0 ; atomk < Natoms ; atomk++){
    
    FTersoff[atomk].x = FTersoff[atomk].y = FTersoff[atomk].z = 0.0;
    for(i=0 ; i < NBond[atomk] ; i++){
      atomi = Bonded[atomk*MaxNBond + i];
      atomj = atomk;
      FTersoff[atomk].x -= calc_dE_Tersoff_ij_drn(atomi,atomj,atomk,0,i,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond);
      FTersoff[atomk].y -= calc_dE_Tersoff_ij_drn(atomi,atomj,atomk,1,i,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond);
      FTersoff[atomk].z -= calc_dE_Tersoff_ij_drn(atomi,atomj,atomk,2,i,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond);
      
      for(j=0 ; j < NBond[atomi] ; j++){
	atomj=Bonded[atomi*MaxNBond + j];
	FTersoff[atomk].x -= calc_dE_Tersoff_ij_drn(atomj,atomi,atomk,0,j,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond);
	FTersoff[atomk].y -= calc_dE_Tersoff_ij_drn(atomj,atomi,atomk,1,j,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond);
	FTersoff[atomk].z -= calc_dE_Tersoff_ij_drn(atomj,atomi,atomk,2,j,angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,NBond);
	
      }
    }
  }
}

void calc_force_arrays(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *angle, double *Fc,int *Neighb,int *NNeighb,double *rijShiled_1,int *NList,int *List,point *R0)
{
  Check_List(particle,R0,Natoms,List,NList,L);
  calc_rij(NBond,Bonded,BondNum,Natoms,particle,L,rij,Fc,BondParams,Neighb,NNeighb,rijShiled_1,NList,List);
  calc_angles(angle,Natoms,NBond,Bonded,rij,particle,L);
}

/****************************************************************/
/* Calculate the total force array                              */

void Calc_Force(point *Fterm, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,point *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag)
{
  int atomi;
  
  // Calculate the bond force term
  
  calc_F_vdW(Fterm,Natoms,particle,rijShiled_1,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb);
  // Add F_vdW to the total force
#pragma omp parallel for
  for(atomi=0 ; atomi < Natoms ; atomi++){
    particle[atomi].Fx = Fterm[atomi].x;
    particle[atomi].Fy = Fterm[atomi].y;
    particle[atomi].Fz = Fterm[atomi].z;
  }
  
  
  /*
  calc_F_coul(Fterm,particle,Natoms,Bmat,qvec,Avec,dqvec,AtomParams,BondParams,rij,rijShiled_1,logfile,L,Neighb,NNeighb,Coulomb_flag);
  
  // Add F_Coulomb to the total force
#pragma omp parallel for
  for(atomi=0 ; atomi < Natoms ; atomi++){
    particle[atomi].Fx = Fterm[atomi].x;
    particle[atomi].Fy = Fterm[atomi].y;
    particle[atomi].Fz = Fterm[atomi].z;
  }
  
  
  calc_F_Tersoff(angle,BondNum,Natoms,rij,L,BondParams,Bonded,particle,Fc,Fterm,NBond);

#pragma omp parallel for
  for(atomi=0 ; atomi < Natoms ; atomi++){
    particle[atomi].Fx = Fterm[atomi].x;
    particle[atomi].Fy = Fterm[atomi].y;
    particle[atomi].Fz = Fterm[atomi].z;
  }
  */
}
