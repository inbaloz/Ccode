#include "declarations.h"
void Init_Atom_Params(AtomParamsStruct *AtomParams)
{
  int i;
  
  for(i=0 ; i < NAtomParams ; i++){
    AtomParams[i].kai     = 0.0;
    AtomParams[i].Etha    = 0.0;
  }
  
  // H
  
  AtomParams[1].kai     = H_Kai;
  AtomParams[1].Etha    = H_Etha;
  
  // B
  
  AtomParams[5].kai     = B_Kai;
  AtomParams[5].Etha    = B_Etha;
  
  // N
  
  AtomParams[7].kai     = N_Kai;
  AtomParams[7].Etha    = N_Etha;
}

void Init_Bond_Params(BondParamsStruct *BondParams)
{
  int i, j;
  int index1, index2;
  
  for(i=0 ; i <  NAtomParams; i++){ 
    for(j=0 ; j < NAtomParams ; j++){
      index1 = i*NAtomParams + j;
     
      BondParams[index1].lambdaW  = 0.0;
      BondParams[index1].alpha    = 0.0;
      BondParams[index1].rvdw     = 0.0;
      BondParams[index1].Epsilon  = 0.0;
      BondParams[index1].gamma    = 0.0;
      BondParams[index1].original_gamma = 0.0;
      BondParams[index1].rcore    = 0.0;
      BondParams[index1].ecore    = 0.0;
      BondParams[index1].acore    = 0.0;
      BondParams[index1].C6       = 0.0;
      BondParams[index1].rvdW_long  = 0.0;
      BondParams[index1].Epsilon_long  = 0.0;
      BondParams[index1].alpha_long = 0.0;
      BondParams[index1].r_eff      = 0.0;
      BondParams[index1].Tersoff_c    = 0.0;
      BondParams[index1].Tersoff_d    = 0.0;
      BondParams[index1].Tersoff_h    = 0.0;
      BondParams[index1].Tersoff_n    = 0.0;
      BondParams[index1].Tersoff_A    = 0.0;
      BondParams[index1].Tersoff_B    = 0.0;
      BondParams[index1].Tersoff_R    = 0.0;
      BondParams[index1].Tersoff_D    = 0.0;
      BondParams[index1].Tersoff_lambda1 = 0.0;
      BondParams[index1].Tersoff_lambda2 = 0.0;
      BondParams[index1].Tersoff_lambda3 = 0.0;
      BondParams[index1].Tersoff_beta    = 0.0;
      BondParams[index1].Tersoff_sqr_c   = 0.0;
      BondParams[index1].Tersoff_sqr_d   = 0.0;      
    }
  }
  
  // H-H
  
  index1 = 1*NAtomParams + 1;

  BondParams[index1].lambdaW  = H_lambdaW;
  BondParams[index1].rvdw     = HH_rvdw;
  BondParams[index1].Epsilon  = H_Epsilon;
  BondParams[index1].Epsilon_long  = H_Epsilon_long;
  BondParams[index1].alpha_long    = H_alpha_long;
  BondParams[index1].rvdW_long     = HH_rvdW_long;
  BondParams[index1].gamma    = H_gamma;
  BondParams[index1].Trans    = H_Trans;
  BondParams[index1].original_gamma    = original_H_gamma;
  BondParams[index1].rcore    = H_rcore;
  BondParams[index1].ecore    = H_ecore;
  BondParams[index1].acore    = H_acore;
  BondParams[index1].C6       = HH_C6;
  BondParams[index1].r_eff    = r_eff_HH;

  //B-B
  
  index1 = 5*NAtomParams + 5;
  
  BondParams[index1].lambdaW  = B_lambdaW;
  BondParams[index1].alpha    = B_alpha;
  BondParams[index1].rvdw     = BB_rvdw;
  BondParams[index1].Epsilon  = B_Epsilon;
  BondParams[index1].Epsilon_long  = B_Epsilon_long;
  BondParams[index1].alpha_long = B_alpha_long;
  BondParams[index1].rvdW_long      = BB_rvdW_long;
  BondParams[index1].gamma    = B_gamma;
  BondParams[index1].pow_gamma  = pow(1.0/B_gamma,3.0);
  BondParams[index1].Trans    = B_Trans;
  BondParams[index1].original_gamma    = original_B_gamma;
  BondParams[index1].rcore    = B_rcore;
  BondParams[index1].ecore    = B_ecore;
  BondParams[index1].acore    = B_acore;
  BondParams[index1].C6       = BB_C6;
  BondParams[index1].r_eff    = r_eff_BB;
  
  BondParams[index1].Tersoff_c       = BB_Tersoff_c;
  BondParams[index1].Tersoff_d       = BB_Tersoff_d;
  BondParams[index1].Tersoff_sqr_c   = BB_Tersoff_sqr_c;
  BondParams[index1].Tersoff_sqr_d   = BB_Tersoff_sqr_d;
  BondParams[index1].Tersoff_h       = BB_Tersoff_h;
  BondParams[index1].Tersoff_n       = BB_Tersoff_n;
  BondParams[index1].Tersoff_A       = BB_Tersoff_A;
  BondParams[index1].Tersoff_B       = BB_Tersoff_B;
  BondParams[index1].Tersoff_lambda1 = BB_Tersoff_lambda1;
  BondParams[index1].Tersoff_lambda2 = BB_Tersoff_lambda2;
  BondParams[index1].Tersoff_lambda3 = BB_Tersoff_lambda3;
  BondParams[index1].Tersoff_beta    = BB_Tersoff_beta;
  BondParams[index1].Tersoff_R       = BB_Tersoff_R;
  BondParams[index1].Tersoff_D       = BB_Tersoff_D;

  //C-C
  
  index1 = 6*NAtomParams + 6;
  
  BondParams[index1].Tersoff_c       = CC_Tersoff_c;
  BondParams[index1].Tersoff_d       = CC_Tersoff_d;
  BondParams[index1].Tersoff_sqr_c   = CC_Tersoff_sqr_c;
  BondParams[index1].Tersoff_sqr_d   = CC_Tersoff_sqr_d;
  BondParams[index1].Tersoff_h       = CC_Tersoff_h;
  BondParams[index1].Tersoff_n       = CC_Tersoff_n;
  BondParams[index1].Tersoff_A       = CC_Tersoff_A;
  BondParams[index1].Tersoff_B       = CC_Tersoff_B;
  BondParams[index1].Tersoff_lambda1 = CC_Tersoff_lambda1;
  BondParams[index1].Tersoff_lambda2 = CC_Tersoff_lambda2;
  BondParams[index1].Tersoff_lambda3 = CC_Tersoff_lambda3;
  BondParams[index1].Tersoff_beta    = CC_Tersoff_beta;
  BondParams[index1].Tersoff_R       = CC_Tersoff_R;
  BondParams[index1].Tersoff_D       = CC_Tersoff_D;
  //B-H
  
  index1 = 1*NAtomParams + 5;
  index2 = 5*NAtomParams + 1;


  BondParams[index1].lambdaW  = BondParams[index2].lambdaW  = sqrt(B_lambdaW*H_lambdaW);
  BondParams[index1].alpha    = BondParams[index2].alpha    = BH_alpha;
  BondParams[index1].rvdw     = BondParams[index2].rvdw     = BH_rvdw;
  BondParams[index1].Epsilon  = BondParams[index2].Epsilon  = BH_Epsilon;
  BondParams[index1].Epsilon_long  = BondParams[index2].Epsilon_long  = BH_Epsilon_long;
  BondParams[index1].alpha_long    = BondParams[index2].alpha_long    = BH_alpha_long;
  BondParams[index1].rvdW_long      = BondParams[index2].rvdW_long      = BH_rvdW_long;
  BondParams[index1].gamma    = BondParams[index2].gamma    = sqrt(B_gamma*H_gamma);
  BondParams[index1].original_gamma    = BondParams[index2].original_gamma    = sqrt(original_B_gamma*original_H_gamma);
  BondParams[index1].rcore    = BondParams[index2].rcore    = sqrt(B_rcore*H_rcore);
  BondParams[index1].ecore    = BondParams[index2].ecore    = sqrt(B_ecore*H_ecore);
  BondParams[index1].acore    = BondParams[index2].acore    = sqrt(B_acore*H_acore);
  BondParams[index1].Trans    = BondParams[index2].Trans    = BH_Trans;
  BondParams[index1].r_eff    = BondParams[index2].r_eff    = r_eff_BH;
      
  //B-N
  
  index1 = 5*NAtomParams + 7;
  index2 = 7*NAtomParams + 5;
  
  BondParams[index1].lambdaW  = BondParams[index2].lambdaW  = sqrt(B_lambdaW*N_lambdaW);
  BondParams[index1].alpha    = BondParams[index2].alpha    = BN_alpha;
  BondParams[index1].rvdw     = BondParams[index2].rvdw     = BN_rvdw;
  BondParams[index1].Epsilon  = BondParams[index2].Epsilon  = BN_Epsilon;
  BondParams[index1].Epsilon_long  = BondParams[index2].Epsilon_long  = BN_Epsilon_long;
  BondParams[index1].alpha_long    = BondParams[index2].alpha_long    = BN_alpha_long;
  BondParams[index1].rvdW_long      = BondParams[index2].rvdW_long      = BN_rvdW_long;
  BondParams[index1].gamma     = BondParams[index2].gamma    = sqrt(B_gamma*N_gamma);
  BondParams[index1].pow_gamma = BondParams[index2].pow_gamma  = pow(1.0/sqrt(B_gamma*N_gamma),3.0);
  BondParams[index1].original_gamma    = BondParams[index2].original_gamma    = sqrt(original_B_gamma*original_N_gamma);
  BondParams[index1].rcore    = BondParams[index2].rcore    = sqrt(B_rcore*N_rcore);
  BondParams[index1].ecore    = BondParams[index2].ecore    = sqrt(B_ecore*N_ecore);
  BondParams[index1].acore    = BondParams[index2].acore    = sqrt(B_acore*N_acore);
  BondParams[index1].C6       = BondParams[index2].C6       = BN_C6;
  BondParams[index1].Trans    = BondParams[index2].Trans    = BN_Trans;
  BondParams[index1].r_eff    = BondParams[index2].r_eff    = r_eff_BN;
  
  BondParams[index1].Tersoff_c       =  BondParams[index2].Tersoff_c       = BN_Tersoff_c;
  BondParams[index1].Tersoff_d       =  BondParams[index2].Tersoff_d       = BN_Tersoff_d;
  BondParams[index1].Tersoff_sqr_c   =  BondParams[index2].Tersoff_sqr_c   = BN_Tersoff_sqr_c; 
  BondParams[index1].Tersoff_sqr_d   =  BondParams[index2].Tersoff_sqr_d   = BN_Tersoff_sqr_d; 
  BondParams[index1].Tersoff_h       =  BondParams[index2].Tersoff_h       = BN_Tersoff_h;
  BondParams[index1].Tersoff_n       =  BondParams[index2].Tersoff_n       = BN_Tersoff_n;
  BondParams[index1].Tersoff_A       =  BondParams[index2].Tersoff_A       = BN_Tersoff_A;
  BondParams[index1].Tersoff_B       =  BondParams[index2].Tersoff_B       = BN_Tersoff_B;
  BondParams[index1].Tersoff_lambda1 =  BondParams[index2].Tersoff_lambda1 = BN_Tersoff_lambda1;
  BondParams[index1].Tersoff_lambda2 =  BondParams[index2].Tersoff_lambda2 = BN_Tersoff_lambda2;
  BondParams[index1].Tersoff_lambda3 =  BondParams[index2].Tersoff_lambda3 = BN_Tersoff_lambda3;
  BondParams[index1].Tersoff_beta    =  BondParams[index2].Tersoff_beta    = BN_Tersoff_beta;
  BondParams[index1].Tersoff_R       =  BondParams[index2].Tersoff_R       = BN_Tersoff_R;
  BondParams[index1].Tersoff_D       =  BondParams[index2].Tersoff_D       = BN_Tersoff_D;

  
  //N-N
  
  index1 = 7*NAtomParams + 7;
  
  BondParams[index1].lambdaW  = N_lambdaW;
  BondParams[index1].alpha    = N_alpha;
  BondParams[index1].rvdw     = NN_rvdw;
  BondParams[index1].Epsilon  = N_Epsilon;
  BondParams[index1].Epsilon_long  = N_Epsilon_long;
  BondParams[index1].alpha_long    = N_alpha_long;
  BondParams[index1].rvdW_long      = NN_rvdW_long;
  BondParams[index1].gamma    = N_gamma;
  BondParams[index1].pow_gamma  = pow(1.0/N_gamma,3.0);
  BondParams[index1].original_gamma = original_N_gamma;
  BondParams[index1].rcore    = N_rcore;
  BondParams[index1].ecore    = N_ecore;
  BondParams[index1].acore    = N_acore;
  BondParams[index1].C6       = NN_C6;
  BondParams[index1].Trans    = N_Trans;
  BondParams[index1].r_eff    = r_eff_NN;

  BondParams[index1].Tersoff_c       = NN_Tersoff_c;
  BondParams[index1].Tersoff_d       = NN_Tersoff_d;
  BondParams[index1].Tersoff_sqr_c   = NN_Tersoff_sqr_c;
  BondParams[index1].Tersoff_sqr_d   = NN_Tersoff_sqr_d;
  BondParams[index1].Tersoff_h       = NN_Tersoff_h;
  BondParams[index1].Tersoff_n       = NN_Tersoff_n;
  BondParams[index1].Tersoff_A       = NN_Tersoff_A;
  BondParams[index1].Tersoff_B       = NN_Tersoff_B;
  BondParams[index1].Tersoff_lambda1 = NN_Tersoff_lambda1;
  BondParams[index1].Tersoff_lambda2 = NN_Tersoff_lambda2;
  BondParams[index1].Tersoff_lambda3 = NN_Tersoff_lambda3;
  BondParams[index1].Tersoff_beta    = NN_Tersoff_beta;
  BondParams[index1].Tersoff_R       = NN_Tersoff_R;
  BondParams[index1].Tersoff_D       = NN_Tersoff_D;
  //N-H
  
  index1 = 1*NAtomParams + 7;
  index2 = 7*NAtomParams + 1;
  
  BondParams[index1].lambdaW  = BondParams[index2].lambdaW  = sqrt(N_lambdaW*H_lambdaW);
  BondParams[index1].alpha    = BondParams[index2].alpha    = NH_alpha;
  BondParams[index1].rvdw     = BondParams[index2].rvdw     = NH_rvdw;
  BondParams[index1].Epsilon  = BondParams[index2].Epsilon  = NH_Epsilon;
  BondParams[index1].Epsilon_long  = BondParams[index2].Epsilon_long  = NH_Epsilon_long;
  BondParams[index1].alpha_long   = BondParams[index2].alpha_long    = NH_alpha_long;
  BondParams[index1].rvdW_long      = BondParams[index2].rvdW_long      = NH_rvdW_long;
  BondParams[index1].gamma    = BondParams[index2].gamma    = sqrt(H_gamma*N_gamma);
  BondParams[index1].original_gamma    = BondParams[index2].original_gamma    = sqrt(original_N_gamma*original_H_gamma);
  BondParams[index1].rcore    = BondParams[index2].rcore    = sqrt(H_rcore*N_rcore);
  BondParams[index1].ecore    = BondParams[index2].ecore    = sqrt(H_ecore*N_ecore);
  BondParams[index1].acore    = BondParams[index2].acore    = sqrt(H_acore*N_acore);
  BondParams[index1].C6       = BondParams[index2].C6       = NH_C6;
  BondParams[index1].Trans    = BondParams[index2].Trans    = NH_Trans;
  BondParams[index1].r_eff    = BondParams[index2].r_eff    = r_eff_NH;
}

