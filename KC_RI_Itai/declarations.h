#ifndef declerations_h
#define declerations_h

// ********************** Libraries **********************

//#include <iostream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <omp.h>
#include <sys/timeb.h>
#include <cstring>
using std::cerr;
using std::endl;
using std::cout;
//using namespace std;
// ********************** constants  **********************

#define LCC          (1.42)                                // CC bond length in Angstroms.
#define rC           (0.5*LCC)
#define PIE          (3.14159265358979323846)
#define TWOPI        (6.28318530717958647692)
#define aut          (2.418884326505e-17)                  // 1 a.u. = 2.418884326505e-17 sec
#define auv          (2.1876912633e6)                      // 1 a.u. = 2.1876912633e6 m/sec       
#define aul          (0.529177e-10)                        // 1 a.u. = 0.529177e-10 m
#define auE          (5.14220642e+11)                      // 1 a.u. = 5.14220642e+11 Volts/m
#define auV          (27.2113845)                          // 1 a.u. = 27.2113845 Volts
#define auI          (6.62361782e-3)                       // 1 a.u. = 6.62361782e-3 A
#define sqrttwo      (1.414213562373095145474621858738828)
#define e            (1.602177e-19)                        // Coulomb
#define h            ( (6.62608e-34) / (1.602177e-19) )    // eV*Sec
#define Bohr         (0.529177)                            // 1 a.u.       = 0.529177 Angstrom
#define ev2au        (0.03675)                             // 0.03675 a.u. = 1 ev
//#define KB           (3.16681520371153e-6)                 // Boltzmann factor in au
//#define Tesla        ((Bohr)*(Bohr)*(1e-20)/(h))         // Tesla to a.u.
#define Tesla        (4.254381183e-6)                      // Tesla to a.u.
#define AKMA2SEC     (0.04888e-12)                         // 1 AKMA = 0.04888e-12 sec
#define SEC2AKMA     ((1.0) / (AKMA2SEC))                  // 1 / AKMA2SEC
#define kB           (0.001987191)                         // Boltzmann's constant in AKMA units (kcal/(mole K)) 0.00198635195
#define eV2kcalmole  (23.0609)                             // Multiply value in eV by this factor to get the corresponding value in kcal/mole.
#define sqr(x)       ((x)*(x))
#define square(x)       ((x)*(x))
#define cube(x)      ((x)*(x)*(x))
#define fourth(x)    ((x)*(x)*(x)*(x))

//#define MIN(A,B) ((A) < (B) ? (A) : (B))
//#define MAX(A,B) ((A) > (B) ? (A) : (B))

//Conjugate-gradients and Steepest-Descent definitions:

#define ITMAXCG            (1000000) // Maximun number of iterations for the conjugate gradient optimization
#define ITMAXSD            (1000000) // Maximum number of iterations for the steepest descent optimization.
#define ITMAXdbrent        (100)     // Maximum number of iterations for dbrent.
//Steepest-Descent definitions:

#define GOLD          (1.618034)
#define GLIMIT        (100.0)
#define TINY          (1e-20)
#define CGOLD         (0.3819660)
#define ZEPS          (1e-10)
#define SIGN(a,b)          ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MAX(a,b)           ((a) > (b) ? (a) : (b))
#define MIN(a,b)           ((a) > (b) ? (b) : (a))

//#define EPS           (1e-12) //1e-12!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define EPS           (1e-12) //1e-12!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define NAtomParams   (10)
//#define BOCutoff      (5.0) // AKMA (Angstrom) 3.25 for hydrocarbons  **need to adjust this parameter!!!!**
#define BOCutoff      (5.0) // AKMA (Angstrom) 3.25 for hydrocarbons  **need to adjust this parameter!!!!**
#define NeighbCutoff  (8.0) // AKMA (Angstrom) 9.75 for hydrocarbons **need to adjust this parameter!!!!**
#define Rcut_List     (8.0) // AKMA (Angstrom) 9.75 for hydrocarbons **need to adjust this parameter!!!!**
#define MaxNBond      (30)
#define MaxNeighb     (1200)
#define MaxList       (1800)
#define Normal_cuttoff (2.0)
//#define Frame_limit   (25) defines how much iteration a new traj.xyz file will open 
#define Frame_limit   (2000)
// *********************** Potential parameters  **********************

// Using the AKMA unit system (see for example: http://www.scripps.edu/rc/softwaredocs/msi/insight2K/charmm_principles/Ch01_intro.FM5.html)

// Lambdas have no units
#define Wd           (75.00) // Steepness of the Fermi damping functions used to smoothen SBO around 0.0 and 2.0 for the new definition of SBO2.
#define C1           (0.10) // Center of the left Fermi damping functions used to smoothen SBO around 0.0 and 2.0 for the new definition of SBO2.
#define C2           (1.90) // (2.0 - C1); Center of the right Fermi damping functions used to smoothen SBO around 0.0 and 2.0 for the new definition of SBO2.
#define lambda1     (50.00)   //p(boc1)
#define lambda2     (9.4569)  //p(boc2)
#define lambda3     (56.6636) //p(coa2)
#define lambda4     (3.0)     //p(trip4)
#define lambda5     (6.5)     //p(trip3)
#define lambda6     (50.0)    //kc2
#define lambda7     (1.0701)  //p(ov/on6)
#define lambda8     (15.0)    //p(trip2)
#define lambda9     (11.9083) //p(ov/on7)
#define lambda10    (13.3822) //p(ov/on8)
#define lambda11    (-24.6710)//p(trip1)
#define swa         (0.0)     //lower Tapper-radius
#define swb         (7.0)    //upper Tapper-radius
#define lambda14    (2.8793)  //not used
#define Pval_6      (33.8667) //(lambda15)(Pval_7)
#define lambda16    (5.8971)  //P(lp1)
#define Pval_9      (1.0563)  //(lambda17)
#define Pval_10     (2.0384)  //(lambda18)
#define lambda19    (6.1431)  //not used
#define lambda20    (6.9290)  //P(pen2)
#define lambda21    (0.3989)  //P(pen3))
#define lambda22    (3.9954)  //P(pen4)
#define lambda23    (-2.4837) //not used
#define Ptor_2      (5.8374)  //(lambda24)
#define Ptor_3      (10.0)    //(lambda25)
#define Ptor_4      (1.8820)  //(lambda26)
#define lambda27    (-1.2327) //not used
#define Pcot_2      (2.1861)  //(lambda28)

#define Wd3          (130.0) // Steepness of the Fermi damping functions used to smoothen f12 close to BO=0. 
#define C3           (0.25)  // center of the Fermi damping functions used to smoothen f12 close to BO=0.

#define lambda29    (1.5591)  //P(vdw1)
#define lambda30    (0.01)    //cuttoff for bond order (*100)
#define lambda31    (0.7151)  //P(coa4)
#define Povun4      (2.7425)  //P(ov/un4) was 2.7425
#define Povun3      (12.5819) //P(ov/un3) was 12.5819
#define Pval_8      (2.1533)  //(lambda34)
#define lambda35    (0.5)     //not used
#define lambda36    (20.0)    //not used
#define lambda37    (5.0)     //Molecular Energy
#define lambda38    (0.0)     //Molecular Energy
#define lambda39    (1.4155)  //P(coa3)

//#define non_bond_cut (swb-swa)
#define non_bond_cut (7.0)
//#define non_bond_cut (50.0)
//#define BO_cut       (0.0001)//this cutoff is responsible for energy discontinuities 
#define BO_cut       (0.0001)
//#define thb_cut      (0.001)//valance angle cutoff, responsible for energy discontinuities
//#define thb_cut      (0.001)//valance angle cutoff, responsible for energy discontinuities 
#define thb_cut      (0.001)//valance angle cutoff

//defining Tapper for long range interactions.

#define swa2 (sqr( swa ))
#define swa3 (cube(swa ))
#define swb2 (sqr( swb ))
#define swb3 (cube(swb ))
#define Tap_7 (20.0 / pow( (swb - swa), 7.0 ))
#define Tap_6 (-70.0 * (swa + swb) / pow( (swb - swa), 7.0 ))
#define Tap_5 (84.0 * (swa2 + 3.0*swa*swb + swb2) / pow( (swb - swa), 7.0 ))
#define Tap_4 (-35.0 * (swa3 + 9.0*swa2*swb + 9.0*swa*swb2 + swb3 ) / pow( (swb - swa), 7.0 ))
#define Tap_3 (140.0 * (swa3*swb + 3.0*swa2*swb2 + swa*swb3 ) / pow( (swb - swa), 7.0 ))
#define Tap_2 (-210.0 * (swa3*swb2 + swa2*swb3) /pow( (swb - swa), 7.0 ) )
#define Tap_1 (140.0 * swa3 * swb3 / pow( (swb - swa), 7.0 ))
#define Tap_0 ((-35.0*swa3*swb2*swb2 + 21.0*swa2*swb3*swb2 + 7.0*swa*swb3*swb3 + swb3*swb3*swb ) / pow( (swb - swa), 7.0 ))


// C6 coeficients in kcal/mol

//********Coefficients for the TS vdW calculation scheme*********

#define  To_Kcal_mol_Angstrom6 (627.503*pow(Bohr,6.0))
#define  alpha0N   7.4   // Free atom dipole polarizability a.u.
#define  alpha0C   7.4   // Free atom dipole polarizability a..u
#define  alpha0B  21.0   // Free atom dipole polarizability a.u.
#define  alpha0H   4.5   // Free atom dipole polarizability a.u.
#define  C60N     29.04  // Free atom C6 coefficient in Hartree*Bohr^6
#define  C60B     119.4 // Free atom C6 coefficient in Hartree*Bohr^6
#define  C60H     7.8   // Free atom C6 coefficient in Hartree*Bohr^6
#define  C60C     7.8   // Free atom C6 coefficient in Hartree*Bohr^6
#define  VeffN    0.86   // Relative Hirshfeld volume for N in h-BN
#define  VeffB    0.768  // Relative Hirshfeld volume for B in h-BN
#define  VeffH    0.62   // Relative Hirshfeld volume for B in h-BN
#define  VeffC    0.62   // Relative Hirshfeld volume for B in h-BN
#define  alphaN   (VeffN * alpha0N) // Bulk dipole polarizability a.u.
#define  alphaB   (VeffB * alpha0B) // Bulk dipole polarizability a.u.
#define  alphaH   (VeffH * alpha0H) // Bulk dipole polarizability a.u.
#define  alphaC   (VeffC * alpha0C) // Bulk dipole polarizability a.u.
#define  NN_C6    (sqr(VeffN) * C60N * To_Kcal_mol_Angstrom6) // Bulk C6 coefficient in Hartree*Bohr^6
#define  BB_C6    (sqr(VeffB) * C60B * To_Kcal_mol_Angstrom6) // Bulk C6 coefficient in Hartree*Bohr^6
#define  HH_C6    (sqr(VeffH) * C60H * To_Kcal_mol_Angstrom6) // Bulk C6 coefficient in Hartree*Bohr^6
#define  CC_C6    (sqr(VeffC) * C60C * To_Kcal_mol_Angstrom6) // Bulk C6 coefficient in Hartree*Bohr^6
#define  BN_C6     ( 2.0*NN_C6 * BB_C6 /(BB_C6*alpha0N/alpha0B + NN_C6*alpha0B/alpha0N) ) // Bulk C6 coefficient in Hartree*Bohr^6
#define  BH_C6     ( 2.0*BB_C6 * HH_C6 /(BB_C6*alpha0H/alpha0B + HH_C6*alpha0B/alpha0H) ) // Bulk C6 coefficient in Hartree*Bohr^6
#define  NH_C6     ( 2.0*NN_C6 * HH_C6 /(NN_C6*alpha0H/alpha0N + HH_C6*alpha0N/alpha0H) ) // Bulk C6 coefficient in Hartree*Bohr^6

#define  CN_C6     ( 2.0*CC_C6 * NN_C6 /(NN_C6*alpha0C/alpha0N + CC_C6*alpha0N/alpha0C) ) // Bulk C6 coefficient in Hartree*Bohr^6
#define  CB_C6     ( 2.0*CC_C6 * BB_C6 /(BB_C6*alpha0C/alpha0B + CC_C6*alpha0B/alpha0C) ) // Bulk C6 coefficient in Hartree*Bohr^6
#define  CH_C6     ( 2.0*CC_C6 * HH_C6 /(HH_C6*alpha0C/alpha0H + CC_C6*alpha0H/alpha0C) ) // Bulk C6 coefficient in Hartree*Bohr^6

#define  r0N       3.34 
#define  r0B       3.89
#define  r0H       3.1
#define  r0C       3.1 
#define  r_eff_NN (2.0*pow(VeffN,1.0/3.0)*r0N*Bohr)
#define  r_eff_BB (2.0*pow(VeffB,1.0/3.0)*r0B*Bohr)
#define  r_eff_HH (2.0*pow(VeffH,1.0/3.0)*r0H*Bohr)
#define  r_eff_CC (2.0*pow(VeffC,1.0/3.0)*r0C*Bohr)
#define  r_eff_BN (0.5*r_eff_BB + 0.5*r_eff_NN)
#define  r_eff_BH (0.5*r_eff_BB + 0.5*r_eff_HH)
#define  r_eff_NH (0.5*r_eff_NN + 0.5*r_eff_HH)

#define  r_eff_CN (0.5*r_eff_CC + 0.5*r_eff_NN)
#define  r_eff_CB (0.5*r_eff_CC + 0.5*r_eff_BB)
#define  r_eff_CH (0.5*r_eff_CC + 0.5*r_eff_HH)

#define  C_vdW    (0.068) //definins the degree of the tranlational component in vdW 
// C6 coeficients in kcal/mol
#define  d_TS     (15.0) // value of 20.0 in the original paper
//#define  d_TS     (20.0) // value of 20.0 in the original paper
#define  Sr_TS    (0.84) 
//**********************************************************************

//#define NN_rvdW_long (3.69) //Old value
#define NN_rvdW_long (3.8)
#define BN_rvdW_long (2.60)
#define BB_rvdW_long (3.1)
#define BH_rvdW_long (2.8)
#define NH_rvdW_long (2.7)
#define HH_rvdW_long (2.7)

#define CH_rvdW_long (2.8)
#define CB_rvdW_long (2.7)
#define CN_rvdW_long (2.7)

// rvdws have units of Angstroms which are also AKMA units of length.

#define HH_rvdw   (2.0*1.3525) //is the 2* correct??????
#define BB_rvdw   (2.0*1.65) //is the 2* correct??????
#define NN_rvdw   (2.0*1.7695) //is the 2* correct??????
#define NH_rvdw   (2.0*1.647)  //is the 2* correct??????
#define BH_rvdw   (2.0*1.501)  //is the 2* correct??????
#define BN_rvdw   (2.0*1.7)    //is the 2* correct??????

#define CH_rvdw   (2.0*1.647)  //is the 2* correct??????
#define CN_rvdw   (2.0*1.501)  //is the 2* correct??????
#define CB_rvdw   (2.0*1.7)    //is the 2* correct??????

//long range rvdW

// Epsilons have units of Kcal/mole which are also AKMA units of energy.
#define H_Epsilon (0.0616) //Dij in some articles
#define B_Epsilon (0.05)   //Dij in some articlse
#define N_Epsilon (0.1375) //Dij in some article
#define BH_Epsilon (0.0566) //Dij in some articles
#define NH_Epsilon (0.0567) //Dij in some articles
#define BN_Epsilon (0.0564) //Dij in some articles

#define CH_Epsilon (0.0566) //Dij in some articles
#define CB_Epsilon (0.0567) //Dij in some articles
#define CN_Epsilon (0.0564) //Dij in some articles

// Epsilons have units of Kcal/mole which are also AKMA units of energy.
#define H_Epsilon_long (0.31) //Dij in some articles
#define B_Epsilon_long (0.46)   //Dij in some articlse
#define N_Epsilon_long (0.21) //Dij in some article
#define BH_Epsilon_long (0.31) //Dij in some articles
#define NH_Epsilon_long (0.25) //Dij in some articles
#define BN_Epsilon_long (0.2) //Dij in some articles

#define CH_Epsilon_long (0.31) //Dij in some articles
#define CN_Epsilon_long (0.25) //Dij in some articles
#define CB_Epsilon_long (0.2) //Dij in some articles

//vdw inner wall parameters:

#define H_rcore  (0.6)
#define B_rcore  (1.4)
#define N_rcore  (1.4)
#define H_ecore  (0.1)
#define B_ecore  (0.1)
#define N_ecore  (0.1)
#define H_acore  (10.0)
#define B_acore  (12.0)
#define N_acore  (10.0)

#define C_rcore  (10.0)
#define C_ecore  (12.0)
#define C_acore  (10.0)

// alphas have no units.

#define H_alpha  (9.3858)
#define B_alpha  (12.4662)
#define N_alpha  (10.0667)
#define BH_alpha (11.2019)
#define NH_alpha (10.5106)
#define BN_alpha (10.7561)

#define CH_alpha (11.2019)
#define CN_alpha (10.5106)
#define CB_alpha (10.7561)

//alpha long range

#define H_alpha_long  (9.0)
#define B_alpha_long  (11.0)
#define N_alpha_long  (12.08)
#define BH_alpha_long (9.0)
#define NH_alpha_long (9.0)
#define BN_alpha_long (10.7)

#define CH_alpha_long (9.0)
#define CN_alpha_long (9.0)
#define CB_alpha_long (10.7)

//Trans
#define H_Trans  (20.0)
//#define B_Trans  (1.8)
#define B_Trans  (0.8) //Old value
//#define N_Trans  (1.20)
#define N_Trans  (0.8) //Old value
#define BH_Trans (20.0)
#define NH_Trans (20.0)
#define BN_Trans (1.8)


#define CH_Trans (20.0)
#define CN_Trans (20.0)
#define CB_Trans (1.8)

// lambsaWs have units of 1/Angstroms which are 1/AKMA units of length.

#define H_lambdaW (5.0013)
#define B_lambdaW (2.6721)
#define N_lambdaW (7.688)
#define C_lambdaW (7.6886)

//*************original parameters********************

#define original_H_gamma (0.8)
#define original_B_gamma (0.7)
#define original_N_gamma (0.69)
#define original_C_gamma (0.69)


// Etas have units of eV - this is used for the calculation of the partial charges of the atoms. For the Coulomb energy we convert to Kcal/mole.
#define H_Etha  (7.0327)
#define B_Etha  (6.7020)
#define N_Etha  (7.0)
#define C_Etha  (7.0)

// Kais have units of eV - this is used for the calculation of the partial charges of the atoms. For the Coulomb energy we convert to Kcal/mole.
#define H_Kai  (10.2)
#define N_Kai  (13.0)
#define B_Kai  (10.0)
#define C_Kai  (10.0)

// Kais have units of 1/Angstrom - this is used for the calculation of the partial charges of the atoms.
#define H_gamma (0.8)
#define B_gamma (0.7)
#define N_gamma (0.69)
#define C_gamma (0.69)

#define Kappa   (14.4)
#define Kappa_eV_kcal (Kappa*eV2kcalmole)
// ***********************Tersoff Params*********************
#define BN_Tersoff_c (25000)
#define BN_Tersoff_d (4.3484)
#define BN_Tersoff_sqr_c (sqr(BN_Tersoff_c))
#define BN_Tersoff_sqr_d (sqr(BN_Tersoff_d))
#define BN_Tersoff_h (-0.89)
#define BN_Tersoff_n (0.72751)
#define BN_Tersoff_A (1380*eV2kcalmole)
#define BN_Tersoff_B (340.0*eV2kcalmole)
#define BN_Tersoff_beta (1.25724e-7)
#define BN_Tersoff_lambda1 (3.568)
#define BN_Tersoff_lambda2 (2.199)
#define BN_Tersoff_lambda3 (0.0)
#define BN_Tersoff_D (0.05)
#define BN_Tersoff_R (1.95)

#define NN_Tersoff_c (17.7959)
#define NN_Tersoff_d (5.9484)
#define NN_Tersoff_sqr_c (sqr(NN_Tersoff_c))
#define NN_Tersoff_sqr_d (sqr(NN_Tersoff_d))
#define NN_Tersoff_h (0.0)
#define NN_Tersoff_n (0.6184432)
#define NN_Tersoff_A (128.868*eV2kcalmole)
#define NN_Tersoff_B (138.778*eV2kcalmole)
#define NN_Tersoff_beta (0.019251)
#define NN_Tersoff_lambda1 (2.8293093)
#define NN_Tersoff_lambda2 (2.6272721)
#define NN_Tersoff_lambda3 (0.0)
#define NN_Tersoff_D (0.1)
#define NN_Tersoff_R (2.0)

#define BB_Tersoff_c (0.52629)
#define BB_Tersoff_d (0.001587)
#define BB_Tersoff_sqr_c (sqr(BB_Tersoff_c))
#define BB_Tersoff_sqr_d (sqr(BB_Tersoff_d))
#define BB_Tersoff_h (0.5)
#define BB_Tersoff_n (3.9929061)
#define BB_Tersoff_A (40.520156*eV2kcalmole)
#define BB_Tersoff_B (43.132016*eV2kcalmole)
#define BB_Tersoff_beta (1.6e-6)
#define BB_Tersoff_lambda1 (2.2372578)
#define BB_Tersoff_lambda2 (2.0774982)
#define BB_Tersoff_lambda3 (0.0)
#define BB_Tersoff_D (0.1)
#define BB_Tersoff_R (2.0)

#define CC_Tersoff_c (38049)
#define CC_Tersoff_d (4.3484)
#define CC_Tersoff_sqr_c (sqr(CC_Tersoff_c))
#define CC_Tersoff_sqr_d (sqr(CC_Tersoff_d))
#define CC_Tersoff_h (-0.93)
#define CC_Tersoff_n (0.72751)
#define CC_Tersoff_A (1393.6*eV2kcalmole)
#define CC_Tersoff_B (430.0*eV2kcalmole)
#define CC_Tersoff_beta (1.5724e-7)
#define CC_Tersoff_lambda1 (3.4879)
#define CC_Tersoff_lambda2 (2.2119)
#define CC_Tersoff_lambda3 (0.0)
#define CC_Tersoff_D (0.15)
#define CC_Tersoff_R (1.95)
// *********************** structures  **********************

typedef struct
{
  int type,layer,lock;
  double charge, dcharge,Xprev, Xnext, Yprev, Ynext, Zprev, Znext;
  double r[3];
  double Mass, Fx, Fy, Fz, Ek;
  double ax, ay, az, vx, vy, vz;
  double RI,Inter;
  double conjx, conjy, conjz;
  double gradUx,gradUy,gradUz,newgradUx,newgradUy,newgradUz;
} atom;
typedef struct
{
  double x, y, z,Normal_length,vector_1x,vector_1y,vector_1z,vector_2x,vector_2y,vector_2z,not_norm_x,not_norm_y,not_norm_z;
} Normal_struct;
typedef struct
{
  double x, y, z;
} point;
typedef struct
{
  int x, y, z;
} Counter;
typedef struct
{
  double r0_Sigma, r0_Pi, r0_Pi_Pi,C6,rvdW_long;
  double Pover, lambdaW, alpha, rvdw, Epsilon,Epsilon_long;
  double alpha_long,Trans,r_eff;
  double original_gamma,gamma,pow_gamma,Povun1,rcore,ecore,acore;
  double Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,Tersoff_A,Tersoff_B,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_beta;
  double Tersoff_R,Tersoff_D,Tersoff_sqr_d,Tersoff_sqr_c;
} BondParamsStruct;
typedef struct
{
  double Povun5, Povun2, Plp2,Val_boc;
  double Etha, kai;
  double Pval_3,Pval_5;
} AtomParamsStruct;

typedef struct
{
  double r[3];
} dr_ij;

/**************************************************************/
/*                     3D dot product                         */

inline double dot_prod(point p1, point p2)
{
  return(p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
}

/**************************************************************/
/*                    3D cross product                        */

inline point cross_prod(point p1, point p2)
{
  point result;

  result.x = p1.y * p2.z - p1.z * p2.y;
  result.y = p1.z * p2.x - p1.x * p2.z;
  result.z = p1.x * p2.y - p1.y * p2.x;

  return(result);
}

/**************************************************************/
/*                        sign(x)                             */

inline double sign(double x)
{
  double sgn;
  
  if(x > 0)       sgn =  1.0;
  else if (x < 0) sgn = -1.0;
  else            sgn =  0.0;
  
  return(sgn);
}
// **********************inline****************************

inline double Sp(double Xij, double Xmin, double Xmax, double &dX)
{
  double cutoff;
  
  double t = (Xij-Xmin) / (Xmax-Xmin);
  if (t <= 0.0) {
    cutoff = 1.0;
    dX = 0.0;
  } else if (t >= 1.0) {
    cutoff = 0.0;
    dX = 0.0;
  } else {
    cutoff = 0.5 * (1.0+cos(t*PIE));
    dX = (-0.5*PIE*sin(t*PIE)) / (Xmax-Xmin);
  }
  return cutoff;
};

  /* ----------------------------------------------------------------------
     LJ cutoff function Sp2
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */
 static inline double powint(const double &x, const int n) {
    double yy,ww;

    if (x == 0.0) return 0.0;
    int nn = (n > 0) ? n : -n;
    ww = x;

    for (yy = 1.0; nn != 0; nn >>= 1, ww *=ww)
      if (nn & 1) yy *= ww;

    return (n > 0) ? yy : 1.0/yy;
  }

inline double Sp2(double Xij, double Xmin, double Xmax, double &dX)
{
  double cutoff;
  
  double t = (Xij-Xmin) / (Xmax-Xmin);
  if (t <= 0.0) {
    cutoff = 1.0;
    dX = 0.0;
  } else if (t >= 1.0) {
    cutoff = 0.0;
    dX = 0.0;
  } else {
    cutoff = (1.0-(t*t*(3.0-2.0*t)));
    dX = 6.0*(t*t-t) / (Xmax-Xmin);
  }
  return cutoff;
};
/* kronecker delta function returning a double */

inline double kronecker(const int a, const int b)
{
  return (a == b) ? 1.0 : 0.0;
};
// ********************** Functions  **********************
//double calc_dP_ij_dr(int atomi,int atomj, int atomk,int dir,atom *particle,int Natoms,double *rij,double r_ij,Normal_struct *Normal,point L);
//void calc_F_vdW(point *FvdW, int Natoms, atom *particle, double *rij, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom)
void read_file(char *filename);
void set_Normal_atom(double *rij,int Natoms,Normal_struct *Normal,atom *particle,int *Normal_atom,point L);
double Calc_Pij(int,int,atom *particle,Normal_struct *Normal,double r_ij,point L);
void calc_dNormal_k(int &atomi,int &atomn,int &dir,dr_ij *dN,atom *particle,Normal_struct *Normal,int *Normal_atom,point &L);
//void calc_dNormal_k(int atomi,int atomn,int dir,dr_ij *dN,atom *particle,Normal_struct *Normal,int *Normal_atom,point L);
double calc_EvdW(int Natoms, double *rijShiled_1, atom *particle, BondParamsStruct *BondParams,Normal_struct *Normal);
void calc_Normal(int Natoms,int *Normal_atom,atom *particle,Normal_struct *Normal,point L);

int GetNatoms(char *filename);
void ReadCoor(char *filename, int Natoms, atom* particle,int Interlayer,point &L);
void read_RUN_parameters(double &tini, double &tend, double &dt, point &L, double &T, int &NConserve, int &NTraj, double &OptTol, int &CalcType,int &Interlayer,int &Periodic_flag,int &nthreads,int & Coulomb_flag);
void Init_Bond_Params(BondParamsStruct *BondParams);
void Init_Atom_Params(AtomParamsStruct *AtomParams);
void calc_angles(double *angle, int Natoms, int *NBond, int *Bonded, double *rij, atom *particle,point L);
double calc_dthetaijk_drn(int atomi, int atomj, int atomk, int atomn, int dir, atom *particle, double *angle, int *BondNum, int Natoms, double *rij,double Thetaijk);
void calc_rij(int *NBond, int *Bonded, int *BondNum, int Natoms, atom *particle, point L, double *rij,double *Fc,BondParamsStruct *BondParams,int *Neighb,int *NNeighb,double *rijShiled_1,int *NList,int *List,int *NBond_once,int *Bonded_once);
void Calc_Charge(atom *particle, int Natoms, double *Bmat, double *qvec, double *Avec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rij, double *rijShiled_1, FILE* logfile,point L);
double calc_Ecoulomb(int Natoms, atom *particle, double *rij, double *Bmat, double *qvec, double *Avec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rijShiled_1, FILE* logfile,int Coulomb_flag,point L);
void Invert(int N, double *A, double *InvA, FILE *logfile);
void MatVec(double *A, double *VecIn, int N, double *VecOut);
//double calc_EvdW(int Natoms, double *rij, atom *particle, BondParamsStruct *BondParams);
double calc_dTap(double r_ij);
double calc_Tap(double r_ij);
//Force
void TORSION_(atom *particle,point *f,int *BondNum,int Natoms,point L,int *Bonded,int *NBond,double &E_TORSION,int *NBond_once,int *Bonded_once);
void relax_cell(int Natoms,int *NBond,AtomParamsStruct *AtomParams,BondParamsStruct *BondParams,int *Bonded,double *Bmat,double *qvec,double *Avec,double *rij,double *rijShiled_1,point &L, atom *particle,Normal_struct *Normal,double *Fc,FILE* logfile,double *angle,int *BondNum,int *Neighb,int *NNeighb,int Coulomb_flag,int Interlayer,int *Normal_atom,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once);
void Print_Conserve(FILE *EFile, int Natoms, int *NBond, int *Bonded, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, atom *particle, point L, double t, double *Bmat, double *qvec, double *Avec, double *rijShiled_1, FILE* logfile, time_t &t0, clock_t &c0, double *rij, double *angle, int *BondNum,Normal_struct *Normal,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double Epot,int *NBond_once,int *Bonded_once);
void Calc_Force(point *Fterm, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double &E_REBO,int *NBond_once,int *Bonded_once);
double calc_dr_ij_drk(int atomi, int atomj, int atomk, int dir, double r_ij, atom *particle, int Natoms,point L);
void calc_force_arrays(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *angle, double *Fc,int *Neighb,int *NNeighb,double *rijShiled_1,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once);
void Propagate(int Natoms, int *NBond, int *Bonded, atom *particle, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, int *BondNum, double *rij, point &L, double dt, double tini, double tfin, double T, int NConserve, int NTraj, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile, time_t &t0, clock_t &c0,double *angle,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,point *R0,int *NList,int *List,int *NBond_once,int *Bonded_once);
void Print_Traj(FILE *TrajFile,FILE *TrajFile_RI,FILE *TrajFile_Inter, int Natoms, double t, atom *particle,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,point L);
void Print_Traj_gv(int Natoms, atom *particle,point L);
void calc_dcharge(int atomk, int dir, atom *particle, int Natoms, double *Bmat, double *qvec, double *Avec, double *dqvec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rij, double *rijShiled_1, FILE* logfile,point L);
double calc_Potential(int Natoms, int *NBond, int *Bonded, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, atom *particle, double *Bmat, double *qvec, double *Avec, double *rijShiled_1, FILE* logfile, double *rij,double *angle, int *BondNum,Normal_struct *Normal,double *Fc,int *Neighb,int *NNeighb,point L,int Coulomb_flag);
void steepest_descent(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams,double *angle, int NTraj, int NConserve, double OptTol, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int Periodic_flag,int *Neighb,int *NNeighb,int Coulomb_flag,int *NList,int *List,point *R0);
void conj_grad(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams,double *angle, double OptTol, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,int NConserve,int NTraj,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int Periodic_flag,int *Neighb,int *NNeighb,int Coulomb_flag,int *NList,int *List,point *R0);
void Quench(int Natoms, int *NBond, int *Bonded, atom *particle, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams,int *BondNum,double *rij, point &L, double dt, double tini, double tfin, double T, int NConserve, int NTraj, double *Bmat, double *qvec, double *dqvec, double *Avec, double *rijShiled_1, FILE* logfile, time_t &t0, clock_t &c0, double *angle, double OptTol,Normal_struct *Normal,int *Normal_atom,int Interlayer, double *Fc,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,int *Neighb,int *NNeighb,int Periodic_flag,int Coulomb_flag,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once);
double R_ij(atom *particle,int atomi,int atomj,point L);
double Fc_(double r_ij,double Tersoff_R,double Tersoff_D);
void UpdateNeighbList(int Dim, atom *particle, int *List, int *NList, point *R0, point L);
void Check_List(atom *particle, point *R0,int Natoms, int *List, int *NList, point L);
void vel(atom *particle, double T, int Natoms, long *idum);
double ran_nrc(long *idum);
void Randomize();
//double ran();
double randv(long *);
point vel_center(atom *particle, double MTot, int Natoms);
double calc_val_angle(int atomi,int atomj, int atomk, int Natoms, int *NBond, int *Bonded,double r_jk, double r_ij, atom *particle,point L,dr_ij &rji,dr_ij &rjk,double inside_acos);
void calc_F_Rep(point *FvdW, int Natoms, atom *particle, double *rijShiled_1, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,dr_ij *dN,double &ERep);
double bondorder(int, int, double *, double, double, point *, int, atom *particle,int *Bonded,int *BondNum,int *NBond,point L);
double R_PBC(atom *particle,int atomi,int atomj,int dir,double L);
double calc_dE_REBO_pot(double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_REBO,int *NBond_once,int *Bonded_once);
double Calc_Potential(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *rijShiled_1, FILE* logfile,double *angle,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double &E_tot,int *NBond_once,int *Bonded_once,int,int,FILE*);
double bondorder_pot(int i, int j, double rij[3], double rijmag, double VA, int vflag_atom,atom *particle,int *Bonded,int *BondNum,int *NBond,point L);
double TORSION_pot(atom *particle,int *BondNum,int Natoms,point L,int *Bonded,int *NBond,double &E_TORSION,int *NBond_once,int *Bonded_once);
double calc_F_Rep_pot(int Natoms, atom *particle, double *rijShiled_1, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,dr_ij *dN,double &ERep_,double &ERepulsion,double &EDisp,double &RI,double &Interlayer_dis);
double CircleOverlap(double d, double r1, double r2);
#endif
