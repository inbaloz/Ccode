#include "declarations.h"

//void DSYSV(int N, int NRHS, double **A, double **B, double **Vout);

/***********************************************************************/
/* Calculate the charge of the different atoms in the molecule using   */ 
/* the EEM method as detailed in: Int. J. Mol. Sci. 8, 572-582 (2007). */
/* We use the screened Coulomb interaction as detailed in:             */
/* J. Mol. Catal. A: Chem. 134, 79-88 (1998)                           */
/***********************************************************************/

//The diagonal part of Bmat and Avec can be calculated once outside.

void Calc_Charge(atom *particle, int Natoms, double *Bmat, double *qvec, double *Avec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rij, double *rijShiled_1, FILE* logfile,point L)
{
  int atomi, atomj, dim, typei, typej, index;
  double Qtot; // Total charge of the molecule.
  double third, gammaij_1,rij_Shiled_1,r_ij;
  
  dim = Natoms + 1;
  Qtot = 0; // We assume that the molecule is neutral.
  third = (1.0 / 3.0);
  
  // Build Bmat of Eq. (5) in the reference above.
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    typei = particle[atomi].type;
    //rijShiled_1[atomi*Natoms + atomi] = -1.0;
    Bmat[atomi*dim + atomi] = 2.0 * AtomParams[typei].Etha;
    
    for(atomj=0 ; atomj < atomi ; atomj++){
      r_ij =R_ij(particle,atomi,atomj,L);
      index     = atomi*Natoms + atomj;
      typej     = particle[atomj].type;
      gammaij_1 = 1.0 / BondParams[typei*NAtomParams + typej].gamma;
      
      rij_Shiled_1 = 1.0 / pow( (cube(r_ij) + cube(gammaij_1))  , third );
      //rijShiled_1[atomj*Natoms + atomi] = rijShiled_1[index];
      
      Bmat[atomj*dim + atomi] = Bmat[atomi*dim + atomj] = Kappa * rij_Shiled_1; // in eV.
    }
  }
  
  
  for(atomi=0 ; atomi < Natoms ; atomi++) Bmat[atomi*dim + (dim-1)] = Bmat[(dim-1)*dim + atomi] = -1.0;
  
  Bmat[(dim-1)*dim + (dim-1)] = 0.0;
  
  // Build Avec of Eq. (5) in the reference above.
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    typei = particle[atomi].type;
    Avec[atomi] = -AtomParams[typei].kai;
  }
  
  Avec[dim-1] = -Qtot;
  
  // Solve for the atomic charges
  
  Invert(dim,Bmat,Bmat,logfile); // At this point Bmat contains Inv(B).
  MatVec(Bmat,Avec,dim,qvec);
  
  /*
  for (atomi=0;atomi < Natoms;atomi++)
    {
      if (particle[atomi].type == 1) qvec[atomi]=0.0;
      else if (particle[atomi].type == 5)qvec[atomi]=0.5;
      else if (particle[atomi].type == 7)qvec[atomi]=-0.5;
    }
  */
  //qvec[6]=0.257505;
  //qvec[8]=0.257505;
  //qvec[9]=0.257505;
  //qvec[7]=-0.063956;
  //qvec[10]=-0.063956;
  //qvec[11]=-0.063956;
  
  //DSYSV(dim,1,Bmat,Avec,qvec);
  
}

void calc_dcharge(int atomk, int dir, atom *particle, int Natoms, double *Bmat, double *qvec, double *Avec, double *dqvec, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, double *rij, double *rijShiled_1, FILE* logfile,point L)
{
  int atomi, atomj, index, dim,typei,typej;
  double drij_dr_k,r_ij,gammaij_1,rij_Shiled_1,third;
  
  //Calc_Charge(particle,Natoms,Bmat,qvec,Avec,AtomParams,BondParams,rij,rijShiled_1,logfile);
  // This needs to be calculated once for each atomic configuration.
  
  dim = Natoms + 1;
  third = (1.0 / 3.0);
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    Avec[atomi] = 0.0; // From this point on Avec holds the R.H.S. of the derivative matrix equation.
    typei = particle[atomi].type;
    for(atomj=0 ; atomj < atomi ; atomj++){
      typej     = particle[atomj].type;
      index = atomi*Natoms + atomj;
      r_ij=R_ij(particle,atomi,atomj,L);
      drij_dr_k = calc_dr_ij_drk(atomi,atomj,atomk,dir,r_ij,particle,Natoms,L);
      gammaij_1 = 1.0 / BondParams[typei*NAtomParams + typej].gamma;
      rij_Shiled_1 = 1.0 / pow( (cube(r_ij) + cube(gammaij_1))  , third );
      Avec[atomi] += fourth(rij_Shiled_1) * sqr(r_ij) * drij_dr_k * qvec[atomj];
    }
    
    for(atomj=atomi+1 ; atomj < Natoms ; atomj++){
      typej     = particle[atomj].type;
      index = atomi*Natoms + atomj;
      r_ij=R_ij(particle,atomi,atomj,L);
      drij_dr_k = calc_dr_ij_drk(atomi,atomj,atomk,dir,r_ij,particle,Natoms,L);
      gammaij_1 = 1.0 / BondParams[typei*NAtomParams + typej].gamma;
      rij_Shiled_1 = 1.0 / pow( (cube(r_ij) + cube(gammaij_1))  , third );
      Avec[atomi] += fourth(rij_Shiled_1) * sqr(r_ij) * drij_dr_k * qvec[atomj];
    }
    
    Avec[atomi] *= Kappa;
  }
  
  Avec[dim-1] = 0.0; // dQ_{tot}/dx_k = 0.
  
  MatVec(Bmat,Avec,dim,dqvec);
  /*
    for (atomi=0;atomi < Natoms;atomi++){
    dqvec[atomi]=0.0;
    }
  */
}
