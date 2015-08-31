
#include "declarations.h"

void Print_Conserve(FILE *EFile, int Natoms, int *NBond, int *Bonded, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, atom *particle, point L, double t, double *Bmat, double *qvec, double *Avec, double *rijShiled_1, FILE* logfile, time_t &t0, clock_t &c0, double *rij, double *angle, int *BondNum,Normal_struct *Normal,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double Epot,int *NBond_once,int *Bonded_once)
{
  int atomi;
  double Ekin=0.0;
  
  //Epot = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
  
#pragma omp parallel for reduction(+:Ekin) 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    //#pragma omp atomic
    Ekin += 0.5 * particle[atomi].Mass * (sqr(particle[atomi].vx) + sqr(particle[atomi].vy) + sqr(particle[atomi].vz));
  }
  
  fprintf(EFile,"%.15f %.16f %.16f %.16f\n",t*AKMA2SEC*1e15,Epot,Ekin,Epot+Ekin);
  fflush(EFile);
}

void Print_Traj(FILE *TrajFile,FILE *TrajFile_RI,FILE *TrajFile_Inter, int Natoms, double t, atom *particle,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,point L)
{
  int atomi;
  char atomtype[3];
  double Term;
  int Natoms_outer=0;
  
  for (atomi=0 ; atomi < Natoms ; atomi++)
    if(particle[atomi].layer==1)Natoms_outer++;
  
  if ((Xmin - Xmax) != 0.0 && (Ymin - Ymax) != 0.0  && (Zmin - Zmax) != 0.0){
    
    Xmin = Xmax = particle[0].r[0];
    Ymin = Ymax = particle[0].r[1];
    Zmin = Zmax = particle[0].r[2];
    
    //building periodic boundry box;

    for (atomi=0 ; atomi < Natoms ; atomi++){
      if (particle[atomi].r[0] < Xmin) Xmin = particle[atomi].r[0];
      if (particle[atomi].r[1] < Ymin) Ymin = particle[atomi].r[1];
      if (particle[atomi].r[2] < Zmin) Zmin = particle[atomi].r[2]; 
    }
    
    for (atomi=0 ; atomi < Natoms ; atomi++){
      if (particle[atomi].r[0] > Xmax) Xmax = particle[atomi].r[0];
      if (particle[atomi].r[1] > Ymax) Ymax = particle[atomi].r[1];
      if (particle[atomi].r[2] > Zmax) Zmax = particle[atomi].r[2];
      
      
    }
    
    Term = 0.5*(L.x - (Xmax - Xmin));
    Xmax += Term;
    Xmin -= Term;
    Term = 0.5*(L.y - (Ymax - Ymin));
    Ymax += Term;
    Ymin -= Term;
    Term = 0.5*(L.z - (Zmax - Zmin));
    Zmax += Term;
    Zmin -= Term;
    
    fprintf(TrajFile,"%i\n%.16f\n",Natoms + 8,t*AKMA2SEC*1e15);
    
    //#pragma omp parallel for private(atomtype)
    for(atomi=0 ; atomi < Natoms ; atomi++){
      if(particle[atomi].type == 1) strcpy(atomtype,"H");
      else if(particle[atomi].type == 5) strcpy(atomtype,"B");
      else if(particle[atomi].type == 7) strcpy(atomtype,"N");
      else if(particle[atomi].type == 0) strcpy(atomtype,"C");
      else{
	cerr<<"Wrong atom type (must be B, H, or N) in Print.c! Ending session.\n";
	exit(0);
      }
      
      fprintf(TrajFile,"%s %.16f %.16f %.16f\n",atomtype,particle[atomi].r[0],particle[atomi].r[1],particle[atomi].r[2]);
    }
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmin,Ymin,Zmin);
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmin,Ymin,Zmax);
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmin,Ymax,Zmin);
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmin,Ymax,Zmax);
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmax,Ymin,Zmin);
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmax,Ymin,Zmax);
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmax,Ymax,Zmin);
    fprintf(TrajFile, "F %.16f %.16f %.16f\n", Xmax,Ymax,Zmax);
  }
  else{
    fprintf(TrajFile,"%i\n%.16f\n",Natoms,t*AKMA2SEC*1e15);
    
    //#pragma omp parallel for firstprivate(Natoms) private(atomtype)
    for(atomi=0 ; atomi < Natoms ; atomi++){
      if(particle[atomi].type == 1) strcpy(atomtype,"H");
      else if(particle[atomi].type == 5) strcpy(atomtype,"B");
      else if(particle[atomi].type == 7) strcpy(atomtype,"N");
      else if(particle[atomi].type == 0) strcpy(atomtype,"C");
      else{
	cerr<<"Wrong atom type (must be B, H, or N) in Print.c! Ending session.\n";
	exit(0);
      }
      
      fprintf(TrajFile,"%s %.16f %.16f %.16f\n",atomtype,particle[atomi].r[0],particle[atomi].r[1],particle[atomi].r[2]);
    }
  }
  
  fprintf(TrajFile_RI,"%i\n%.16f\n",Natoms_outer,t*AKMA2SEC*1e15);
  fprintf(TrajFile_Inter,"%i\n%.16f\n",Natoms_outer,t*AKMA2SEC*1e15);
  
  //#pragma omp parallel for firstprivate(Natoms) private(atomtype)
  for(atomi=0 ; atomi < Natoms ; atomi++){
    if(particle[atomi].layer!=1)continue;
    
    if(particle[atomi].type == 1) strcpy(atomtype,"H");
    else if(particle[atomi].type == 5) strcpy(atomtype,"B");
    else if(particle[atomi].type == 7) strcpy(atomtype,"N");
    else if(particle[atomi].type == 0) strcpy(atomtype,"C");
    else{
      cerr<<"Wrong atom type (must be B, H, or N) in Print.c! Ending session.\n";
      exit(0);
    }
    
    fprintf(TrajFile_RI,"%s %.16f %.16f %.16f %.16f\n",atomtype,particle[atomi].r[0],particle[atomi].r[1],particle[atomi].r[2],particle[atomi].RI);
    fprintf(TrajFile_Inter,"%s %.16f %.16f %.16f %.16f\n",atomtype,particle[atomi].r[0],particle[atomi].r[1],particle[atomi].r[2],particle[atomi].Inter);
  }
  
  
  
  fflush(TrajFile);
  fflush(TrajFile_RI);
  fflush(TrajFile_Inter);
}

void Print_Traj_gv(int Natoms, atom *particle,point L)
{
  FILE *Traj_gv;
  FILE *Traj_jmol;
  Traj_gv = fopen("Optimized_Molecule_ReaxFF.xyz","w");
  Traj_jmol = fopen("Optimized_Molecule_jmol.xyz","w"); 
  int atomi;
  char atomtype[3];
  
  fprintf(Traj_jmol,"%i \n\n",Natoms);
  
  fprintf(Traj_gv,"%.16f %.16f %.16f\n",L.x,L.y,L.z);
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    
    if(particle[atomi].type == 1) strcpy(atomtype,"H");
    else if(particle[atomi].type == 5) strcpy(atomtype,"B");
    else if(particle[atomi].type == 7) strcpy(atomtype,"N");
    else if(particle[atomi].type == 0) strcpy(atomtype,"C");
    else{
      cerr<<"Wrong atom type (must be B, H, or N) in Print.c! Ending session.\n";
      exit(0);
    }
    
    fprintf(Traj_gv,"%i %i %.16f %.16f %.16f\n",particle[atomi].type,particle[atomi].layer,particle[atomi].r[0],particle[atomi].r[1],particle[atomi].r[2]);
    fprintf(Traj_jmol,"%s %.16f %.16f %.16f\n",atomtype,particle[atomi].r[0],particle[atomi].r[1],particle[atomi].r[2]);
  }
  
  fflush(Traj_gv);
  fclose(Traj_gv);
  fflush(Traj_jmol);
  fclose(Traj_jmol);
}
