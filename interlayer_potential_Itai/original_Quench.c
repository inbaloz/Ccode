#include "declarations.h"

/****************************************************************************/

void Velocity_Verlet_Quench(double dt, double dt2, point *accel, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point &L, point *Fterm , AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile, double *angle,Normal_struct *Normal,int *Normal_atom,int Interlayer, double *Fc,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,int *Neighb,int *NNeighb,int Periodic_flag,int Coulomb_flag,int *NList,int *List,point *R0,double alpha,double deltat,int &count,double &E_REBO,int *NBond_once,int *Bonded_once,int &counter)
{
  int atomi;
  //**
  int vx_unit,vz_unit,vy_unit;
  double P,v_norm,F_norm,Term1;
  //** 
  P      = 0.0;
  v_norm = 0.0;
  F_norm = 0.0;
  
#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    particle[atomi].r[0] += particle[atomi].vx*deltat + 0.5 * accel[atomi].x * sqr(deltat);
    particle[atomi].r[1] += particle[atomi].vy*deltat + 0.5 * accel[atomi].y * sqr(deltat);
    particle[atomi].r[2] += particle[atomi].vz*deltat + 0.5 * accel[atomi].z * sqr(deltat);
    
    particle[atomi].r[0] -= L.x * floor((particle[atomi].r[0] - Xmin) / (L.x + TINY));
    particle[atomi].r[1] -= L.y * floor((particle[atomi].r[1] - Ymin) / (L.y + TINY));
    particle[atomi].r[2] -= L.z * floor((particle[atomi].r[2] - Zmin) / (L.z + TINY));
    
  }
  
  if(!(counter%2000)){
    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once,Coulomb_flag);
    
    if (Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    if(Periodic_flag == 1)    
      relax_cell(Natoms,NBond,AtomParams,BondParams,Bonded,Bmat,qvec,Avec,rij,rijShiled_1,L,particle,Normal,Fc,logfile,angle,BondNum,Neighb,NNeighb,Coulomb_flag,Interlayer,Normal_atom,NList,List,R0,NBond_once,Bonded_once);
  }
  
  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once,Coulomb_flag);
  
  if (Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  Calc_Force(Fterm,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,E_REBO,NBond_once,Bonded_once);
 
#pragma omp parallel for reduction(+:P,v_norm,F_norm)
  for(atomi=0 ; atomi < Natoms ; atomi++){
    P     += particle[atomi].Fx*particle[atomi].vx+particle[atomi].Fy*particle[atomi].vy+particle[atomi].Fz*particle[atomi].vz;
    v_norm+= sqrt(sqr(particle[atomi].vx)+sqr(particle[atomi].vy)+sqr(particle[atomi].vz));
    F_norm+= sqrt(sqr(particle[atomi].Fx)+sqr(particle[atomi].Fy)+sqr(particle[atomi].Fz));
    accel[atomi].x = particle[atomi].Fx / 1.0;
    accel[atomi].y = particle[atomi].Fy / 1.0;
    accel[atomi].z = particle[atomi].Fz / 1.0;
    
  }
  
  if(F_norm==0.0)Term1=0.0;
  else Term1=alpha/(F_norm)*v_norm;
  
  if(P > 0){
#pragma omp parallel  for
    for(atomi=0 ; atomi < Natoms ; atomi++){
      particle[atomi].vx = (1-alpha)*particle[atomi].vx + particle[atomi].Fx*Term1;
      particle[atomi].vy = (1-alpha)*particle[atomi].vy + particle[atomi].Fy*Term1;
      particle[atomi].vz = (1-alpha)*particle[atomi].vz + particle[atomi].Fz*Term1;
    }
    
    if(count > 5){
      deltat=MIN(deltat*1.1,10*dt);
      alpha*=0.99;
    }
  }
  else{
#pragma omp parallel  for
    for(atomi=0 ; atomi < Natoms ; atomi++){
      particle[atomi].vx=particle[atomi].vy=particle[atomi].vz=0;
    }   
    alpha   = 0.1;
    deltat *= 0.5;
    count   = 0;
  }
#pragma omp parallel  for
  for(atomi=0 ; atomi < Natoms ; atomi++){
    particle[atomi].vx += accel[atomi].x*deltat;
    particle[atomi].vy += accel[atomi].y*deltat;
    particle[atomi].vz += accel[atomi].z*deltat;
  }
  
  count++;
}

/****************************************************************************/

void Quench(int Natoms, int *NBond, int *Bonded, atom *particle, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams,int *BondNum,double *rij, point &L, double dt, double tini, double tfin, double T, int NConserve, int NTraj, double *Bmat, double *qvec, double *dqvec, double *Avec, double *rijShiled_1, FILE* logfile, time_t &t0, clock_t &c0, double *angle, double OptTol,Normal_struct *Normal,int *Normal_atom,int Interlayer, double *Fc,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,int *Neighb,int *NNeighb,int Periodic_flag,int Coulomb_flag,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once)
{
  char FileName[50];
  int atomi, counter, file_counter, frame_counter, converged;
  double dt2, Eold, Enew;
  point *Fterm = new point[Natoms];
  point *accel = new point[Natoms];
  double alpha, deltat;
  int count;
  double E_REBO;
  alpha = 0.1;
  deltat= dt;
  count = 0;
  
  FILE *TrajFile = fopen("Traj-1.xyz","w");
  
  if(!TrajFile){
    cerr<<"cannot open trajectory file for writing\n";
    exit(0);
  }
  
  FILE *EFile    = fopen("E.dat","w");
  
  if(!EFile){
    cerr<<"cannot open energy file for writing\n";
    exit(0);
  }
  
  converged = 0;
  frame_counter = 0;
  file_counter = 1;
  
  //dt2 = sqr(dt)*10000000;

  dt2 =sqr(dt);
  
  if(Periodic_flag == 1)    
    relax_cell(Natoms,NBond,AtomParams,BondParams,Bonded,Bmat,qvec,Avec,rij,rijShiled_1,L,particle,Normal,Fc,logfile,angle,BondNum,Neighb,NNeighb,Coulomb_flag,Interlayer,Normal_atom,NList,List,R0,NBond_once,Bonded_once);
  
  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once,Coulomb_flag);
  
  if (Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  //Eold = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
  
  Calc_Force(Fterm,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,E_REBO,NBond_once,Bonded_once);

  Eold = E_REBO;
  
#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++){ // Calculate initial accelerations.
    accel[atomi].x = particle[atomi].Fx / 1.0;
    accel[atomi].y = particle[atomi].Fy / 1.0;
    accel[atomi].z = particle[atomi].Fz / 1.0;
  }
  
  fprintf(EFile,"%i %.16f\n",0,Eold);
  fflush(EFile);
  
  Print_Traj(TrajFile,Natoms,0,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
  frame_counter++;
  
  // Main propagation loop.
  
  for(counter=1 ; (!converged) ; counter++){
    
    Velocity_Verlet_Quench(dt,dt2,accel,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,Fterm,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Neighb,NNeighb,Periodic_flag, Coulomb_flag,NList,List,R0,alpha,deltat,count,E_REBO,NBond_once,Bonded_once,counter);
    
    Enew = E_REBO;
    
    //if((2.0 * fabs(Enew-Eold)) <= (OptTol * ( fabs(Eold) + fabs(Enew) + EPS ))) converged = 1;
    if(counter >= (OptTol-1.0)) converged = 1;
    else Eold = Enew;
    
    if(!(counter%NConserve)){
      fprintf(EFile,"%i %.16f\n",counter,Enew);
      fflush(EFile);
    }
    
    if(!(counter%NTraj)){
      Print_Traj(TrajFile,Natoms,counter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      Print_Traj_gv(Natoms,particle,L);
      
    }
  }
  
  fprintf(EFile,"%i %.16f\n",counter,Enew);
  fflush(EFile);

  Print_Traj(TrajFile,Natoms,counter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
  Print_Traj_gv(Natoms,particle,L);
  
  fclose(EFile);
  fclose(TrajFile);
  
  delete [] Fterm;
  delete [] accel;
}
