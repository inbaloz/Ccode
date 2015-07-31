#include "declarations.h"
void Randomize();
void Vel(atom *particle, double T, int Natoms, long *idum);

void Velocity_Verlet(double dt, double dt2, point *accel, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point &L, point *Fterm ,AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,double Xmin,double Ymin,double Zmin,Normal_struct *Normal ,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,point *R0,int *NList,int *List,double &E_REBO,int *NBond_once,int *Bonded_once)
{
  int atomi;
  //double axnew, aynew, aznew;
  
#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    particle[atomi].r[0] += particle[atomi].vx * dt + 0.5 * accel[atomi].x * dt2;
    particle[atomi].r[1] += particle[atomi].vy * dt + 0.5 * accel[atomi].y * dt2;
    particle[atomi].r[2] += particle[atomi].vz * dt + 0.5 * accel[atomi].z * dt2;
    
    particle[atomi].r[0] -= L.x * floor((particle[atomi].r[0] - Xmin) / (L.x + TINY));
    particle[atomi].r[1] -= L.y * floor((particle[atomi].r[1] - Ymin) / (L.y + TINY));
    particle[atomi].r[2] -= L.z * floor((particle[atomi].r[2] - Zmin) / (L.z + TINY));
  }
  
  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once,Coulomb_flag);
  
  if(Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  //cerr<<"hey00"<<endl;
  Calc_Force(Fterm,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,E_REBO,NBond_once,Bonded_once);
  //cerr<<"E_REBO_02"<<E_REBO<<endl;
#pragma omp parallel
  {
    double axnew, aynew, aznew;
#pragma omp for
    for(atomi=0 ; atomi < Natoms ; atomi++){
      axnew = particle[atomi].Fx / particle[atomi].Mass;
      aynew = particle[atomi].Fy / particle[atomi].Mass;
      aznew = particle[atomi].Fz / particle[atomi].Mass;
      
      particle[atomi].vx += 0.5 * (accel[atomi].x + axnew) * dt;
      particle[atomi].vy += 0.5 * (accel[atomi].y + aynew) * dt;
      particle[atomi].vz += 0.5 * (accel[atomi].z + aznew) * dt;
      
      accel[atomi].x = axnew;
      accel[atomi].y = aynew;
      accel[atomi].z = aznew;
    }
  }
}

void Propagate(int Natoms, int *NBond, int *Bonded, atom *particle, AtomParamsStruct *AtomParams, BondParamsStruct *BondParams, int *BondNum, double *rij, point &L, double dt, double tini, double tfin, double T, int NConserve, int NTraj, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile, time_t &t0, clock_t &c0,double *angle,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,Normal_struct *Normal ,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,point *R0,int *NList,int *List,int *NBond_once,int *Bonded_once)
{
  char FileName[50];
  int atomi, counter, file_counter, frame_counter;
  double t, dt2;
  double E_REBO;
  point *Fterm = new point[Natoms];
  point *accel = new point[Natoms];
  
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
  
  frame_counter = 1;
  file_counter = 1;
  
  dt2 = sqr(dt);
  
  Randomize();
  long idum = -random();
  
  vel(particle,T,Natoms,&idum); // Initialize velocities

  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once,Coulomb_flag);
  
  if(Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  Calc_Force(Fterm,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,E_REBO,NBond_once,Bonded_once);

  
  
#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++){ // Calculate initial accelerations.
    accel[atomi].x = particle[atomi].Fx / particle[atomi].Mass;
    accel[atomi].y = particle[atomi].Fy / particle[atomi].Mass;
    accel[atomi].z = particle[atomi].Fz / particle[atomi].Mass;
  }
  
  // Main propagation loop.

  for(counter = 0, t=tini ; t < tfin ; t += dt, counter++){
    
    Velocity_Verlet(dt,dt2,accel,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,Fterm,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Xmin,Ymin,Zmin,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,R0,NList,List,E_REBO,NBond_once,Bonded_once);
    
    //cerr<<"E_REBO_01"<<E_REBO<<endl;
    if(!(counter%NConserve))
      Print_Conserve(EFile,Natoms,NBond,Bonded,AtomParams,BondParams,particle,L,t,Bmat,qvec,Avec,rijShiled_1,logfile,t0,c0,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,Coulomb_flag,E_REBO,NBond_once,Bonded_once);
    //cerr<<"hey"<<endl;
    if(!(counter%NTraj)){
      if(!(frame_counter%Frame_limit)){
	fclose(TrajFile);
	file_counter++;
	sprintf(FileName,"Traj-%d.xyz",file_counter);
	TrajFile = fopen(FileName , "w");
	
	if(!TrajFile){
	  cerr<<"cannot open optimization file for writing\n";
	  exit(0);
	}
      }
      
      Print_Traj(TrajFile,Natoms,t,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      frame_counter++;
    }
  }
  
  fclose(EFile);
  fclose(TrajFile);
  
  delete [] Fterm;
  delete [] accel;
}
