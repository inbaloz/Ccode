#include "pair_airebo.h"
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
struct timeb t_initial,t_final,caliper[3];


int main()
{
  /********** Variable declarations **********/
  
  time_t  t0; /* time_t is defined on <time.h> and <sys/types.h> as long */
  clock_t c0; /* clock_t is defined on <time.h> and <sys/types.h> as int */
  int i,j,atomi,Interlayer,Periodic_flag, Coulomb_flag;
  int CalcType, Natoms, NConserve, NTraj,nthreads;
  int *NBond, *Bonded,*NBond_once, *Bonded_once, *BondNum,*List,*NList,Max_NBond,Max_NNeighb,Max_NList;
  int *Normal_atom, *NNeighb, *Neighb,p;
  double dt,tini, tend, T, OptTol;
  double *rij;
  double  *angle,*rijShiled_1,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Term;           //Coulomb term array.
  double *Bmat, *qvec, *dqvec, *Avec, *Fc; //Coulomb term arrays.
  char atomtype[3];
  point L;
  point *R0;
  point *Fbond;
  FILE *logfile;
  long idum;
  Randomize();
  idum = -random();
  
  AtomParamsStruct AtomParams[NAtomParams];
  BondParamsStruct BondParams[NAtomParams*NAtomParams];
    
  /********** Initialization **********/
  
  logfile = fopen("log.dat","w");
  
  Init_Atom_Params(AtomParams);
  Init_Bond_Params(BondParams);
    
  Natoms              = GetNatoms("Coords.xyz");
  atom *particle      = new atom[Natoms];
  atom *Temp_atom     = new atom[Natoms];
  NBond               = new int[Natoms];
  NBond_once               = new int[Natoms];
  BondNum             = new int[1];
  NNeighb             = new int[Natoms];
  
  Neighb              = new int[Natoms*MaxNeighb];
  
  List                = new int[Natoms*MaxList];
  NList               = new int[Natoms];
  //r_ij_Neighb         = new int[Natoms*MaxNeighb];
  //**
  Fc                  = new double[Natoms*MaxNeighb];
  //**
  nC             = new double[Natoms];
  nH             = new double[Natoms];
  Bonded              = new int[Natoms*MaxNBond];
  Bonded_once              = new int[Natoms*MaxNBond];
  rij                 = new double[Natoms*MaxNBond];
  rijShiled_1         = new double[Natoms*MaxNeighb];
  //rijShiled_1         = new double[1];
  Fbond               = new point[Natoms];
  R0                  = new point[Natoms];
  
  // Arrays needed for the Normal of the atoms
  Normal_struct *Normal = new Normal_struct[Natoms];
  // Normal_struct *Normal = new point[Natoms];
  Normal_atom = new int[Natoms*3];
   
  // Arrays required by the Coulomb term calculations.
  //Bmat             = new double[(Natoms+1)*(Natoms+1)];
  Bmat             = new double[1];
  qvec             = new double[Natoms+1];
  dqvec            = new double[Natoms+1];
  Avec             = new double[Natoms+1];
  //cerr<<"yo"<<endl;
  // Arrays required by the valence angle term calculations.
  //angle            = new double[Natoms*MaxNBond*MaxNBond];
  angle            = new double[1];

 
  
  read_RUN_parameters(tini,tend,dt,L,T,NConserve,NTraj,OptTol,CalcType,Interlayer,Periodic_flag,nthreads,Coulomb_flag);
 
  spline_init();
  

  
  
  ReadCoor("Coords.xyz",Natoms,particle,Interlayer,L);
  fprintf(logfile,"yo yo 1");  
  read_file("CH.airebo");
  fprintf(logfile,"yo yo 2");  
    
  Xmin = Xmax = particle[0].r[0];
  Ymin = Ymax = particle[0].r[1];
  Zmin = Zmax = particle[0].r[2];
 
  //building periodic boundry box; 
  for (atomi=0 ; atomi < Natoms ; atomi++){
    R0[atomi].x=particle[atomi].r[0];
    R0[atomi].y=particle[atomi].r[1];
    R0[atomi].z=particle[atomi].r[2];
    
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
  
  /////start timeing****
  
  if (nthreads != 0)
    omp_set_num_threads(nthreads);
  
  int milestone(0), rc;
  rc = ftime(&caliper[milestone]);
  
  for(atomi=0 ; atomi < Natoms ; atomi++){
    NBond[atomi]   = 0;
    NNeighb[atomi] = 0;
    NList[atomi]   = 0;
    particle[atomi].RI=0.0;
    particle[atomi].Inter=0.0;
    dqvec[atomi]         =0.0;
  }
  //**check number of neihbors***
  /*
    double Tersoff_R,Tersoff_D,r_ij;
    int atomj,typej,typei;
    
    
    //#pragma omp parrallel for private(typei,typej,Tersoff_R,Tersoff_D,atomaj,r_ij)
    for(atomi=0 ; atomi < Natoms ; atomi++){
    typei = particle[atomi].type;
    for(atomj=0 ; atomj < atomi  ; atomj++){
    typej=particle[atomj].type;
    
    r_ij=R_ij(particle,atomi,atomj,L);
    
    //if(r_ij < Tersoff_R+Tersoff_D){
    if(r_ij < sqrt(rcmaxsq[typei][typej])){
    NBond[atomi]++;
    NBond[atomj]++;
    if(NBond[atomi] > Max_NBond)Max_NBond=NBond[atomi];
    if(NBond[atomj] > Max_NBond)Max_NBond=NBond[atomj];
    }
    
    if((r_ij < non_bond_cut) && (particle[atomi].layer != particle[atomj].layer)){
    NNeighb[atomi]++;
    if(NNeighb[atomi] > Max_NNeighb) Max_NNeighb=NNeighb[atomi];
    }
    
    if(r_ij < Rcut_List){
    NList[atomi]++;
    if(NList[atomi] > Max_NList)Max_NList=NList[atomi];
    }
    
    }
    }      
    
    
    fprintf(logfile,"number of Neigbors at initial config\n Max_NBond=%i Max_NNeighb=%i Max_NList=%i\n",Max_NBond,Max_NNeighb,Max_NList);
    fflush(logfile);
  */
  //*****************************
  
  UpdateNeighbList(Natoms,particle,List,NList,R0,L);
  calc_rij(NBond,Bonded,BondNum,Natoms,particle,L,rij,Fc,BondParams,Neighb,NNeighb,rijShiled_1,NList,List,NBond_once,Bonded_once);
  
  
  if (Interlayer == 1)
    set_Normal_atom(rij,Natoms,Normal,particle,Normal_atom,L);//The normal is calculating according to an atoms three nearest neighbors, Sp2 bonding is required appart from hydrogen atoms.
  else if (Interlayer == 0){
#pragma omp parallel for
    for(atomi=0; atomi < Natoms; atomi++){
      particle[atomi].layer = 0;
      Normal[atomi].x = 0.0;
      Normal[atomi].y = 0.0;
      Normal[atomi].z = 0.0;
    }
  }
  else if(Interlayer == 2){
#pragma omp parallel for  
    for(atomi=0; atomi < Natoms; atomi++){
      Normal[atomi].x = 0.0;
      Normal[atomi].y = 0.0;
      Normal[atomi].z = 1.0;
    }
  }
 
  //********
  //Calc_Charge(particle,Natoms,Bmat,qvec,Avec,AtomParams,BondParams,rij,rijShiled_1,logfile,L);
  /********** Calculate **********/

  //for (atomi=0; atomi < Natoms; atomi++){
  //if(particle[atomi].layer==1)
  // particle[atomi].r[0]-=0.715;
  //}
  //**
  //ReadCoor("Temp.xyz",Natoms,particle,Interlayer);
  //**
  cerr<<"hey 00"<<endl;
  if(CalcType == 0){
    vel(particle,T,Natoms,&idum); 
    Propagate(Natoms,NBond,Bonded,particle,AtomParams,BondParams,BondNum,rij,L,dt,tini,tend,T,NConserve,NTraj,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,t0,c0,angle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,R0,NList,List,NBond_once,Bonded_once);
  }
  else if(CalcType == 1){
    double Energy,Energy_;
    FILE *Single_point;
    FILE *TrajFile = fopen("Single_point_traj.xyz","w");
    FILE *TrajFile_RI = fopen("Traj-RI.xyz","w");
    FILE *TrajFile_Inter = fopen("Traj-Inter.xyz","w");
    
    point *Fterm = new point[Natoms];
    Single_point = fopen("Single_point_Energy.dat","w");
    
    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once);
    
    if(Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    Energy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Energy_,NBond_once,Bonded_once,0,0,Single_point);
    //Calc_Force(Fterm,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Energy,NBond_once,Bonded_once);
    
    fprintf(Single_point,"ToT_Energy= %.16f\nEnergy_per_atom= %.16f\n",Energy,Energy/Natoms);
    
    Print_Traj(TrajFile,TrajFile_RI,TrajFile_Inter,Natoms,0,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
    
    fclose(TrajFile_RI);
    fclose(TrajFile_Inter);
    fclose(TrajFile);
    fclose(Single_point);
    delete [] Fterm;
  }
  else if(CalcType == 2){
    //steepest_descent(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,NTraj,NConserve,OptTol,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Normal,Normal_atom,Interlayer,Fc,Periodic_flag,Neighb,NNeighb,Coulomb_flag,NList,List,R0); 
  }
  else if(CalcType == 3){
    //conj_grad(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,OptTol,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,NConserve, NTraj,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Normal,Normal_atom,Interlayer,Fc,Periodic_flag,Neighb,NNeighb,Coulomb_flag,NList,List,R0);
  }
  else if(CalcType ==4){
    vel(particle,T,Natoms,&idum); 
    Quench( Natoms,NBond,Bonded,particle,AtomParams,BondParams,BondNum,rij,L,dt,tini,tend,T,NConserve,NTraj,Bmat,qvec,dqvec,Avec,rijShiled_1,logfile,t0,c0,angle,OptTol,Normal,Normal_atom,Interlayer,Fc,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,Neighb,NNeighb,Periodic_flag,Coulomb_flag,NList,List,R0,NBond_once,Bonded_once);
  }
  else{
    cerr<<"Unrecognized requested calculation type! Ending session.\n";
    exit(0);
  }
  
  
    //set_Normal_atom(rij,Natoms,Normal,particle,Normal_atom,L);
    //calc_rij(NBond,Bonded,BondNum,Natoms,particle,L,rij,Fc,BondParams,Neighb,NNeighb,rijShiled_1,NList,List);  
    //calc_Normal(Natoms,Normal_atom,particle,Normal,L);  
    //Calc_Force(Fterm,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Min_Epot);
  
  milestone++;
  rc = ftime(&caliper[milestone]);
  
  double delta_t(0);
      int hours,minutes,seconds,time;
      for (int i=0;i<milestone;++i)
	{
	  delta_t = (caliper[i+1].time-caliper[i].time);
	  delta_t += (caliper[i+1].millitm-caliper[i].millitm)/1000.0;
	  //printf("time=%-.3f seconds using %i threads\n",delta_t,omp_get_max_threads());
	  time    = int(delta_t);
	  hours   = (time/3600);
	  time    = time%3600;
	  minutes = time/60;
	  seconds  = time%60;
	  //seconds = ((time%3600)/60)%60;
	  //printf(logfile,"time=%i:%i:%i h:m:s using %i threads\n",hours,minutes,seconds,omp_get_max_threads());
	  fprintf(logfile,"time=%i:%i:%i h:m:s using %i threads\n",hours,minutes,seconds,omp_get_max_threads());
	}      
      
      //**********************************************

      fclose(logfile);
      delete [] particle;
      delete [] Temp_atom;
      delete [] NBond;
      delete [] Bonded;
      delete [] BondNum;
      delete [] rij;
      delete [] Fbond;
      delete [] angle;
      delete [] rijShiled_1;
      delete [] NBond_once;
      delete [] Bonded_once;
      // Deallocate Coulomb term related arrays
      //delete [] Fterm;
      delete [] Bmat;
      delete [] qvec;
      delete [] dqvec;
      delete [] Avec;
}


/**********************************************************************/
/* Calculate dr_{ij}_dr_k                                             */
void calc_rij(int *NBond, int *Bonded, int *BondNum, int Natoms, atom *particle, point L, double *rij,double *Fc,BondParamsStruct *BondParams,int *Neighb,int *NNeighb,double *rijShiled_1,int *NList,int *List,int *NBond_once,int *Bonded_once)
{
  //#pragma omp parallel firstprivate(Natoms) shared(Natoms,NBond,Bonded,particle,L,rij,Fc,BondNum)
  //#pragma omp parallel
  //{
  //{
  int atomi, atomj,typei,typej,j;
  double rx, ry, rz,r_ij;
  double Tersoff_R,Tersoff_D;
  int NBondj, NBondi,NBondi_once;

  double dS,rsq;;
  
#pragma omp parallel for private(atomj)
  for(atomi=0 ; atomi < Natoms ; atomi++){
    NBond[atomi]   = 0;
    NBond_once[atomi]   = 0;
    NNeighb[atomi] = 0;
    nC[atomi]=0;
    nH[atomi]=0;
    for(atomj=0 ; atomj < MaxNBond ; atomj++){
      Bonded[atomi*MaxNBond + atomj] = -1;
      Bonded_once[atomi*MaxNBond + atomj] = -1;
      //rij[atomi*MaxNBond + atomj] = 0.0;
    }
    for(atomj=0 ; atomj < MaxNeighb ; atomj++){
      Neighb[atomi*MaxNeighb + atomj] = -1;
      //rijShiled_1[atomi*MaxNeighb + atomj] = 0.0;
    }
  }

  
#pragma omp parallel for private(rx,ry,rz,typei,typej,Tersoff_R,Tersoff_D,atomj,r_ij,j,atomi,NBondj, NBondi, NBondi_once,dS) 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    
    //rij[atomi*Natoms + atomi] = 0.0;
    typei = particle[atomi].type;
      for(j=0 ; j < NList[atomi] ; j++){
	atomj=List[atomi*MaxList + j];
	typej=particle[atomj].type;
	
	r_ij=R_ij(particle,atomi,atomj,L);
	//rsq=sqr(r_ij);
	//rij[atomi*Natoms + atomj] = sqrt(sqr(particle[atomi].r[0] - particle[atomj].r[0] - L.x * rint(rx / (L.x + TINY))) + sqr(particle[atomi].r[1] - particle[atomj].r[1] - L.y * rint(ry / (L.y + TINY))) + sqr(particle[atomi].r[2] - particle[atomj].r[2] - L.z * rint(rz / (L.z + TINY))));
      //rij[atomj*Natoms + atomi] = rij[atomi*Natoms + atomj];
      
      if(r_ij < sqrt(rcmaxsq[typei][typej])){

		
#pragma omp critical
	{
	  NBond[atomi]++;
	  NBondi = NBond[atomi];
	  NBond_once[atomi]++;
	  NBondi_once = NBond_once[atomi];
	  //#pragma omp atomic
	  NBond[atomj]++;
	  NBondj = NBond[atomj];
	}

	Bonded_once[atomi*MaxNBond + NBondi_once-1] = atomj;
	Bonded[atomi*MaxNBond + NBondi-1] = atomj;
	Bonded[atomj*MaxNBond + NBondj-1] = atomi;
	
	
		
	//#pragma omp critical
	//rij[atomi*MaxNBond + NBondi-1] = rij[atomj*MaxNBond + NBondj-1] = r_ij;
	rij[atomi*MaxNBond + NBondi_once-1] = r_ij;
	//#pragma omp critical
	//Fc[atomi*MaxNBond + NBondi-1]  = Fc[atomj*MaxNBond + NBondj-1]  = Fc_(r_ij,Tersoff_R,Tersoff_D);
#pragma omp critical
	{
	  if (typej == 0){
	    nC[atomi] += Sp(r_ij,rcmin[typei][typej],rcmax[typei][typej],dS);
	    nC[atomj] += Sp(r_ij,rcmin[typej][typei],rcmax[typej][typei],dS);
	  }
	  else{
	    nH[atomi] += Sp(r_ij,rcmin[typei][typej],rcmax[typei][typej],dS);
	    nH[atomj] += Sp(r_ij,rcmin[typej][typei],rcmax[typej][typei],dS);
	  }
	}
	
	if(NBond[atomi] == MaxNBond ){
	  cerr<<"Exceeded maximum number of bonded atoms for atom "<<atomi<<"! Ending session.\n";
	  exit(0);
	}
	
	/*
	if(NBond[atomj] == MaxNBond ){
	  cerr<<"Exceeded maximum number of bonded atoms for atom "<<atomj<<"! Ending session.\n";
	  exit(0);
	}
	*/
      }
      
      if((r_ij < non_bond_cut) && (particle[atomi].layer != particle[atomj].layer)){
      //if((r_ij < non_bond_cut)){
	Neighb[atomi*MaxNeighb + NNeighb[atomi]]      = atomj;
	rijShiled_1[atomi*MaxNeighb + NNeighb[atomi]] = r_ij;
	NNeighb[atomi]++;
	if(NNeighb[atomi] == MaxNeighb){
	  cerr<<"Exceeded maximum number of neighbors for atom "<<atomi<<"! Ending session.\n";
	  exit(0);
	}
      }
      
      }      
  }
  
}
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

/****************************************************************/
/*  dThetaijk_drn                                               */

// Calculating the derivative with respect to movments of atom-n of
// the valence angle {ijk}.

double calc_dthetaijk_drn(int atomi, int atomj, int atomk, int atomn, int dir, atom *particle, double *angle, int *BondNum, int Natoms, double *rij, point L,double Thetaijk,double r_jk, double r_ij,dr_ij &rji,dr_ij &rjk,double cosijk)
{
  int indexji, indexjk;
  double ri, rj, rk, rjk_1, rji_1;
  double Length;
  
  indexji = atomj*Natoms + atomi;
  indexjk = atomj*Natoms + atomk;
  
  if((fabs(Thetaijk) > EPS) && (fabs(Thetaijk - PIE) > EPS)){
    if(atomn == atomi){
      
      rjk_1 = 1.0 /r_jk;
      rji_1 = 1.0 /r_ij;
      
      return( (1.0 / sqrt(1.0 - sqr(cosijk))) * rji_1 * (rji.r[dir] * rji_1 * cosijk - rjk.r[dir] * rjk_1) );
    }
    else if(atomn == atomk){
      
      rjk_1 = 1.0 / r_jk;
      rji_1 = 1.0 / r_ij;
      
      return( (1.0 / sqrt(1.0 - sqr(cosijk))) * rjk_1 * (rjk.r[dir] * rjk_1 * cosijk -  rji.r[dir] * rji_1) );
    }
    else if(atomn == atomj){
      
      rjk_1 = 1.0 / r_jk;
      rji_1 = 1.0 / r_ij;
      
      return( - (1.0 / sqrt(1.0 - sqr(cosijk))) * (rji_1 * (rji.r[dir] * rji_1 * cosijk - rjk.r[dir] * rjk_1) + rjk_1 * (rjk.r[dir] * rjk_1 * cosijk - rji.r[dir] * rji_1)) );
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

void calc_dE_REBO(point *F_REBO,double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_REBO_,int *NBond_once,int *Bonded_once)
{
  double E_REBO=0.0;
#pragma omp parallel
  {
    int index,atomk,typek,typej,typei,index2;
    double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
    double Xsi,ETersoff,dTheta_ijk_drn,dg_drn,dXsi_drn,Fc_ij,Fc_jk,Tersoff_R,Tersoff_D;
    double dr_ij_drn,dr_jk_drn,dFc_ij_drn,dFc_jk_drn,Fa,Fr,dFa_drn,dFr_drn,db_ij_drn,n,dE_Tersoff_ij_drn;
    double Tersoff_sqr_c,Tersoff_sqr_d,sqr_Tersoff_h_Cos_Theta;
    //i=BondNum[atomj*Natoms + atomi];
    double inside_acos,cos_ijk;
    //cerr<<"min[0][0]2 ="<<rcmin[0][0]<<endl;
    //exit(0);
    
    int i,j,k,m,ii,itype,jtype,atomi,atomj;
    // tagint itag,jtag;
    int itag ,jtag;
    double delx,dely,delz,evdwl,fpair,xtmp,ytmp,ztmp;
    double rsq,rij,wij;
    double Qij,Aij,alphaij,VR,pre,dVRdi,VA,term,bij,dVAdi,dVA;
    double dwij,del[3];
    //double E_REBO;
    
    int vflag_atom=0;
    evdwl = 0.0;
    
    //double **x = atom->x;
    //double **f = atom->f;
    //int *type = atom->type;
    //tagint *tag = atom->tag;
    //int nlocal = atom->nlocal;
    //int newton_pair = force->newton_pair;
    
    //Natoms = list->Natoms;
    //ilist = list->ilist;
  
    // two-body interactions from REBO neighbor list, skip half of them
#pragma omp for //reduction(+:E_REBO)
    for (atomi = 0; atomi < Natoms; atomi++) {
      itag=atomi;
      itype = particle[atomi].type;
      xtmp = particle[atomi].r[0];
      ytmp = particle[atomi].r[1];
      ztmp = particle[atomi].r[2];
      //REBO_neighs = REBO_firstneigh[i];
      //cerr<<"atomi= "<<atomi<<endl;
      for (j = 0; j < NBond_once[atomi]; j++) {
	atomj = Bonded_once[atomi*MaxNBond + j];
	//cerr<<"atomj= "<<atomj<<endl;
	//j = REBO_neighs[k];
	//jtag = tag[j];


	/*
	if (itag > jtag) {
	  if ((itag+jtag) % 2 == 0) continue;
	} else if (itag < jtag) {
	  if ((itag+jtag) % 2 == 1) continue;
	} else {
	  if (particle[atomj].r[2] < ztmp) continue;
	  if (particle[atomj].r[2] == ztmp && particle[atomj].r[1] < ytmp) continue;
	  if (particle[atomj].r[2] == ztmp && particle[atomj].r[1] == ytmp && particle[atomj].r[0] < xtmp) continue;
	}
	*/
	//jtype = map[type[j]];
	jtype = particle[atomj].type;
	
	//delx = particle[atomi].r[0] - particle[atomj].r[0];
	//dely = particle[atomi].r[1] - particle[atomj].r[1];
	//delz = particle[atomi].r[2] - particle[atomj].r[2];
	delx = R_PBC(particle,atomi,atomj,0,L.x);
	dely = R_PBC(particle,atomi,atomj,1,L.y);	
	delz = R_PBC(particle,atomi,atomj,2,L.z);

	
	rsq = delx*delx + dely*dely + delz*delz;
	rij = sqrt(rsq);
	
	wij = Sp(rij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
	
	if (wij <= TOL) continue;

	
	//Qij= BondParams[index].REBO_Q;
	Qij = Q[itype][jtype];
	//Aij = A[itype][jtype];
	Aij = A_ij[itype][jtype];
	//Aij= BondParams[index].REBO_A;
	//Aij = A[itype][jtype];
	//alphaij= BondParams[index].REBO_alpha;
	alphaij = alpha[itype][jtype];
	
	VR = wij*(1.0+(Qij/rij)) * Aij*exp(-alphaij*rij);
	//printf("VR=%.16f wij=%.16f Qij=%.16f Aij=%.16f alphaij=%.16f rij=%.16f\n",VR,wij,Qij,Aij,alphaij,rij);
	pre = wij*Aij * exp(-alphaij*rij);
	dVRdi = pre * ((-alphaij)-(Qij/rsq)-(Qij*alphaij/rij));
	dVRdi += VR/wij * dwij;
	
	VA = dVA = 0.0;
	
	
	//for (m = 0; m < 3; m++) {
	//term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rij);
	//VA += term;
	//dVA += -Beta[itype][jtype][m] * term;
	//}
	
	
	//term = -wij * BIJc1[itype][jtype][0] * exp(-Beta1[itype][jtype][0]*rij) -wij * BIJc2[index] * exp(-Beta2[index]*rij) - wij * BIJc3[index] * exp(-Beta3[index]*rij);
	//VA   = term;
	//dVA  = - Beta1[index] * term - Beta2[index] * term - Beta3[index] * term;
	
	for (m = 0; m < 3; m++) {
	  term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rij);
	  VA += term;
	  //printf("VA=%.16f wij=%.16f  BIJc=%.16f Beta=%.16f rij=%.16f term=%.16f\n",VA,wij,BIJc[itype][jtype][m] ,Beta[itype][jtype][m],rij,term);
	  dVA += -Beta[itype][jtype][m] * term;
	}
	
	dVA += VA/wij * dwij;
	
	del[0] = delx;
	del[1] = dely;
	del[2] = delz;
	
	bij = bondorder(atomi,atomj,del,rij,VA,F_REBO,vflag_atom,particle,Bonded,BondNum,NBond,L);
	//bij=1;
	dVAdi = bij*dVA;
	
	fpair = -(dVRdi+dVAdi) / rij;
	//fpair = -(dVAdi) / rij;
	//fpair = -(dVRdi) / rij;
#pragma omp critical
	{
	F_REBO[atomi].x += delx*fpair;
	F_REBO[atomi].y += dely*fpair;
	F_REBO[atomi].z += delz*fpair;
	F_REBO[atomj].x -= delx*fpair;
	F_REBO[atomj].y -= dely*fpair;
	F_REBO[atomj].z -= delz*fpair;
	E_REBO+=VR + bij*VA;
	}

	
	
      
	//printf("E_REBO=%.16f VR=%.16f VA=%.16f bij=%.16f\n",0.5*E_REBO,VR,VA,bij);
	//E_REBO += VR;
	//E_REBO+=bij*VA;
	//cerr<<E_REBO<<endl;
	//if (eflag) evdwl = VR + bij*VA;
	//if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
  //cerr<<"E_REBO_03"<<E_REBO<<endl;
  E_REBO_=E_REBO;
}


double bondorder(int i, int j, double rij[3], double rijmag, double VA, point *f, int vflag_atom,atom *particle,int *Bonded,int *BondNum,int *NBond,point L){

  int atomi,atomj,k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int itype,jtype,ktype,ltype,ntype;
  double rik[3],rjl[3],rkn[3],rji[3],rki[3],rlj[3],rknmag,dNki,dwjl,bij;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ,Nki,Nlj,dS;
  double lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosjikdrj[3],dcosijldrj[3],dcosijldrl[3];
  double Tij;
  double r32[3],r32mag,cos321,r43[3],r13[3];
  double dNlj;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double dcut321,PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;

  //double **x = atom->x;
  //int *type = atom->type;

  atomi = i;
  atomj = j;
  itype = particle[atomi].type;
  jtype = particle[atomj].type;

  //index      = itype*NAtomParams + jtype;
  wij = Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  //printf("jtype=%i itype=%i wij=%.16f\n",jtype,itype,wij);  
  //  wij = Sp(rijmag,rcmin[index],rcmax[index],dwij);
  NijC = nC[i]-(wij*kronecker(jtype,0));
  NijH = nH[i]-(wij*kronecker(jtype,1));
  NjiC = nC[j]-(wij*kronecker(itype,0));
  NjiH = nH[j]-(wij*kronecker(itype,1));

  //printf("jtype=%i itype=%i wij=%.16f nC[i]=%.16f nH[i]=%.16f nC[j]=%.16f nH[j]=%.16f \n",jtype,itype,wij,nC[i],nH[i],nC[j],nH[j]);
  
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  
  //REBO_neighs = REBO_firstneigh[i];
  
    for(k=0 ; k < NBond[atomi] ; k++){ // Go over all j neighbors k != i.
      //atomk = REBO_neighs[k];
      atomk = Bonded[atomi*MaxNBond + k];
    if (atomk != atomj) {
      ktype = particle[atomk].type;
      //ktype = particle[atomk].type;
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) *
        ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
        (wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp = Etmp+(wik*g*exp(lamdajik));
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS)); //parameter Nmin Nmax
    }
  }

  PijS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2);
  pij = pow(1.0+Etmp+PijS,-0.5);
  tmp = -0.5*cube(pij);
  //printf("pij=%.16f Etmp=%.16f PijS=%.16f \n",pij,Etmp,PijS);
  // pij forces
  
  //REBO_neighs = REBO_firstneigh[i];
  for (k = 0; k < NBond[atomi]; k++) {
    //atomk = REBO_neighs[k];
    atomk = Bonded[atomi*MaxNBond + k];
    if (atomk != atomj) {
      ktype = particle[atomk].type;
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) *
        ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      cosjik = (rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2]) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      dcosjikdri[0] = ((rij[0]+rik[0])/(rijmag*rikmag)) -
        (cosjik*((rij[0]/(rijmag*rijmag))+(rik[0]/(rikmag*rikmag))));
      dcosjikdri[1] = ((rij[1]+rik[1])/(rijmag*rikmag)) -
        (cosjik*((rij[1]/(rijmag*rijmag))+(rik[1]/(rikmag*rikmag))));
      dcosjikdri[2] = ((rij[2]+rik[2])/(rijmag*rikmag)) -
        (cosjik*((rij[2]/(rijmag*rijmag))+(rik[2]/(rikmag*rikmag))));
      dcosjikdrk[0] = (-rij[0]/(rijmag*rikmag)) +
        (cosjik*(rik[0]/(rikmag*rikmag)));
      dcosjikdrk[1] = (-rij[1]/(rijmag*rikmag)) +
        (cosjik*(rik[1]/(rikmag*rikmag)));
      dcosjikdrk[2] = (-rij[2]/(rijmag*rikmag)) +
        (cosjik*(rik[2]/(rikmag*rikmag)));
      dcosjikdrj[0] = (-rik[0]/(rijmag*rikmag)) +
        (cosjik*(rij[0]/(rijmag*rijmag)));
      dcosjikdrj[1] = (-rik[1]/(rijmag*rikmag)) +
        (cosjik*(rij[1]/(rijmag*rijmag)));
      dcosjikdrj[2] = (-rik[2]/(rijmag*rikmag)) +
        (cosjik*(rij[2]/(rijmag*rijmag)));

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wik*dgdc*exp(lamdajik));
      fj[0] = -tmp2*dcosjikdrj[0];
      fj[1] = -tmp2*dcosjikdrj[1];
      fj[2] = -tmp2*dcosjikdrj[2];
      fi[0] = -tmp2*dcosjikdri[0];
      fi[1] = -tmp2*dcosjikdri[1];
      fi[2] = -tmp2*dcosjikdri[2];
      fk[0] = -tmp2*dcosjikdrk[0];
      fk[1] = -tmp2*dcosjikdrk[1];
      fk[2] = -tmp2*dcosjikdrk[2];

      tmp2 = VA*.5*(tmp*wik*g*exp(lamdajik)*4.0*kronecker(itype,1));
      fj[0] -= tmp2*(-rij[0]/rijmag);
      fj[1] -= tmp2*(-rij[1]/rijmag);
      fj[2] -= tmp2*(-rij[2]/rijmag);
      fi[0] -= tmp2*((-rik[0]/rikmag)+(rij[0]/rijmag));
      fi[1] -= tmp2*((-rik[1]/rikmag)+(rij[1]/rijmag));
      fi[2] -= tmp2*((-rik[2]/rikmag)+(rij[2]/rijmag));
      fk[0] -= tmp2*(rik[0]/rikmag);
      fk[1] -= tmp2*(rik[1]/rikmag);
      fk[2] -= tmp2*(rik[2]/rikmag);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwik*g*exp(lamdajik))/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ktype]*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwik)/rikmag;
      fi[0] -= tmp2*rik[0];
      fi[1] -= tmp2*rik[1];
      fi[2] -= tmp2*rik[2];
      fk[0] += tmp2*rik[0];
      fk[1] += tmp2*rik[1];
      fk[2] += tmp2*rik[2];
#pragma omp critical
	{
      f[atomi].x += fi[0]; f[atomi].y += fi[1]; f[atomi].z += fi[2];
      f[atomj].x += fj[0]; f[atomj].y += fj[1]; f[atomj].z += fj[2];
      f[atomk].x += fk[0]; f[atomk].y += fk[1]; f[atomk].z += fk[2];
	}
      /*
      if (vflag_atom) {
        rji[0] = -rij[0]; rji[1] = -rij[1]; rji[2] = -rij[2];
        rki[0] = -rik[0]; rki[1] = -rik[1]; rki[2] = -rik[2];
        v_tally3(atomi,atomj,atomk,fj,fk,rji,rki);
	}
      */
    }
  }

  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  //REBO_neighs = REBO_firstneigh[j];
  
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml != atomi) {
      ltype = particle[atoml].type;
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) *
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0)) +
        nH[atoml]-(wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp = Etmp+(wjl*g*exp(lamdaijl));
      //printf("Etmp=%.16f wjl=%.16f g=%.16f lamdaijl=%.16f\n",Etmp,wjl,g,lamdaijl);
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl));
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  //NjiC =1.0;
  //NjiH =1.0;
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2);
  //printf("PjiS=%.16f NjiC=%.16f NjiH=%.16f dN2[0]=%.16f dN2[1]=%.16f   jtype=%i itype=%i\n",PjiS,NjiC,NjiH,dN2[0],dN2[1],jtype,itype);
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmp = -0.5*cube(pji);
  //printf("pji=%.16f Etmp=%.16f PjiS=%.16f \n",pji,Etmp,PjiS);
  
  //REBO_neighs = REBO_firstneigh[j];
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml != atomi) {
      ltype = particle[atoml].type;
      
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) *
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      cosijl = (-1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2]))) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      dcosijldri[0] = (-rjl[0]/(rijmag*rjlmag)) -
        (cosijl*rij[0]/(rijmag*rijmag));
      dcosijldri[1] = (-rjl[1]/(rijmag*rjlmag)) -
        (cosijl*rij[1]/(rijmag*rijmag));
      dcosijldri[2] = (-rjl[2]/(rijmag*rjlmag)) -
        (cosijl*rij[2]/(rijmag*rijmag));
      dcosijldrj[0] = ((-rij[0]+rjl[0])/(rijmag*rjlmag)) +
        (cosijl*((rij[0]/square(rijmag))-(rjl[0]/(rjlmag*rjlmag))));
      dcosijldrj[1] = ((-rij[1]+rjl[1])/(rijmag*rjlmag)) +
        (cosijl*((rij[1]/square(rijmag))-(rjl[1]/(rjlmag*rjlmag))));
      dcosijldrj[2] = ((-rij[2]+rjl[2])/(rijmag*rjlmag)) +
        (cosijl*((rij[2]/square(rijmag))-(rjl[2]/(rjlmag*rjlmag))));
      dcosijldrl[0] = (rij[0]/(rijmag*rjlmag))+(cosijl*rjl[0]/(rjlmag*rjlmag));
      dcosijldrl[1] = (rij[1]/(rijmag*rjlmag))+(cosijl*rjl[1]/(rjlmag*rjlmag));
      dcosijldrl[2] = (rij[2]/(rijmag*rjlmag))+(cosijl*rjl[2]/(rjlmag*rjlmag));

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      tmp2 = VA*.5*(tmp*wjl*dgdc*exp(lamdaijl));
      fi[0] = -tmp2*dcosijldri[0];
      fi[1] = -tmp2*dcosijldri[1];
      fi[2] = -tmp2*dcosijldri[2];
      fj[0] = -tmp2*dcosijldrj[0];
      fj[1] = -tmp2*dcosijldrj[1];
      fj[2] = -tmp2*dcosijldrj[2];
      fl[0] = -tmp2*dcosijldrl[0];
      fl[1] = -tmp2*dcosijldrl[1];
      fl[2] = -tmp2*dcosijldrl[2];

      tmp2 = VA*.5*(tmp*wjl*g*exp(lamdaijl)*4.0*kronecker(jtype,1));
      fi[0] -= tmp2*(rij[0]/rijmag);
      fi[1] -= tmp2*(rij[1]/rijmag);
      fi[2] -= tmp2*(rij[2]/rijmag);
      fj[0] -= tmp2*((-rjl[0]/rjlmag)-(rij[0]/rijmag));
      fj[1] -= tmp2*((-rjl[1]/rjlmag)-(rij[1]/rijmag));
      fj[2] -= tmp2*((-rjl[2]/rjlmag)-(rij[2]/rijmag));
      fl[0] -= tmp2*(rjl[0]/rjlmag);
      fl[1] -= tmp2*(rjl[1]/rjlmag);
      fl[2] -= tmp2*(rjl[2]/rjlmag);

      // coordination forces

      // dwik forces

      tmp2 = VA*.5*(tmp*dwjl*g*exp(lamdaijl))/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // PIJ forces

      tmp2 = VA*.5*(tmp*dN2[ltype]*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];

      // dgdN forces

      tmp2 = VA*.5*(tmp*tmp3*dwjl)/rjlmag;
      fj[0] -= tmp2*rjl[0];
      fj[1] -= tmp2*rjl[1];
      fj[2] -= tmp2*rjl[2];
      fl[0] += tmp2*rjl[0];
      fl[1] += tmp2*rjl[1];
      fl[2] += tmp2*rjl[2];
#pragma omp critical
	{
      f[atomi].x += fi[0]; f[atomi].y += fi[1]; f[atomi].z += fi[2];
      f[atomj].x += fj[0]; f[atomj].y += fj[1]; f[atomj].z += fj[2];
      f[atoml].x += fl[0]; f[atoml].y += fl[1]; f[atoml].z += fl[2];
	}
      /*
      if (vflag_atom) {
        rlj[0] = -rjl[0]; rlj[1] = -rjl[1]; rlj[2] = -rjl[2];
        v_tally3(atomi,atomj,atoml,fi,fl,rij,rlj);
	}
      */
    }
  }

  // evaluate Nij conj

  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3);

  // piRC forces
  
  //REBO_neighs_i = REBO_firstneigh[i];
  for (k = 0; k < NBond[atomi]; k++) {
    atomk = Bonded[atomi*MaxNBond + k];
    if (atomk !=atomj) {
      ktype = particle[atomk].type;
      
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
        (wik*kronecker(itype,1));
      SpN = Sp(Nki,Nmin,Nmax,dNki);

      tmp2 = VA*dN3[0]*dwik/rikmag;
      #pragma omp critical
	{
      f[atomi].x -= tmp2*rik[0];
      f[atomi].y -= tmp2*rik[1];
      f[atomi].z -= tmp2*rik[2];
      f[atomk].x += tmp2*rik[0];
      f[atomk].y += tmp2*rik[1];
      f[atomk].z += tmp2*rik[2];
	}
      //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

	tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)/rikmag;
	#pragma omp critical
	{
      f[atomi].x -= tmp2*rik[0];
      f[atomi].y -= tmp2*rik[1];
      f[atomi].z -= tmp2*rik[2];
      f[atomk].x += tmp2*rik[0];
      f[atomk].y += tmp2*rik[1];
      f[atomk].z += tmp2*rik[2];
	}
      //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);
      atomk = Bonded[atomi*MaxNBond + k];
      if (fabs(dNki) > TOL) {
        //REBO_neighs_k = REBO_firstneigh[atomk];
        for (n = 0; n < NBond[atomk]; n++) {
          atomn = Bonded[atomk*MaxNBond + n];
          if (atomn != atomi) {
            ntype = particle[atomn].type;
	    rkn[0] = R_PBC(particle,atomk,atomn,0,L.x);
	    rkn[1] = R_PBC(particle,atomk,atomn,1,L.y);
	    rkn[2] = R_PBC(particle,atomk,atomn,2,L.z);
            rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
            Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

            tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)/rknmag;
	    #pragma omp critical
	{
            f[atomk].x -= tmp2*rkn[0];
            f[atomk].y -= tmp2*rkn[1];
            f[atomk].z -= tmp2*rkn[2];
            f[atomn].x += tmp2*rkn[0];
            f[atomn].y += tmp2*rkn[1];
            f[atomn].z += tmp2*rkn[2];
	}
            //if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
          }
        }
      }
    }
  }

  // piRC forces

  //REBO_neighs = REBO_firstneigh[atomj];
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml !=atomi) {
      ltype = particle[atoml].type;
     
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] -
        (wjl*kronecker(jtype,1));
      SpN = Sp(Nlj,Nmin,Nmax,dNlj);

      tmp2 = VA*dN3[1]*dwjl/rjlmag;
       #pragma omp critical
	{
      f[atomj].x -= tmp2*rjl[0];
      f[atomj].y -= tmp2*rjl[1];
      f[atomj].z -= tmp2*rjl[2];
      f[atoml].x += tmp2*rjl[0];
      f[atoml].y += tmp2*rjl[1];
      f[atoml].z += tmp2*rjl[2];
	}
      //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)/rjlmag;
      #pragma omp critical
	{
      f[atomj].x -= tmp2*rjl[0];
      f[atomj].y -= tmp2*rjl[1];
      f[atomj].z -= tmp2*rjl[2];
      f[atoml].x += tmp2*rjl[0];
      f[atoml].y += tmp2*rjl[1];
      f[atoml].z += tmp2*rjl[2];
	}
      //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

      if (fabs(dNlj) > TOL) {
        //REBO_neighs_l = REBO_firstneigh[atoml];
        for (n = 0; n <NBond[atoml]; n++) {
	  atomn = Bonded[atoml*MaxNBond + n];
          //atomn = REBO_neighs_l[n];
          if (atomn != atomj) {
            ntype = particle[atomn].type;
	    rln[0] = R_PBC(particle,atoml,atomn,0,L.x);
	    rln[1] = R_PBC(particle,atoml,atomn,1,L.y);
	    rln[2] = R_PBC(particle,atoml,atomn,2,L.z);
            rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
            Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);

            tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)/rlnmag;
	    #pragma omp critical
	{
            f[atoml].x -= tmp2*rln[0];
            f[atoml].y -= tmp2*rln[1];
            f[atoml].z -= tmp2*rln[2];
            f[atomn].x += tmp2*rln[0];
            f[atomn].y += tmp2*rln[1];
            f[atomn].z += tmp2*rln[2];
	}
            //if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
          }
        }
      }
    }
  }

  Tij = 0.0;
  dN3[0] = 0.0;
  dN3[1] = 0.0;
  dN3[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3);
  Etmp = 0.0;

  if (fabs(Tij) > TOL) {
    atom2 = atomi;
    atom3 = atomj;
    r32[0] = R_PBC(particle,atom3,atom2,0,L.x);
    r32[1] = R_PBC(particle,atom3,atom2,1,L.y);
    r32[2] = R_PBC(particle,atom3,atom2,2,L.z);
    r32mag = sqrt((r32[0]*r32[0])+(r32[1]*r32[1])+(r32[2]*r32[2]));
    r23[0] = -r32[0];
    r23[1] = -r32[1];
    r23[2] = -r32[2];
    r23mag = r32mag;
    //REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < NBond[atomi]; k++) {
      atomk = Bonded[atomi*MaxNBond + k];
      atom1 = atomk;
      ktype = particle[atomk].type;
      if (atomk != atomj) {
	r21[0] = R_PBC(particle,atom2,atom1,0,L.x);
	r21[1] = R_PBC(particle,atom2,atom1,1,L.y);
	r21[2] = R_PBC(particle,atom2,atom1,2,L.z);
        r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
        cos321 = -1.0*((r21[0]*r32[0])+(r21[1]*r32[1])+(r21[2]*r32[2])) /
          (r21mag*r32mag);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);
        Sp2(cos321,thmin,thmax,dcut321);
        sin321 = sqrt(1.0 - cos321*cos321);
        sink2i = 1.0/(sin321*sin321);
        rik2i = 1.0/(r21mag*r21mag);
        if (sin321 != 0.0) {
          rr = (r23mag*r23mag)-(r21mag*r21mag);
          //rjk[0] = r21[0]-r23[0];
          //rjk[1] = r21[1]-r23[1];
          //rjk[2] = r21[2]-r23[2];

	  rjk[0] = R_PBC(particle,atomj,atomk,0,L.x);
          rjk[1] = R_PBC(particle,atomj,atomk,1,L.y);
          rjk[2] = R_PBC(particle,atomj,atomk,2,L.z);
	  
          rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
          rijrik = 2.0*r23mag*r21mag;
          rik2 = r21mag*r21mag;
          dctik = (-rr+rjk2)/(rijrik*rik2);
          dctij = (rr+rjk2)/(rijrik*r23mag*r23mag);
          dctjk = -2.0/rijrik;
          w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
          rijmag = r32mag;
          rikmag = r21mag;
          rij2 = r32mag*r32mag;
          rik2 = r21mag*r21mag;
          costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
          tspjik = Sp2(costmp,thmin,thmax,dtsjik);
          dtsjik = -dtsjik;

          //REBO_neighs_j = REBO_firstneigh[j];
          for (l = 0; l < NBond[atomj]; l++) {
            //atoml = REBO_neighs_j[l];
	    atoml = Bonded[atomj*MaxNBond + l];
            atom4 = atoml;
            ltype = particle[atoml].type;
            if (!(atoml == atomi || atoml == atomk)) {
	      r34[0] = R_PBC(particle,atom3,atom4,0,L.x);
	      r34[1] = R_PBC(particle,atom3,atom4,1,L.y);
	      r34[2] = R_PBC(particle,atom3,atom4,2,L.z);
              r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
              cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) /
                (r32mag*r34mag);
              cos234 = MIN(cos234,1.0);
              cos234 = MAX(cos234,-1.0);
              sin234 = sqrt(1.0 - cos234*cos234);
              sinl2i = 1.0/(sin234*sin234);
              rjl2i = 1.0/(r34mag*r34mag);

              if (sin234 != 0.0) {
                w34 = Sp(r34mag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dw34);
                rr = (r23mag*r23mag)-(r34mag*r34mag);
                //ril[0] = r23[0]+r34[0];
                //ril[1] = r23[1]+r34[1];
                //ril[2] = r23[2]+r34[2];

		ril[0] = R_PBC(particle,atomi,atoml,0,L.x);
                ril[1] = R_PBC(particle,atomi,atoml,1,L.y);
                ril[2] = R_PBC(particle,atomi,atoml,2,L.z);
		
                ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
                rijrjl = 2.0*r23mag*r34mag;
                rjl2 = r34mag*r34mag;
                dctjl = (-rr+ril2)/(rijrjl*rjl2);
                dctji = (rr+ril2)/(rijrjl*r23mag*r23mag);
                dctil = -2.0/rijrjl;
                rjlmag = r34mag;
                rjl2 = r34mag*r34mag;
                costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
                tspijl = Sp2(costmp,thmin,thmax,dtsijl);
                dtsijl = -dtsijl;
                prefactor = VA*Tij;

                cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
                cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
                cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
                cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
                cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
                cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

                cwnum = (cross321[0]*cross234[0]) +
                  (cross321[1]*cross234[1]) + (cross321[2]*cross234[2]);
                cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
                om1234 = cwnum/cwnom;
                cw = om1234;
                Etmp += ((1.0-square(om1234))*w21*w34) *
                  (1.0-tspjik)*(1.0-tspijl);

                dt1dik = (rik2i)-(dctik*sink2i*cos321);
                dt1djk = (-dctjk*sink2i*cos321);
                dt1djl = (rjl2i)-(dctjl*sinl2i*cos234);
                dt1dil = (-dctil*sinl2i*cos234);
                dt1dij = (2.0/(r23mag*r23mag))-(dctij*sink2i*cos321) -
                  (dctji*sinl2i*cos234);

                dt2dik[0] = (-r23[2]*cross234[1])+(r23[1]*cross234[2]);
                dt2dik[1] = (-r23[0]*cross234[2])+(r23[2]*cross234[0]);
                dt2dik[2] = (-r23[1]*cross234[0])+(r23[0]*cross234[1]);

                dt2djl[0] = (-r23[1]*cross321[2])+(r23[2]*cross321[1]);
                dt2djl[1] = (-r23[2]*cross321[0])+(r23[0]*cross321[2]);
                dt2djl[2] = (-r23[0]*cross321[1])+(r23[1]*cross321[0]);

                dt2dij[0] = (r21[2]*cross234[1])-(r34[2]*cross321[1]) -
                  (r21[1]*cross234[2])+(r34[1]*cross321[2]);
                dt2dij[1] = (r21[0]*cross234[2])-(r34[0]*cross321[2]) -
                  (r21[2]*cross234[0])+(r34[2]*cross321[0]);
                dt2dij[2] = (r21[1]*cross234[0])-(r34[1]*cross321[0]) -
                  (r21[0]*cross234[1])+(r34[0]*cross321[1]);

                aa = (prefactor*2.0*cw/cwnom)*w21*w34 *
                  (1.0-tspjik)*(1.0-tspijl);
                aaa1 = -prefactor*(1.0-square(om1234)) *
                  (1.0-tspjik)*(1.0-tspijl);
                aaa2 = aaa1*w21*w34;
                at2 = aa*cwnum;

                fcijpc = (-dt1dij*at2)+(aaa2*dtsjik*dctij*(1.0-tspijl)) +
                  (aaa2*dtsijl*dctji*(1.0-tspjik));
                fcikpc = (-dt1dik*at2)+(aaa2*dtsjik*dctik*(1.0-tspijl));
                fcjlpc = (-dt1djl*at2)+(aaa2*dtsijl*dctjl*(1.0-tspjik));
                fcjkpc = (-dt1djk*at2)+(aaa2*dtsjik*dctjk*(1.0-tspijl));
                fcilpc = (-dt1dil*at2)+(aaa2*dtsijl*dctil*(1.0-tspjik));

                F23[0] = (fcijpc*r23[0])+(aa*dt2dij[0]);
                F23[1] = (fcijpc*r23[1])+(aa*dt2dij[1]);
                F23[2] = (fcijpc*r23[2])+(aa*dt2dij[2]);

                F12[0] = (fcikpc*r21[0])+(aa*dt2dik[0]);
                F12[1] = (fcikpc*r21[1])+(aa*dt2dik[1]);
                F12[2] = (fcikpc*r21[2])+(aa*dt2dik[2]);

                F34[0] = (fcjlpc*r34[0])+(aa*dt2djl[0]);
                F34[1] = (fcjlpc*r34[1])+(aa*dt2djl[1]);
                F34[2] = (fcjlpc*r34[2])+(aa*dt2djl[2]);

                F31[0] = (fcjkpc*rjk[0]);
                F31[1] = (fcjkpc*rjk[1]);
                F31[2] = (fcjkpc*rjk[2]);

                F24[0] = (fcilpc*ril[0]);
                F24[1] = (fcilpc*ril[1]);
                F24[2] = (fcilpc*ril[2]);

                f1[0] = -F12[0]-F31[0];
                f1[1] = -F12[1]-F31[1];
                f1[2] = -F12[2]-F31[2];
                f2[0] = F23[0]+F12[0]+F24[0];
                f2[1] = F23[1]+F12[1]+F24[1];
                f2[2] = F23[2]+F12[2]+F24[2];
                f3[0] = -F23[0]+F34[0]+F31[0];
                f3[1] = -F23[1]+F34[1]+F31[1];
                f3[2] = -F23[2]+F34[2]+F31[2];
                f4[0] = -F34[0]-F24[0];
                f4[1] = -F34[1]-F24[1];
                f4[2] = -F34[2]-F24[2];

                // coordination forces

                tmp2 = VA*Tij*((1.0-(om1234*om1234))) *
                  (1.0-tspjik)*(1.0-tspijl)*dw21*w34/r21mag;
                f2[0] -= tmp2*r21[0];
                f2[1] -= tmp2*r21[1];
                f2[2] -= tmp2*r21[2];
                f1[0] += tmp2*r21[0];
                f1[1] += tmp2*r21[1];
                f1[2] += tmp2*r21[2];

                tmp2 = VA*Tij*((1.0-(om1234*om1234))) *
                  (1.0-tspjik)*(1.0-tspijl)*w21*dw34/r34mag;
                f3[0] -= tmp2*r34[0];
                f3[1] -= tmp2*r34[1];
                f3[2] -= tmp2*r34[2];
                f4[0] += tmp2*r34[0];
                f4[1] += tmp2*r34[1];
                f4[2] += tmp2*r34[2];
#pragma omp critical
	{
                f[atom1].x += f1[0]; f[atom1].y += f1[1];
                f[atom1].z += f1[2];
                f[atom2].x += f2[0]; f[atom2].y += f2[1];
                f[atom2].z += f2[2];
                f[atom3].x += f3[0]; f[atom3].y += f3[1];
                f[atom3].z += f3[2];
                f[atom4].x += f4[0]; f[atom4].y += f4[1];
                f[atom4].z += f4[2];
	}
		/*
                if (vflag_atom) {
                  r13[0] = -rjk[0]; r13[1] = -rjk[1]; r13[2] = -rjk[2];
                  r43[0] = -r34[0]; r43[1] = -r34[1]; r43[2] = -r34[2];
                  v_tally4(atom1,atom2,atom3,atom4,f1,f2,f4,r13,r23,r43);
		  }
		*/
              }
            }
          }
        }
      }
    }

    // Tij forces now that we have Etmp

    //REBO_neighs = REBO_firstneigh[i];
    for (k = 0; k < NBond[atomi]; k++) {
      atomk = Bonded[atomi*MaxNBond + k];
      //atomk = REBO_neighs[k];
      if (atomk != atomj) {
        ktype = particle[atomk].type;
	rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
	rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
	rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
        rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
        wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dwik);
        Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
          (wik*kronecker(itype,1));
        SpN = Sp(Nki,Nmin,Nmax,dNki);

        tmp2 = VA*dN3[0]*dwik*Etmp/rikmag;
	#pragma omp critical
	{
        f[atomi].x -= tmp2*rik[0];
        f[atomi].y -= tmp2*rik[1];
        f[atomi].z -= tmp2*rik[2];
        f[atomk].x += tmp2*rik[0];
        f[atomk].y += tmp2*rik[1];
        f[atomk].z += tmp2*rik[2];
	}
        //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

        tmp2 = VA*dN3[2]*(2.0*NconjtmpI*dwik*SpN)*Etmp/rikmag;
	#pragma omp critical
	{
        f[atomi].x -= tmp2*rik[0];
        f[atomi].y -= tmp2*rik[1];
        f[atomi].z -= tmp2*rik[2];
        f[atomk].x += tmp2*rik[0];
        f[atomk].y += tmp2*rik[1];
        f[atomk].z += tmp2*rik[2];
	}
        //if (vflag_atom) v_tally2(atomi,atomk,-tmp2,rik);

        if (fabs(dNki) > TOL) {
          //REBO_neighs_k = REBO_firstneigh[atomk];
          for (n = 0; n < NBond[atomk]; n++) {
            atomn = Bonded[atomk*MaxNBond + n];
            ntype = particle[atomn].type;
            if (atomn != atomi) {
	      rkn[0] = R_PBC(particle,atomk,atomn,0,L.x);
	      rkn[1] = R_PBC(particle,atomk,atomn,1,L.y);
	      rkn[2] = R_PBC(particle,atomk,atomn,2,L.z);
              rknmag = sqrt((rkn[0]*rkn[0])+(rkn[1]*rkn[1])+(rkn[2]*rkn[2]));
              Sp(rknmag,rcmin[ktype][ntype],rcmax[ktype][ntype],dwkn);

              tmp2 = VA*dN3[2]*(2.0*NconjtmpI*wik*dNki*dwkn)*Etmp/rknmag;
	      #pragma omp critical
	{
              f[atomk].x -= tmp2*rkn[0];
              f[atomk].y -= tmp2*rkn[1];
              f[atomk].z -= tmp2*rkn[2];
              f[atomn].x += tmp2*rkn[0];
              f[atomn].y += tmp2*rkn[1];
              f[atomn].z += tmp2*rkn[2];
	}
              //if (vflag_atom) v_tally2(atomk,atomn,-tmp2,rkn);
            }
          }
        }
      }
    }

    // Tij forces

    //REBO_neighs = REBO_firstneigh[j];
    for (l = 0; l < NBond[atomj]; l++) {
      atoml = Bonded[atomj*MaxNBond + l];
      if (atoml != atomi) {
        ltype = particle[atoml].type;
      	rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
	rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
	rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
        rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
        wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dwjl);
        Nlj = nC[atoml]-(wjl*kronecker(jtype,0))+nH[atoml] -
          (wjl*kronecker(jtype,1));
        SpN = Sp(Nlj,Nmin,Nmax,dNlj);

        tmp2 = VA*dN3[1]*dwjl*Etmp/rjlmag;
	#pragma omp critical
	{
        f[atomj].x -= tmp2*rjl[0];
        f[atomj].y -= tmp2*rjl[1];
        f[atomj].z -= tmp2*rjl[2];
        f[atoml].x += tmp2*rjl[0];
        f[atoml].y += tmp2*rjl[1];
        f[atoml].z += tmp2*rjl[2];
	}
        //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

        tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*dwjl*SpN)*Etmp/rjlmag;
	#pragma omp critical
	{
        f[atomj].x -= tmp2*rjl[0];
        f[atomj].y -= tmp2*rjl[1];
        f[atomj].z -= tmp2*rjl[2];
        f[atoml].x += tmp2*rjl[0];
        f[atoml].y += tmp2*rjl[1];
        f[atoml].z += tmp2*rjl[2];
	}
        //if (vflag_atom) v_tally2(atomj,atoml,-tmp2,rjl);

        if (fabs(dNlj) > TOL) {
          //REBO_neighs_l = REBO_firstneigh[atoml];
          for (n = 0; n < NBond[atoml]; n++) {
	    atomn = Bonded[atoml*MaxNBond + n];
            //atomn = REBO_neighs_l[n];
            ntype = particle[atomn].type;
            if (atomn !=atomj) {
	      rln[0] = R_PBC(particle,atoml,atomn,0,L.x);
	      rln[1] = R_PBC(particle,atoml,atomn,1,L.y);
	      rln[2] = R_PBC(particle,atoml,atomn,2,L.z);
              rlnmag = sqrt((rln[0]*rln[0])+(rln[1]*rln[1])+(rln[2]*rln[2]));
              Sp(rlnmag,rcmin[ltype][ntype],rcmax[ltype][ntype],dwln);
	      
              tmp2 = VA*dN3[2]*(2.0*NconjtmpJ*wjl*dNlj*dwln)*Etmp/rlnmag;
	      #pragma omp critical
	{
	      f[atoml].x -= tmp2*rln[0];
              f[atoml].y -= tmp2*rln[1];
              f[atoml].z -= tmp2*rln[2];
              f[atomn].x += tmp2*rln[0];
              f[atomn].y += tmp2*rln[1];
              f[atomn].z += tmp2*rln[2];
	}
              //if (vflag_atom) v_tally2(atoml,atomn,-tmp2,rln);
            }
          }
        }
      }
    }
  }

  bij = (0.5*(pij+pji))+piRC+(Tij*Etmp);
  //printf("bij=%.16f pij=%.16f pji=%.16f piRC=%.16f Tij=%.16f Etmp=%.16f\n",bij,pij,pji,piRC,Tij,Etmp);
  return bij;


 }


double gSpline(double costh, double Nij, int typei, double *dgdc, double *dgdN)
{
  double coeffs[6],dS,g1,g2,dg1,dg2,cut,g;
  int i,j;

  i = 0;
  j = 0;
  g = 0.0;
  cut = 0.0;
  dS = 0.0;
  dg1 = 0.0;
    dg2 = 0.0;
   *dgdc = 0.0;
   *dgdN = 0.0;

  // central atom is Carbon

   if (typei == 0) {
    if (costh < gCdom[0]) costh = gCdom[0]; //gCdom[i] parameter
    if (costh > gCdom[4]) costh = gCdom[4];
    if (Nij >= NCmax) {
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1])  {
          for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j]; //gC2[i] parameter
         } 
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      g = g2;
      *dgdc = dg2;
      *dgdN = 0.0;
    }
    if (Nij <= NCmin) {                        //NCmin parameter
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
          for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
        }
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      g = g1;
      *dgdc = dg1;
      *dgdN = 0.0;
    }
    if (Nij > NCmin && Nij < NCmax) { //NCmax parameter
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
          for (j = 0; j < 6; j++) coeffs[j] = gC1[i][j];
        }
      }
      g1 = Sp5th(costh,coeffs,&dg1);
      for (i = 0; i < 4; i++) {
        if (costh >= gCdom[i] && costh <= gCdom[i+1]) {
          for (j = 0; j < 6; j++) coeffs[j] = gC2[i][j];
        }
      }
      g2 = Sp5th(costh,coeffs,&dg2);
      cut = Sp(Nij,NCmin,NCmax,dS);
      g = g2+cut*(g1-g2);
      *dgdc = dg2+(cut*(dg1-dg2));
      *dgdN = dS*(g1-g2);
    }
  }

  // central atom is Hydrogen

  if (typei == 1) {
    if (costh < gHdom[0]) costh = gHdom[0]; //gHdom parameter
    if (costh > gHdom[3]) costh = gHdom[3];
    for (i = 0; i < 3; i++) {
      if (costh >= gHdom[i] && costh <= gHdom[i+1]) {
        for (j = 0; j < 6; j++) coeffs[j] = gH[i][j];
      }
    }
    g = Sp5th(costh,coeffs,&dg1);
    *dgdN = 0.0;
    *dgdc = dg1;
  }

  return g;
}

double PijSpline(double NijC, double NijH, int typei, int typej,double dN2[2])
{
  int x,y,i,done;
  double Pij,coeffs[16];
  
  for (i = 0; i < 16; i++) coeffs[i]=0.0;
  
  x = 0;
  y = 0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  done = 0;

  // if inputs are out of bounds set them back to a point in bounds

  if (typei == 0 && typej == 0) {
    if (NijC < pCCdom[0][0]) NijC=pCCdom[0][0];
    if (NijC > pCCdom[0][1]) NijC=pCCdom[0][1];
    if (NijH < pCCdom[1][0]) NijH=pCCdom[1][0];
    if (NijH > pCCdom[1][1]) NijH=pCCdom[1][1];
    //printf("NijC=%.16f NijH=%.16f\n",NijC,NijH);
    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCCf[(int) NijC][(int) NijH];
      dN2[0] = PCCdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCCdfdy[(int) NijC][(int) NijH];
      done = 1;
      //printf("Pij=%.16f NijC=%.16f NijH=%.16f dN[0]=%.16f dN[1]=%.16f \n",Pij,NijC,NijH,dN2[0],dN2[1]);
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCC[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
      //printf("Pij=%.16f NijC=%.16f NijH=%.16f coeffs=%.16f dN2[0]==%.16f dN2[1]=%.16f \n",Pij,NijC,NijH,coeffs,dN2[0],dN2[1]);
    }
  }

  // if inputs are out of bounds set them back to a point in bounds

   if (typei == 0 && typej == 1){
     if (NijC < pCHdom[0][0]) NijC=pCHdom[0][0];
     if (NijC > pCHdom[0][1]) NijC=pCHdom[0][1];
      if (NijH < pCHdom[1][0]) NijH=pCHdom[1][0];
      if (NijH > pCHdom[1][1]) NijH=pCHdom[1][1];

    if (fabs(NijC-floor(NijC)) < TOL && fabs(NijH-floor(NijH)) < TOL) {
      Pij = PCHf[(int) NijC][(int) NijH];
      dN2[0] = PCHdfdx[(int) NijC][(int) NijH];
      dN2[1] = PCHdfdy[(int) NijC][(int) NijH];
      done = 1;
    }
    if (done == 0) {
      x = (int) (floor(NijC));
      y = (int) (floor(NijH));
      for (i = 0; i<16; i++) coeffs[i] = pCH[x][y][i];
      Pij = Spbicubic(NijC,NijH,coeffs,dN2);
    }
  }

  if (typei == 1 && typej == 0) {
    Pij = 0.0;
    dN2[0] = 0.0;
    dN2[1] = 0.0;
  }


  if (typei == 1 && typej == 1) {
    Pij = 0.0;
    dN2[0] = 0.0;
    dN2[1] = 0.0;
  }
  return Pij;
}
double Sp5th(double x, double coeffs[6], double *df)
{
  double f, d;
  const double x2 = x*x;
  const double x3 = x2*x;
  
  f  = coeffs[0];
  f += coeffs[1]*x;
  d  = coeffs[1];
  f += coeffs[2]*x2;
  d += 2.0*coeffs[2]*x;
  f += coeffs[3]*x3;
  d += 3.0*coeffs[3]*x2;
  f += coeffs[4]*x2*x2;
  d += 4.0*coeffs[4]*x3;
  f += coeffs[5]*x2*x3;
  d += 5.0*coeffs[5]*x2*x2;
  
  *df = d;
  return f;
}
double Spbicubic(double x, double y,double coeffs[16], double df[2])
{
  double f,xn,yn,xn1,yn1,c;
  int i,j;

  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;

  xn = 1.0;
  for (i = 0; i < 4; i++) {
    yn = 1.0;
    for (j = 0; j < 4; j++) {
      c = coeffs[i*4+j];

      f += c*xn*yn;
      if (i > 0) df[0] += c * ((double) i) * xn1 * yn;
      if (j > 0) df[1] += c * ((double) j) * xn * yn1;

      yn1 = yn;
      yn *= y;
    }
    xn1 = xn;
    xn *= x;
  }

  return f;
}

double piRCSpline(double Nij, double Nji, double Nijconj,int typei, int typej, double dN3[3])
{
  int x,y,z,i,done;
  double piRC,coeffs[64];
  x=0;
  y=0;
  z=0;
  i=0;

  done=0;

  for (i=0; i<64; i++) coeffs[i]=0.0;

  if (typei==0 && typej==0) {
    //if the inputs are out of bounds set them back to a point in bounds
    if (Nij<piCCdom[0][0]) Nij=piCCdom[0][0];
    if (Nij>piCCdom[0][1]) Nij=piCCdom[0][1];
    if (Nji<piCCdom[1][0]) Nji=piCCdom[1][0];
    if (Nji>piCCdom[1][1]) Nji=piCCdom[1][1];
    if (Nijconj<piCCdom[2][0]) Nijconj=piCCdom[2][0];
    if (Nijconj>piCCdom[2][1]) Nijconj=piCCdom[2][1];

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
        fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCCf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCCdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCCdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCCdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCCdom[0][1]; i++)
        if (Nij>=(double) i && Nij<=(double) i+1 || Nij==(double) i) x=i;
      for (i=0; i<piCCdom[1][1]; i++)
        if (Nji>=(double) i && Nji<=(double) i+1 || Nji==(double) i) y=i;
      for (i=0; i<piCCdom[2][1]; i++)
        if (Nijconj>=(double) i && Nijconj<=(double) i+1 ||
            Nijconj==(double) i) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCC[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }


  // CH interaction

  if (typei==0 && typej==1 || typei==1 && typej==0) {
    // if the inputs are out of bounds set them back to a point in bounds

    if (Nij<piCHdom[0][0] || Nij>piCHdom[0][1] ||
        Nji<piCHdom[1][0] || Nji>piCHdom[1][1] ||
        Nijconj<piCHdom[2][0] || Nijconj>piCHdom[2][1]) {
      if (Nij<piCHdom[0][0]) Nij=piCHdom[0][0];
      if (Nij>piCHdom[0][1]) Nij=piCHdom[0][1];
      if (Nji<piCHdom[1][0]) Nji=piCHdom[1][0];
      if (Nji>piCHdom[1][1]) Nji=piCHdom[1][1];
      if (Nijconj<piCHdom[2][0]) Nijconj=piCHdom[2][0];
      if (Nijconj>piCHdom[2][1]) Nijconj=piCHdom[2][1];
    }

    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
        fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piCHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piCHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piCHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piCHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }

    if (done==0) {
      for (i=0; i<piCHdom[0][1]; i++)
        if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piCHdom[1][1]; i++)
        if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piCHdom[2][1]; i++)
        if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piCH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }

  if (typei==1 && typej==1) {
    if (Nij<piHHdom[0][0] || Nij>piHHdom[0][1] ||
        Nji<piHHdom[1][0] || Nji>piHHdom[1][1] ||
        Nijconj<piHHdom[2][0] || Nijconj>piHHdom[2][1]) {
      Nij=0.0;
      Nji=0.0;
      Nijconj=0.0;
    }
    if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
        fabs(Nijconj-floor(Nijconj))<TOL) {
      piRC=piHHf[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[0]=piHHdfdx[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[1]=piHHdfdy[(int) Nij][(int) Nji][(int) Nijconj];
      dN3[2]=piHHdfdz[(int) Nij][(int) Nji][(int) Nijconj];
      done=1;
    }
    if (done==0) {
      for (i=0; i<piHHdom[0][1]; i++)
        if (Nij>=i && Nij<=i+1) x=i;
      for (i=0; i<piHHdom[1][1]; i++)
        if (Nji>=i && Nji<=i+1) y=i;
      for (i=0; i<piHHdom[2][1]; i++)
        if (Nijconj>=i && Nijconj<=i+1) z=i;

      for (i=0; i<64; i++) coeffs[i]=piHH[x][y][z][i];
      piRC=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
    }
  }

  return piRC;
}

double TijSpline(double Nij, double Nji,double Nijconj, double dN3[3])
{
  int x,y,z,i,done;
  double Tijf,coeffs[64];

  x=0;
  y=0;
  z=0;
  i=0;
  Tijf=0.0;
  done=0;
  for (i=0; i<64; i++) coeffs[i]=0.0;

  //if the inputs are out of bounds set them back to a point in bounds

  if (Nij<Tijdom[0][0]) Nij=Tijdom[0][0];
  if (Nij>Tijdom[0][1]) Nij=Tijdom[0][1];
  if (Nji<Tijdom[1][0]) Nji=Tijdom[1][0];
  if (Nji>Tijdom[1][1]) Nji=Tijdom[1][1];
  if (Nijconj<Tijdom[2][0]) Nijconj=Tijdom[2][0];
  if (Nijconj>Tijdom[2][1]) Nijconj=Tijdom[2][1];

  if (fabs(Nij-floor(Nij))<TOL && fabs(Nji-floor(Nji))<TOL &&
      fabs(Nijconj-floor(Nijconj))<TOL) {
    Tijf=Tf[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[0]=Tdfdx[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[1]=Tdfdy[(int) Nij][(int) Nji][(int) Nijconj];
    dN3[2]=Tdfdz[(int) Nij][(int) Nji][(int) Nijconj];
    done=1;
  }

  if (done==0) {
    for (i=0; i<Tijdom[0][1]; i++)
      if (Nij>=i && Nij<=i+1) x=i;
    for (i=0; i<Tijdom[1][1]; i++)
      if (Nji>=i && Nji<=i+1) y=i;
    for (i=0; i<Tijdom[2][1]; i++)
      if (Nijconj>=i && Nijconj<=i+1) z=i;

    for (i=0; i<64; i++) coeffs[i]=Tijc[x][y][z][i];
    Tijf=Sptricubic(Nij,Nji,Nijconj,coeffs,dN3);
  }

  return Tijf;
}

double Sptricubic(double x, double y, double z,
                              double coeffs[64], double df[3])
{
  double f,ir,jr,kr,xn,yn,zn,xn1,yn1,zn1,c;
  int i,j,k;

  f = 0.0;
  df[0] = 0.0;
  df[1] = 0.0;
  df[2] = 0.0;

  xn = 1.0;
  for (i = 0; i < 4; i++) {
    ir = (double) i;
    yn = 1.0;
    for (j = 0; j < 4; j++) {
      jr = (double) j;
      zn = 1.0;
      for (k = 0; k < 4; k++) {
        kr = (double) k;
        c = coeffs[16*i+4*j+k];
        f += c*xn*yn*zn;
        if (i > 0) df[0] += c * ir * xn1 * yn * zn;
        if (j > 0) df[1] += c * jr * xn * yn1 * zn;
        if (k > 0) df[2] += c * kr * xn * yn * zn1;
        zn1 = zn;
        zn *= z;
      }
      yn1 = yn;
      yn *= y;
    }
    xn1 = xn;
    xn *= x;
  }

  return f;
}

/*
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
*/
void calc_force_arrays(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *angle, double *Fc,int *Neighb,int *NNeighb,double *rijShiled_1,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once)
{
  Check_List(particle,R0,Natoms,List,NList,L);
  calc_rij(NBond,Bonded,BondNum,Natoms,particle,L,rij,Fc,BondParams,Neighb,NNeighb,rijShiled_1,NList,List,NBond_once,Bonded_once);
  //calc_angles(angle,Natoms,NBond,Bonded,rij,particle,L);
}

/****************************************************************/
/* Calculate the total force array                              */

void Calc_Force(point *Fterm, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double &E_tot,int *NBond_once,int *Bonded_once)
{
  int atomi;
  double E_TORSION,E_KC,E_REBO;
  
  clock_t start , diff;
  int msec;
  
  dr_ij *dN = new dr_ij[1];
 
  E_TORSION = 0.0; 
  E_REBO    = 0.0;
  E_KC      = 0.0;
  E_tot     = 0.0;
  // Calculate the bond force term
#pragma omp parallel for
  for(atomi=0 ; atomi < Natoms ; atomi++){
    Fterm[atomi].x=0.0;
    Fterm[atomi].y=0.0;
    Fterm[atomi].z=0.0;
  }
  
  //start = clock();
  
  calc_dE_REBO(Fterm,angle,BondNum,Natoms,L,BondParams,Bonded,particle,Fc,NBond,E_REBO,NBond_once,Bonded_once);
  
  //diff = clock() - start;
  //msec = diff * 1000 / CLOCKS_PER_SEC;
  //printf("dE_REBO Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
  
  if (Coulomb_flag==1 || Coulomb_flag==2)
    TORSION_(particle,Fterm,BondNum,Natoms,L,Bonded,NBond,E_TORSION,NBond_once,Bonded_once);
  
  //diff = clock() - start;
  //msec = diff * 1000 / CLOCKS_PER_SEC;
  
  //printf("Torsion Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
  
  //start = clock();
  
  if (Coulomb_flag==2)
    calc_F_Rep(Fterm,Natoms,particle,rijShiled_1,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,dN,E_KC);
  
  //diff = clock() - start;
  //msec = diff * 1000 / CLOCKS_PER_SEC;
  //printf("vdW Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);
  
  
#pragma omp parallel for
  for(atomi=0 ; atomi < Natoms ; atomi++){
    particle[atomi].Fx = Fterm[atomi].x;
    particle[atomi].Fy = Fterm[atomi].y;
    particle[atomi].Fz = Fterm[atomi].z;
  }
  
  E_tot = E_REBO + E_TORSION + E_KC;
  delete [] dN;
}


void read_file(char *filename)
{
  int i,j,k,l,limit;
  char s[MAXLINE];

  // REBO Parameters (AIREBO)

  double rcmin_CC,rcmin_CH,rcmin_HH,rcmax_CC,rcmax_CH,
    rcmax_HH,rcmaxp_CC,rcmaxp_CH,rcmaxp_HH;
  double Q_CC,Q_CH,Q_HH,alpha_CC,alpha_CH,alpha_HH,A_CC,A_CH,A_HH;
  double BIJc_CC1,BIJc_CC2,BIJc_CC3,BIJc_CH1,BIJc_CH2,BIJc_CH3,
    BIJc_HH1,BIJc_HH2,BIJc_HH3;
  double Beta_CC1,Beta_CC2,Beta_CC3,Beta_CH1,Beta_CH2,Beta_CH3,Beta_HH1,Beta_HH2,Beta_HH3;
  double rho_CC,rho_CH,rho_HH;

  // LJ Parameters (AIREBO)
  
  double rcLJmin_CC,rcLJmin_CH,rcLJmin_HH,rcLJmax_CC,rcLJmax_CH,
    rcLJmax_HH,bLJmin_CC;
  double bLJmin_CH,bLJmin_HH,bLJmax_CC,bLJmax_CH,bLJmax_HH,
    epsilon_CC,epsilon_CH,epsilon_HH;
  double sigma_CC,sigma_CH,sigma_HH,epsilonT_CCCC,epsilonT_CCCH,epsilonT_HCCH;
  
  // read file on proc 0
  
  
  //FILE *fp = open_potential(filename);
  FILE *fp;

  /*
  if (fp == NULL) {
    char str[128];
    sprintf(str,"Cannot open AIREBO potential file %s",filename);
    error->one(FLERR,str);
  }
  */
  if ((fp = fopen(filename,"r"))== NULL) {
    char str[128];
    sprintf(str,"Cannot open AIREBO potental file %s",filename);
    exit(0);
  }
  
  // skip initial comment lines
  
  while (1) {
    fgets(s,MAXLINE,fp);
    if (s[0] != '#') break;
  }
    
    // read parameters
    
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcmaxp_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&smin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Nmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmin);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&NCmax);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Q_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&alpha_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&A_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&BIJc_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CC3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_CH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH1);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH2);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&Beta_HH3);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rho_HH);
    //cerr<<"BIJc_CH1 ="<<BIJc_CH1<<endl;
    //cerr<<"rcmin_CC ="<<rcmin_CC<<endl;
    // LJ parameters

    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&rcLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmin_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&bLJmax_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilon_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_CH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&sigma_HH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCC);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_CCCH);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%lg",&epsilonT_HCCH);

    // gC spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    // number-1 = # of domains for the spline

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);
  
    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gCdom[i]);
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&gC1[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);
    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&gC2[i][j]);
      }
    }
  
    // gH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit; i++) {
      fgets(s,MAXLINE,fp);
      sscanf(s,"%lg",&gHdom[i]);
    }

    fgets(s,MAXLINE,fp);

    for (i = 0; i < limit-1; i++) {
      for (j = 0; j < 6; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&gH[i][j]);
      }
    }

    // pCC spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&pCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) pCCdom[0][1]; i++) {
      for (j = 0; j < (int) pCCdom[1][1]; j++) {
        for (k = 0; k < 16; k++) {
          fgets(s,MAXLINE,fp);
          sscanf(s,"%lg",&pCC[i][j][k]);
        }
      }
    }

    // pCH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/2; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&pCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) pCHdom[0][1]; i++) {
      for (j = 0; j < (int) pCHdom[1][1]; j++) {
        for (k = 0; k < 16; k++) {
          fgets(s,MAXLINE,fp);
          sscanf(s,"%lg",&pCH[i][j][k]);
        }
      }
    }

    // piCC cpline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&piCCdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) piCCdom[0][1]; i++) {
      for (j = 0; j < (int) piCCdom[1][1]; j++) {
        for (k = 0; k < (int) piCCdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&piCC[i][j][k][l]);
          }
        }
      }
    }
  
    // piCH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&piCHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) piCHdom[0][1]; i++) {
      for (j = 0; j < (int) piCHdom[1][1]; j++) {
        for (k = 0; k < (int) piCHdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&piCH[i][j][k][l]);
          }
        }
      }
    }

    // piHH spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&piHHdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) piHHdom[0][1]; i++) {
      for (j = 0; j < (int) piHHdom[1][1]; j++) {
        for (k = 0; k < (int) piHHdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&piHH[i][j][k][l]);
          }
        }
      }
    }

    // Tij spline

    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);
    fgets(s,MAXLINE,fp);

    fgets(s,MAXLINE,fp);
    sscanf(s,"%d",&limit);

    for (i = 0; i < limit/2; i++) {
      for (j = 0; j < limit/3; j++) {
        fgets(s,MAXLINE,fp);
        sscanf(s,"%lg",&Tijdom[i][j]);
      }
    }
    fgets(s,MAXLINE,fp);

    for (i = 0; i < (int) Tijdom[0][1]; i++) {
      for (j = 0; j < (int) Tijdom[1][1]; j++) {
        for (k = 0; k < (int) Tijdom[2][1]; k++) {
          for (l = 0; l < 64; l = l+1) {
            fgets(s,MAXLINE,fp);
            sscanf(s,"%lg",&Tijc[i][j][k][l]);
          }
        }
      }
    }

    fclose(fp);
  

  // store read-in values in arrays

  // REBO

    rcmin[0][0] = rcmin_CC;
    
    rcmin[0][1] = rcmin_CH;
    rcmin[1][0] = rcmin[0][1];
    rcmin[1][1] = rcmin_HH;

    rcmax[0][0] = rcmax_CC;
    rcmax[0][1] = rcmax_CH;
    rcmax[1][0] = rcmax[0][1];
    rcmax[1][1] = rcmax_HH;

    rcmaxsq[0][0] = rcmax[0][0]*rcmax[0][0];
    rcmaxsq[1][0] = rcmax[1][0]*rcmax[1][0];
    rcmaxsq[0][1] = rcmax[0][1]*rcmax[0][1];
    rcmaxsq[1][1] = rcmax[1][1]*rcmax[1][1];

    rcmaxp[0][0] = rcmaxp_CC;
    rcmaxp[0][1] = rcmaxp_CH;
    rcmaxp[1][0] = rcmaxp[0][1];
    rcmaxp[1][1] = rcmaxp_HH;

    Q[0][0] = Q_CC;
    Q[0][1] = Q_CH;
    Q[1][0] = Q[0][1];
    Q[1][1] = Q_HH;

    alpha[0][0] = alpha_CC;
    alpha[0][1] = alpha_CH;
    alpha[1][0] = alpha[0][1];
    alpha[1][1] = alpha_HH;

    A_ij[0][0] = A_CC;
    A_ij[0][1] = A_CH;
    A_ij[1][0] = A_ij[0][1];
    A_ij[1][1] = A_HH;

    rho[0][0] = rho_CC;
    rho[0][1] = rho_CH;
    rho[1][0] = rho[0][1];
    rho[1][1] = rho_HH;

    BIJc[0][0][0] = BIJc_CC1;
    BIJc[0][0][1] = BIJc_CC2;
    BIJc[0][0][2] = BIJc_CC3;
    BIJc[0][1][0] = BIJc_CH1;
    BIJc[0][1][1] = BIJc_CH2;
    BIJc[0][1][2] = BIJc_CH3;
    BIJc[1][0][0] = BIJc_CH1;
    BIJc[1][0][1] = BIJc_CH2;
    BIJc[1][0][2] = BIJc_CH3;
    BIJc[1][1][0] = BIJc_HH1;
    BIJc[1][1][1] = BIJc_HH2;
    BIJc[1][1][2] = BIJc_HH3;

    Beta[0][0][0] = Beta_CC1;
    Beta[0][0][1] = Beta_CC2;
    Beta[0][0][2] = Beta_CC3;
    Beta[0][1][0] = Beta_CH1;
    Beta[0][1][1] = Beta_CH2;
    Beta[0][1][2] = Beta_CH3;
    Beta[1][0][0] = Beta_CH1;
    Beta[1][0][1] = Beta_CH2;
    Beta[1][0][2] = Beta_CH3;
    Beta[1][1][0] = Beta_HH1;
    Beta[1][1][1] = Beta_HH2;
    Beta[1][1][2] = Beta_HH3;

    // LJ

    rcLJmin[0][0] = rcLJmin_CC;
    rcLJmin[0][1] = rcLJmin_CH;
    rcLJmin[1][0] = rcLJmin[0][1];
    rcLJmin[1][1] = rcLJmin_HH;

    rcLJmax[0][0] = rcLJmax_CC;
    rcLJmax[0][1] = rcLJmax_CH;
    rcLJmax[1][0] = rcLJmax[0][1];
    rcLJmax[1][1] = rcLJmax_HH;

    rcLJmaxsq[0][0] = rcLJmax[0][0]*rcLJmax[0][0];
    rcLJmaxsq[1][0] = rcLJmax[1][0]*rcLJmax[1][0];
    rcLJmaxsq[0][1] = rcLJmax[0][1]*rcLJmax[0][1];
    rcLJmaxsq[1][1] = rcLJmax[1][1]*rcLJmax[1][1];

    bLJmin[0][0] = bLJmin_CC;
    bLJmin[0][1] = bLJmin_CH;
    bLJmin[1][0] = bLJmin[0][1];
    bLJmin[1][1] = bLJmin_HH;

    bLJmax[0][0] = bLJmax_CC;
    bLJmax[0][1] = bLJmax_CH;
    bLJmax[1][0] = bLJmax[0][1];
    bLJmax[1][1] = bLJmax_HH;

    epsilon[0][0] = epsilon_CC;
    epsilon[0][1] = epsilon_CH;
    epsilon[1][0] = epsilon[0][1];
    epsilon[1][1] = epsilon_HH;

    sigma[0][0] = sigma_CC;
    sigma[0][1] = sigma_CH;
    sigma[1][0] = sigma[0][1];
    sigma[1][1] = sigma_HH;

    // torsional

    thmin = -1.0;
    thmax = -0.995;
    epsilonT[0][0] = epsilonT_CCCC;
    epsilonT[0][1] = epsilonT_CCCH;
    epsilonT[1][0] = epsilonT[0][1];
    epsilonT[1][1] = epsilonT_HCCH;
  


}
void scale_coords(int dir,int Natoms,point &L,atom *particle, double factor){
  int i;
  
  for(i=0;i<Natoms;i++) particle[i].r[dir] *= factor;
  
  if (dir==0)
    L.x *= factor;
  else if (dir==1)
    L.y *= factor;
  else
    L.z *= factor;
  //  *pbc *= factor;
  return;
}


void relax_cell(int Natoms,int *NBond,AtomParamsStruct *AtomParams,BondParamsStruct *BondParams,int *Bonded,double *Bmat,double *qvec,double *Avec,double *rij,double *rijShiled_1,point &L, atom *particle,Normal_struct *Normal,double *Fc,FILE* logfile,double *angle,int *BondNum,int *Neighb,int *NNeighb,int Coulomb_flag,int Interlayer,int *Normal_atom,int *NList,int *List,point *R0,int *NBond_once,int *Bonded_once){
  double origEnergy,origEnergy_;
  double newEnergy,newEnergy_;
  double scaleup = 1.001;
  double scaledown = 1./scaleup;
  ///static const char message[]="\rrelax box size %c: %.3f    ";
  
  int nchange=0, maxchange=1000,dir;
  //try to increase the volume
  //printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);

  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once);
  if(Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  origEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,origEnergy_,NBond_once,Bonded_once,1,0,logfile);
  
  for(dir=0;dir<3;dir++){
    scale_coords(dir,Natoms,L,particle,scaleup);
    //newEnergy=E_interlayer()+E_intralayer();
    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once);
    if(Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    newEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,newEnergy_,NBond_once,Bonded_once,1,0,logfile);
    
    if(newEnergy < origEnergy){ //keep increasing the volume
      origEnergy = newEnergy;
      while(nchange<maxchange){
	//printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);
	nchange++;
	scale_coords(dir,Natoms,L,particle,scaleup);
	calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once);
	if(Interlayer == 1)
	  calc_Normal(Natoms,Normal_atom,particle,Normal,L);
	
	newEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,newEnergy_,NBond_once,Bonded_once,1,0,logfile);
	
	if(newEnergy < origEnergy){
	  origEnergy = newEnergy;
	}
	else{ //step back and exit
	  //printf(message,dir,pow(scaleup,--nchange)); fflush(stdout);
	  nchange--;
	  scale_coords(dir,Natoms,L,particle,scaledown);
	  break;
	}
      }
    }
    else { //try to reduce the volume
      origEnergy = newEnergy;
      while(nchange>-maxchange){
	//printf(message,dir,pow(scaleup,--nchange)); fflush(stdout);
	nchange--;
	scale_coords(dir,Natoms,L,particle,scaledown);
	calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0,NBond_once,Bonded_once);
	if(Interlayer == 1)
	  calc_Normal(Natoms,Normal_atom,particle,Normal,L);
	
	newEnergy = Calc_Potential(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,newEnergy_,NBond_once,Bonded_once,1,0,logfile);
	
	if(newEnergy < origEnergy){
	  origEnergy = newEnergy;
	}
	else{ //step back and exit
	  //printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);
	  nchange++;
	  scale_coords(dir,Natoms,L,particle,scaleup);
	  break;
	}
      }
    }
  }
  //printf("\n"); fflush(stdout);
  return;
}

void TORSION_(atom *particle,point *f,int *BondNum,int Natoms,point L,int *Bonded,int *NBond,double &E_TORSION_,int *NBond_once,int *Bonded_once)
{
  double E_TORSION=0.0;
#pragma omp parallel
  {
    int i,j,k,l,ii,inum;
    //tagint itag,jtag;
    int itag,jtag;
    double evdwl,fpair,xtmp,ytmp,ztmp;
    double cos321;
    double w21,dw21,cos234,w34,dw34;
    double cross321[3],cross321mag,cross234[3],cross234mag;
    double w23,dw23,cw2,ekijl,Ec;
    double cw,cwnum,cwnom;
    double rij,rij2,rik,rjl,tspjik,dtsjik,tspijl,dtsijl,costmp,fcpc;
    double sin321,sin234,rjk2,rik2,ril2,rjl2;
    double rjk,ril;
    double Vtors;
    double dndij[3],tmpvec[3],dndik[3],dndjl[3];
    double dcidij,dcidik,dcidjk,dcjdji,dcjdjl,dcjdil;
    double dsidij,dsidik,dsidjk,dsjdji,dsjdjl,dsjdil;
    double dxidij,dxidik,dxidjk,dxjdji,dxjdjl,dxjdil;
    double ddndij,ddndik,ddndjk,ddndjl,ddndil,dcwddn,dcwdn,dvpdcw,Ftmp[3];
    double del32[3],rsq,r32,del23[3],del21[3],r21;
    double deljk[3],del34[3],delil[3],delkl[3],r23,r34;
    double fi[3],fj[3],fk[3],fl[3];
    int itype,jtype,ktype,ltype,kk,ll,jj;
    int *ilist,*REBO_neighs_i,*REBO_neighs_j;
    
    //double **x = atom->x;
    //double **f = atom->f;
    //int *type = atom->type;
    //tagint *tag = atom->tag;
    
    //inum = list->inum;
    //ilist = list->ilist;
    //map[type[k]]
    //REBO_numneigh[i]
    //for (j = 0; j < NBond[atomi]; j++) {
    // atomj = Bonded[atomi*MaxNBond + j];
    // itag = atomj;
    //}
    
#pragma omp for //reduction(+:E_TORSION)
    for (ii = 0; ii < Natoms; ii++) {
      //cerr<<"hey1"<<endl;
      //i = ilist[ii];
      i=ii;
      //itag = tag[i];
      itag = ii;
      itype = particle[i].type;
      if (itype != 0) continue;
      xtmp = particle[i].r[0];
      ytmp = particle[i].r[1];
      ztmp = particle[i].r[2];
      //REBO_neighs_i = REBO_firstneigh[i];
      //cerr<<"hey2"<<endl;
      for (jj = 0; jj < NBond_once[ii]; jj++) {
	//j = REBO_neighs_i[jj];
	j = Bonded_once[i*MaxNBond + jj];
	//jtag = tag[j];
	jtag = j;
	/*
	  if (itag > jtag) {
	  if ((itag+jtag) % 2 == 0) continue;
	  } else if (itag < jtag) {
	  if ((itag+jtag) % 2 == 1) continue;
	  } else {
	  if (particle[j].r[2] < ztmp) continue;
	  if (particle[j].r[2] == ztmp && particle[j].r[1] < ytmp) continue;
	  if (particle[j].r[2] == ztmp && particle[j].r[1] == ytmp && particle[j].r[0] < xtmp) continue;
	  }
	*/
      jtype = particle[j].type;
      if (jtype != 0) continue;
      
      del32[0]  = R_PBC(particle,j,i,0,L.x);
      del32[1]  = R_PBC(particle,j,i,1,L.y);
      del32[2]  = R_PBC(particle,j,i,2,L.z);
      rsq = del32[0]*del32[0] + del32[1]*del32[1] + del32[2]*del32[2];
      r32 = sqrt(rsq);
      del23[0] = -del32[0];
      del23[1] = -del32[1];
      del23[2] = -del32[2];
      r23 = r32;
      w23 = Sp(r23,rcmin[itype][jtype],rcmax[itype][jtype],dw23);
      //cerr<<"hey3"<<endl;
      for (kk = 0; kk < NBond[i]; kk++) {
        //k = REBO_neighs_i[kk];
	k = Bonded[i*MaxNBond + kk];
        ktype = particle[k].type;
        if (k == j) continue;
	del21[0] = R_PBC(particle,i,k,0,L.x);
	del21[1] = R_PBC(particle,i,k,1,L.y);
	del21[2] = R_PBC(particle,i,k,2,L.z);
        rsq = del21[0]*del21[0] + del21[1]*del21[1] + del21[2]*del21[2];
        r21 = sqrt(rsq);
        cos321 = - ((del21[0]*del32[0]) + (del21[1]*del32[1]) +
                    (del21[2]*del32[2])) / (r21*r32);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);
        sin321 = sqrt(1.0 - cos321*cos321);
        if (sin321 < TOL) continue;

        //deljk[0] = del21[0]-del23[0];
        //deljk[1] = del21[1]-del23[1];
        //deljk[2] = del21[2]-del23[2];
	
	deljk[0] = R_PBC(particle,j,k,0,L.x);
        deljk[1] = R_PBC(particle,j,k,1,L.y);
        deljk[2] = R_PBC(particle,j,k,2,L.z);
        rjk2 = deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2];
        rjk=sqrt(rjk2);
        rik2 = r21*r21;
        w21 = Sp(r21,rcmin[itype][ktype],rcmax[itype][ktype],dw21);

        rij = r32;
        rik = r21;
        rij2 = r32*r32;
        rik2 = r21*r21;
        costmp = 0.5*(rij2+rik2-rjk2)/rij/rik;
        tspjik = Sp2(costmp,thmin,thmax,dtsjik);
        dtsjik = -dtsjik;
	//cerr<<"hey4"<<endl;
        //REBO_neighs_j = REBO_firstneigh[j];
        for (ll = 0; ll <NBond[j]; ll++) {
          //l = REBO_neighs_j[ll];
	  l = Bonded[j*MaxNBond + ll];
          ltype = particle[l].type;
          if (l == i || l == k) continue;
	  del34[0] = R_PBC(particle,j,l,0,L.x);
	  del34[1] = R_PBC(particle,j,l,1,L.y);
	  del34[2] = R_PBC(particle,j,l,2,L.z);
          rsq = del34[0]*del34[0] + del34[1]*del34[1] + del34[2]*del34[2];
          r34 = sqrt(rsq);
          cos234 = (del32[0]*del34[0] + del32[1]*del34[1] +
                    del32[2]*del34[2]) / (r32*r34);
          cos234 = MIN(cos234,1.0);
          cos234 = MAX(cos234,-1.0);
          sin234 = sqrt(1.0 - cos234*cos234);
          if (sin234 < TOL) continue;
          w34 = Sp(r34,rcmin[jtype][ltype],rcmax[jtype][ltype],dw34);
          //delil[0] = del23[0] + del34[0];
          //delil[1] = del23[1] + del34[1];
          //delil[2] = del23[2] + del34[2];
	  delil[0] = R_PBC(particle,i,l,0,L.x);
	  delil[1] = R_PBC(particle,i,l,1,L.y);
	  delil[2] = R_PBC(particle,i,l,2,L.z);
          ril2 = delil[0]*delil[0] + delil[1]*delil[1] + delil[2]*delil[2];
          ril=sqrt(ril2);
          rjl2 = r34*r34;

          rjl = r34;
          rjl2 = r34*r34;
          costmp = 0.5*(rij2+rjl2-ril2)/rij/rjl;
          tspijl = Sp2(costmp,thmin,thmax,dtsijl);
          dtsijl = -dtsijl; //need minus sign
          cross321[0] = (del32[1]*del21[2])-(del32[2]*del21[1]);
          cross321[1] = (del32[2]*del21[0])-(del32[0]*del21[2]);
          cross321[2] = (del32[0]*del21[1])-(del32[1]*del21[0]);
          cross321mag = sqrt(cross321[0]*cross321[0]+
                             cross321[1]*cross321[1]+
                             cross321[2]*cross321[2]);
          cross234[0] = (del23[1]*del34[2])-(del23[2]*del34[1]);
          cross234[1] = (del23[2]*del34[0])-(del23[0]*del34[2]);
          cross234[2] = (del23[0]*del34[1])-(del23[1]*del34[0]);
          cross234mag = sqrt(cross234[0]*cross234[0]+
                             cross234[1]*cross234[1]+
                             cross234[2]*cross234[2]);
          cwnum = (cross321[0]*cross234[0]) +
            (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
          cwnom = r21*r34*r32*r32*sin321*sin234;
          cw = cwnum/cwnom;

          cw2 = (.5*(1.0-cw));
          ekijl = epsilonT[ktype][ltype];
          Ec = 256.0*ekijl/405.0;
          Vtors = (Ec*(powint(cw2,5)))-(ekijl/10.0);

          //if (eflag) evdwl = Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);

	  //evdwl =
#pragma omp critical
	  E_TORSION += Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);
	  //printf("E_TORSION= %.16f\n", E_TORSION);
          dndij[0] = (cross234[1]*del21[2])-(cross234[2]*del21[1]);
          dndij[1] = (cross234[2]*del21[0])-(cross234[0]*del21[2]);
          dndij[2] = (cross234[0]*del21[1])-(cross234[1]*del21[0]);

          tmpvec[0] = (del34[1]*cross321[2])-(del34[2]*cross321[1]);
          tmpvec[1] = (del34[2]*cross321[0])-(del34[0]*cross321[2]);
          tmpvec[2] = (del34[0]*cross321[1])-(del34[1]*cross321[0]);

          dndij[0] = dndij[0]+tmpvec[0];
          dndij[1] = dndij[1]+tmpvec[1];
          dndij[2] = dndij[2]+tmpvec[2];

          dndik[0] = (del23[1]*cross234[2])-(del23[2]*cross234[1]);
          dndik[1] = (del23[2]*cross234[0])-(del23[0]*cross234[2]);
          dndik[2] = (del23[0]*cross234[1])-(del23[1]*cross234[0]);

          dndjl[0] = (cross321[1]*del23[2])-(cross321[2]*del23[1]);
          dndjl[1] = (cross321[2]*del23[0])-(cross321[0]*del23[2]);
          dndjl[2] = (cross321[0]*del23[1])-(cross321[1]*del23[0]);

          dcidij = ((r23*r23)-(r21*r21)+(rjk*rjk))/(2.0*r23*r23*r21);
          dcidik = ((r21*r21)-(r23*r23)+(rjk*rjk))/(2.0*r23*r21*r21);
          dcidjk = (-rjk)/(r23*r21);
          dcjdji = ((r23*r23)-(r34*r34)+(ril*ril))/(2.0*r23*r23*r34);
          dcjdjl = ((r34*r34)-(r23*r23)+(ril*ril))/(2.0*r23*r34*r34);
          dcjdil = (-ril)/(r23*r34);

          dsidij = (-cos321/sin321)*dcidij;
          dsidik = (-cos321/sin321)*dcidik;
          dsidjk = (-cos321/sin321)*dcidjk;

          dsjdji = (-cos234/sin234)*dcjdji;
          dsjdjl = (-cos234/sin234)*dcjdjl;
          dsjdil = (-cos234/sin234)*dcjdil;

          dxidij = (r21*sin321)+(r23*r21*dsidij);
          dxidik = (r23*sin321)+(r23*r21*dsidik);
          dxidjk = (r23*r21*dsidjk);

          dxjdji = (r34*sin234)+(r23*r34*dsjdji);
          dxjdjl = (r23*sin234)+(r23*r34*dsjdjl);
          dxjdil = (r23*r34*dsjdil);

          ddndij = (dxidij*cross234mag)+(cross321mag*dxjdji);
          ddndik = dxidik*cross234mag;
          ddndjk = dxidjk*cross234mag;
          ddndjl = cross321mag*dxjdjl;
          ddndil = cross321mag*dxjdil;
          dcwddn = -cwnum/(cwnom*cwnom);
          dcwdn = 1.0/cwnom;
          dvpdcw = (-1.0)*Ec*(-.5)*5.0*powint(cw2,4) *
            w23*w21*w34*(1.0-tspjik)*(1.0-tspijl);

          Ftmp[0] = dvpdcw*((dcwdn*dndij[0])+(dcwddn*ddndij*del23[0]/r23));
          Ftmp[1] = dvpdcw*((dcwdn*dndij[1])+(dcwddn*ddndij*del23[1]/r23));
          Ftmp[2] = dvpdcw*((dcwdn*dndij[2])+(dcwddn*ddndij*del23[2]/r23));
          fi[0] = Ftmp[0];
          fi[1] = Ftmp[1];
          fi[2] = Ftmp[2];
          fj[0] = -Ftmp[0];
          fj[1] = -Ftmp[1];
          fj[2] = -Ftmp[2];

          Ftmp[0] = dvpdcw*((dcwdn*dndik[0])+(dcwddn*ddndik*del21[0]/r21));
          Ftmp[1] = dvpdcw*((dcwdn*dndik[1])+(dcwddn*ddndik*del21[1]/r21));
          Ftmp[2] = dvpdcw*((dcwdn*dndik[2])+(dcwddn*ddndik*del21[2]/r21));
          fi[0] += Ftmp[0];
          fi[1] += Ftmp[1];
          fi[2] += Ftmp[2];
          fk[0] = -Ftmp[0];
          fk[1] = -Ftmp[1];
          fk[2] = -Ftmp[2];

          Ftmp[0] = (dvpdcw*dcwddn*ddndjk*deljk[0])/rjk;
          Ftmp[1] = (dvpdcw*dcwddn*ddndjk*deljk[1])/rjk;
          Ftmp[2] = (dvpdcw*dcwddn*ddndjk*deljk[2])/rjk;
          fj[0] += Ftmp[0];
          fj[1] += Ftmp[1];
          fj[2] += Ftmp[2];
          fk[0] -= Ftmp[0];
          fk[1] -= Ftmp[1];
          fk[2] -= Ftmp[2];

          Ftmp[0] = dvpdcw*((dcwdn*dndjl[0])+(dcwddn*ddndjl*del34[0]/r34));
          Ftmp[1] = dvpdcw*((dcwdn*dndjl[1])+(dcwddn*ddndjl*del34[1]/r34));
          Ftmp[2] = dvpdcw*((dcwdn*dndjl[2])+(dcwddn*ddndjl*del34[2]/r34));
          fj[0] += Ftmp[0];
          fj[1] += Ftmp[1];
          fj[2] += Ftmp[2];
          fl[0] = -Ftmp[0];
          fl[1] = -Ftmp[1];
          fl[2] = -Ftmp[2];

          Ftmp[0] = (dvpdcw*dcwddn*ddndil*delil[0])/ril;
          Ftmp[1] = (dvpdcw*dcwddn*ddndil*delil[1])/ril;
          Ftmp[2] = (dvpdcw*dcwddn*ddndil*delil[2])/ril;
          fi[0] += Ftmp[0];
          fi[1] += Ftmp[1];
          fi[2] += Ftmp[2];
          fl[0] -= Ftmp[0];
          fl[1] -= Ftmp[1];
          fl[2] -= Ftmp[2];

          // coordination forces

          fpair = Vtors*dw21*w23*w34*(1.0-tspjik)*(1.0-tspijl) / r21;
          fi[0] -= del21[0]*fpair;
          fi[1] -= del21[1]*fpair;
          fi[2] -= del21[2]*fpair;
          fk[0] += del21[0]*fpair;
          fk[1] += del21[1]*fpair;
          fk[2] += del21[2]*fpair;

          fpair = Vtors*w21*dw23*w34*(1.0-tspjik)*(1.0-tspijl) / r23;
          fi[0] -= del23[0]*fpair;
          fi[1] -= del23[1]*fpair;
          fi[2] -= del23[2]*fpair;
          fj[0] += del23[0]*fpair;
          fj[1] += del23[1]*fpair;
          fj[2] += del23[2]*fpair;

          fpair = Vtors*w21*w23*dw34*(1.0-tspjik)*(1.0-tspijl) / r34;
          fj[0] -= del34[0]*fpair;
          fj[1] -= del34[1]*fpair;
          fj[2] -= del34[2]*fpair;
          fl[0] += del34[0]*fpair;
          fl[1] += del34[1]*fpair;
          fl[2] += del34[2]*fpair;

          // additional cut off function forces

          fcpc = -Vtors*w21*w23*w34*dtsjik*(1.0-tspijl);
          fpair = fcpc*dcidij/rij;
          fi[0] += fpair*del23[0];
          fi[1] += fpair*del23[1];
          fi[2] += fpair*del23[2];
          fj[0] -= fpair*del23[0];
          fj[1] -= fpair*del23[1];
          fj[2] -= fpair*del23[2];

          fpair = fcpc*dcidik/rik;
          fi[0] += fpair*del21[0];
          fi[1] += fpair*del21[1];
          fi[2] += fpair*del21[2];
          fk[0] -= fpair*del21[0];
          fk[1] -= fpair*del21[1];
          fk[2] -= fpair*del21[2];

          fpair = fcpc*dcidjk/rjk;
          fj[0] += fpair*deljk[0];
          fj[1] += fpair*deljk[1];
          fj[2] += fpair*deljk[2];
          fk[0] -= fpair*deljk[0];
          fk[1] -= fpair*deljk[1];
          fk[2] -= fpair*deljk[2];

          fcpc = -Vtors*w21*w23*w34*(1.0-tspjik)*dtsijl;
          fpair = fcpc*dcjdji/rij;
          fi[0] += fpair*del23[0];
          fi[1] += fpair*del23[1];
          fi[2] += fpair*del23[2];
          fj[0] -= fpair*del23[0];
          fj[1] -= fpair*del23[1];
          fj[2] -= fpair*del23[2];

          fpair = fcpc*dcjdjl/rjl;
          fj[0] += fpair*del34[0];
          fj[1] += fpair*del34[1];
          fj[2] += fpair*del34[2];
          fl[0] -= fpair*del34[0];
          fl[1] -= fpair*del34[1];
          fl[2] -= fpair*del34[2];

          fpair = fcpc*dcjdil/ril;
          fi[0] += fpair*delil[0];
          fi[1] += fpair*delil[1];
          fi[2] += fpair*delil[2];
          fl[0] -= fpair*delil[0];
          fl[1] -= fpair*delil[1];
          fl[2] -= fpair*delil[2];

          // sum per-atom forces into atom force array
#pragma omp critical
	  {
          f[i].x += fi[0]; f[i].y += fi[1]; f[i].z += fi[2];
          f[j].x += fj[0]; f[j].y += fj[1]; f[j].z += fj[2];
          f[k].x += fk[0]; f[k].y += fk[1]; f[k].z += fk[2];
          f[l].x += fl[0]; f[l].y += fl[1]; f[l].z += fl[2];
	  }
          //if (evflag) {
	  // delkl[0] = delil[0] - del21[0];
	  // delkl[1] = delil[1] - del21[1];
	  // delkl[2] = delil[2] - del21[2];
	  // ev_tally4(i,j,k,l,evdwl,fi,fj,fk,delil,del34,delkl);
	  // }
	  
        }
      }
    }
  }
  }
  //#pragma omp critical
  E_TORSION_=E_TORSION;
}

/****************************************************************/
double calc_dP_ij_dr(int atomi,int atomj, int atomk,int dir,atom *particle,int Natoms,double *rij,double r_ij,Normal_struct *Normal,point L,int *Normal_atom,int Interlayer,double rx,double ry,double rz,double inside_acos,double acos_inside_acos,double sqrt_inside,double sin_angle,double cos_angle,double Pij_term,double dR_ij,dr_ij *dN)
{
  
  double angle,dPij,dinside_acos,dangle;
  //double dR_ij;
    
  //dR_ij = calc_dr_ij_drk(atomi,atomj,atomk,dir,r_ij,particle,Natoms,L);
  
  //rx = particle[atomi].r[0] - particle[atomj].r[0];
  //ry = particle[atomi].r[1] - particle[atomj].r[1];
  //rz = particle[atomi].r[2] - particle[atomj].r[2];
  
  //rx -= L.x * rint(rx / (L.x + TINY));
  //ry -= L.y * rint(ry / (L.y + TINY));
  //rz -= L.z * rint(rz / (L.z + TINY));
  
  //inside_acos = (rx*Normal[atomi].x + ry*Normal[atomi].y + rz*Normal[atomi].z) / (r_ij);
  if (fabs(inside_acos) < 1.0){
    if (atomk == atomi){
      if (dir == 0)dinside_acos = (Normal[atomi].x)/(r_ij) - Pij_term*dR_ij;
      if (dir == 1)dinside_acos = (Normal[atomi].y)/(r_ij) - Pij_term*dR_ij;
      if (dir == 2)dinside_acos = (Normal[atomi].z)/(r_ij) - Pij_term*dR_ij;
    }
    else if (atomk == atomj){
      if (dir == 0)dinside_acos = (-Normal[atomi].x)/(r_ij) - Pij_term*dR_ij;
      if (dir == 1)dinside_acos = (-Normal[atomi].y)/(r_ij) - Pij_term*dR_ij;
      if (dir == 2)dinside_acos = (-Normal[atomi].z)/(r_ij) - Pij_term*dR_ij;
    }
    else{
      if (Interlayer == 1)
	calc_dNormal_k(atomi,atomk,dir,dN,particle,Normal,Normal_atom,L);
      dinside_acos = (rx*dN[0].r[0] + ry*dN[0].r[1] + rz*dN[0].r[2])/r_ij ;
    }
    dangle = -1.0/sqrt_inside*dinside_acos;
    dPij = dR_ij*sin_angle + r_ij*cos_angle*dangle;  
  }
  else {
    dPij = 0.0;
  }
    
  return(dPij);
  //return(dN[0].r[0]+dN[0].r[1]+dN[0].r[2]);
}

//*********************************
void calc_F_vdW(point *FvdW, int Natoms, atom *particle, double *rijShiled_1, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,double &EvdW)
{
  
  EvdW=0.0;
#pragma omp parallel
  {
    int atomk, atomi, typek, typei,atomj,typej, index,j;
    double alpha, rvdW, term1, term2,_term2, dEvdW_drij;
    double dTap, r_ij, Tap,C6,term6,term7;
    double R,dr_ij_drk,C6_r_ij;
    double reff,r_ik,r_jk;
    int n,dir;
    
#pragma omp parallel for 
    for(atomk=0 ; atomk < Natoms ; atomk++)
      FvdW[atomk].x = FvdW[atomk].y = FvdW[atomk].z = 0.0;
    
#pragma omp for
    for(atomi=0 ; atomi < Natoms ; atomi++){
      typei = particle[atomi].type;
      for (j=0 ; j < NNeighb[atomi]; j++){
	atomk = atomj = Neighb[atomi*MaxNeighb + j] ;	  
	typek = typej = particle[atomj].type;
		
	index = typej*NAtomParams + typei;
	//r_ij  = R_ij(particle,atomi,atomj, L);
	r_ij = rijShiled_1[atomi*MaxNeighb + j];
	Tap   = calc_Tap(r_ij);
	dTap  = calc_dTap(r_ij);  
	
	//interlayer vdW
	reff    = BondParams[index].r_eff;	    
	R       = BondParams[index].rvdW_long;
	C6      = BondParams[index].C6;
	
	term2   = exp(-d_TS*(r_ij/(Sr_TS*reff) - 1.0));
	_term2  = 1.0/(1.0 + term2);
	
	C6_r_ij = C6/(pow(r_ij,6.0)+0.1);
	
	term6   = (-(pow(_term2,2.0)*term2*(-d_TS/(Sr_TS*reff))*C6_r_ij ) + _term2*(-1)*C6*pow(pow(r_ij,6)+0.1,-2)*6*pow(r_ij,5));
	term7   = (-(_term2)*C6_r_ij)*dTap;

	EvdW += Tap*(-(1.0/(1.0 + exp(-d_TS*(r_ij/(Sr_TS*reff)- 1.0))))*C6/(pow(r_ij,6.0)+0.1));
	
	//********dirivative in the x direction******
	
	dir = 0;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dEvdW_drij = Tap*(- (term6*dr_ij_drk)) + dr_ij_drk*term7;
	FvdW[atomk].x -= dEvdW_drij;
	FvdW[atomi].x += dEvdW_drij;
	
	//********dirivative in the y direction******
	
	dir = 1;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dEvdW_drij = Tap*(- (term6*dr_ij_drk)) + dr_ij_drk*term7;
	
	FvdW[atomk].y -= dEvdW_drij;
	FvdW[atomi].y += dEvdW_drij;
	//********dirivative in the z direction******
	
	dir = 2;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	dEvdW_drij = Tap*(- (term6*dr_ij_drk)) + dr_ij_drk*term7;
	FvdW[atomk].z -= dEvdW_drij;
	FvdW[atomi].z += dEvdW_drij;
      }
    }
  }

}

void calc_F_Rep(point *FvdW, int Natoms, atom *particle, double *rijShiled_1, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,dr_ij *dN,double &ERep_)
{
  ERep_=0;
  
  double ERep=0.0;
#pragma omp parallel
  {
    
    int atomk, atomi, typek, typei,atomj,typej, index,j;
    double alpha, rvdW, Epsilon, term1, dEvdW_drij;
    double dTap, r_ij, Tap,Epsilon_C_vdW_exp,term3,term4,term5,term7;
    double dP_ij_drk,dP_ji_drk,R,Pij,Pji,gamma,dr_ij_drk,exp_term1;
    double exp_Pij_gamma, exp_Pji_gamma,Pij_gamma,Pji_gamma,reff,r_ik,r_jk;
    double rx_ij,ry_ij,rz_ij,rx_ji,ry_ji,rz_ji,inside_acos_ij,inside_acos_ji,Pij_term1,Pji_term1,Pij_term,Pji_term,dr_ji_drk;
    double angle_ij, sqrt_inside_ij, sin_angle_ij,cos_angle_ij,angle_ji,sqrt_inside_ji,sin_angle_ji,cos_angle_ji,dinside_acos,dangle;
    double C0_kc,C2_kc,C4_kc,A_kc,fji,fij,term6,dfji,dfij,A_term;
    double ERep_Tmp,dP_ij_term,dP_ji_term,d_Norm_length,Pow_normal_i,Pow_normal_j;
    int n,dir;
 
    dP_ij_drk = 0.0;
    dP_ji_drk = 0.0;

    C0_kc      = 15.71/1000;
    C2_kc      = 12.29/1000;
    C4_kc      = 4.933/1000;
    alpha   = 3.629;
    Epsilon = 3.03/1000;
    R       = 3.34;
    gamma   = 0.578;
    A_kc    = 10.238/1000;
    
#pragma omp for //schedule(dynamic,100) //reduction(+:ERep)
    for(atomi=0 ; atomi < Natoms ; atomi++){
      typei = particle[atomi].type;
      for (j=0 ; j < NNeighb[atomi]; j++){
	atomk = atomj = Neighb[atomi*MaxNeighb + j] ;	  
	typek = typej = particle[atomj].type;	
	
	index = typej*NAtomParams + typei;
	//r_ij  = R_ij(particle,atomi,atomj, L);
	r_ij = rijShiled_1[atomi*MaxNeighb + j];
	Tap   = calc_Tap(r_ij);
	dTap  = calc_dTap(r_ij);  
	
	reff    = BondParams[index].r_eff;	    
	
	dP_ij_drk=0.0;
	dP_ji_drk=0.0;
	Pij_term1=0.0;
	Pji_term1=0.0;
	dN[0].r[0]=0;
	dN[0].r[1]=0;
	dN[0].r[2]=0;
	rx_ij = particle[atomi].r[0] - particle[atomj].r[0];
	ry_ij = particle[atomi].r[1] - particle[atomj].r[1];
	rz_ij = particle[atomi].r[2] - particle[atomj].r[2];
	rx_ij -= L.x * rint(rx_ij / (L.x + TINY));
	ry_ij -= L.y * rint(ry_ij / (L.y + TINY));
	rz_ij -= L.z * rint(rz_ij / (L.z + TINY));
	
	rx_ji = -rx_ij;
	ry_ji = -ry_ij;
	rz_ji = -rz_ij;
	
	inside_acos_ij = (rx_ij*Normal[atomi].x + ry_ij*Normal[atomi].y + rz_ij*Normal[atomi].z) / (r_ij);
	inside_acos_ji = (rx_ji*Normal[atomj].x + ry_ji*Normal[atomj].y + rz_ji*Normal[atomj].z) / (r_ij);
	Pij_term1 = inside_acos_ij/r_ij;
	Pji_term1 = inside_acos_ji/r_ij;
		
	angle_ij  = acos(inside_acos_ij);
	sqrt_inside_ij=sqrt(1.0-sqr(inside_acos_ij));
	sin_angle_ij=sin(angle_ij);
	cos_angle_ij=cos(angle_ij);
	angle_ji  = acos(inside_acos_ji);
	sqrt_inside_ji=sqrt(1.0-sqr(inside_acos_ji));
	sin_angle_ji=sin(angle_ji);
	cos_angle_ji=cos(angle_ji);
	
	dP_ij_term =  r_ij*cos_angle_ij*(-1.0/sqrt_inside_ij);
	dP_ji_term =  r_ij*cos_angle_ji*(-1.0/sqrt_inside_ji);
	Pow_normal_i = pow(Normal[atomi].Normal_length,-3.0);
	Pow_normal_j = pow(Normal[atomj].Normal_length,-3.0);
	  
	if (fabs(inside_acos_ij) < 1.0){
	  Pij=r_ij*sin(angle_ij);
	}
	else Pij=0;
	
	if (fabs(inside_acos_ji) < 1.0){
	  Pji=r_ij*sin(angle_ji);
	}
	else Pji=0;
	//****
	
	Pij_gamma = Pij/gamma; //new
	Pji_gamma = Pji/gamma; //new
	exp_Pij_gamma = exp(-pow((Pij_gamma),2.0)); //new
	exp_Pji_gamma = exp(-pow((Pji_gamma),2.0)); //new
	
	exp_term1 = exp(-alpha*(r_ij-R));
	
	fij=0.0;
	fji=0.0;
	
	fij = C0_kc+ C2_kc*pow((Pij_gamma),2.0)+ C4_kc*pow((Pij_gamma),4.0);
	fji = C0_kc+ C2_kc*pow((Pji_gamma),2.0)+ C4_kc*pow((Pji_gamma),4.0);
	
	fij*=exp_Pij_gamma;
	fji*=exp_Pji_gamma;
	
	A_term=A_kc*pow((r_ij/R),-6);
	
	Epsilon_C_vdW_exp = Epsilon + fij + fji;
	
	ERep_Tmp = (exp_term1*(Epsilon_C_vdW_exp) - A_term);
	
#pragma omp critical
	ERep += Tap*(ERep_Tmp);
	
	dfij    = 2.0*C2_kc*Pij_gamma/gamma+4.0*C4_kc*pow(Pij_gamma,3)/gamma;
        dfji    = 2.0*C2_kc*Pji_gamma/gamma+4.0*C4_kc*pow(Pji_gamma,3)/gamma;
	//***
	term3   = exp_term1*(-alpha)*(Epsilon_C_vdW_exp);
	term4   = exp_term1*((-2.0)*Pij_gamma/gamma*fij+exp_Pij_gamma*dfij);
	term5   = exp_term1*((-2.0)*Pji_gamma/gamma*fji+exp_Pji_gamma*dfji);
	term6   = 6*A_kc*pow(r_ij/R,-7)/R;

	//********dirivative in the x direction******
	
	dir = 0;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	if (fabs(inside_acos_ij) < 1.0){
	  dinside_acos = (-Normal[atomi].x)/(r_ij) - Pij_term1*dr_ij_drk;
	  dangle = -1.0/sqrt_inside_ij*dinside_acos;
	  dP_ij_drk = dr_ij_drk*sin_angle_ij + r_ij*cos_angle_ij*dangle;
	}
	
	if (fabs(inside_acos_ji) < 1.0){
	  dinside_acos = (Normal[atomj].x)/(r_ij) - Pji_term1*dr_ij_drk;
	  dangle = -1.0/sqrt_inside_ji*dinside_acos;
	  dP_ji_drk = dr_ij_drk*sin_angle_ji + r_ij*cos_angle_ji*dangle;
	}
	
	dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk + term6*dr_ij_drk) + dTap*(ERep_Tmp)*dr_ij_drk;
	
	FvdW[atomk].x -= dEvdW_drij;
	FvdW[atomi].x += dEvdW_drij;
	
	//********dirivative in the y direction******
	
	dir = 1;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	if (fabs(inside_acos_ij) < 1.0){
	  dinside_acos = (-Normal[atomi].y)/(r_ij) - Pij_term1*dr_ij_drk;
	  dangle = -1.0/sqrt_inside_ij*dinside_acos;
	  dP_ij_drk = dr_ij_drk*sin_angle_ij + r_ij*cos_angle_ij*dangle;
	}
		
	if (fabs(inside_acos_ji) < 1.0){
	  dinside_acos = (Normal[atomj].y)/(r_ij) - Pji_term1*dr_ij_drk;
	  dangle = -1.0/sqrt_inside_ji*dinside_acos;
	  dP_ji_drk = dr_ij_drk*sin_angle_ji + r_ij*cos_angle_ji*dangle;
	}
		
	dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk + term6*dr_ij_drk) + dTap*(ERep_Tmp)*dr_ij_drk;
	FvdW[atomk].y -= dEvdW_drij;
	FvdW[atomi].y += dEvdW_drij;

	//********dirivative in the z direction******
	
	dir = 2;
	dr_ij_drk = calc_dr_ij_drk(atomj,atomi,atomk,dir,r_ij,particle,Natoms,L);
	
	if (fabs(inside_acos_ij) < 1.0){
	  dinside_acos = (-Normal[atomi].z)/(r_ij) - Pij_term1*dr_ij_drk;
	  dangle = -1.0/sqrt_inside_ij*dinside_acos;
	  dP_ij_drk = dr_ij_drk*sin_angle_ij + r_ij*cos_angle_ij*dangle;
	}
	
	if (fabs(inside_acos_ji) < 1.0){
	  dinside_acos = (Normal[atomj].z)/(r_ij) - Pji_term1*dr_ij_drk;
	  dangle = -1.0/sqrt_inside_ji*dinside_acos;
	  dP_ji_drk = dr_ij_drk*sin_angle_ji + r_ij*cos_angle_ji*dangle;
	}
	
	dEvdW_drij = Tap*(term3*dr_ij_drk + term4*dP_ij_drk + term5*dP_ji_drk + term6*dr_ij_drk) + dTap*(ERep_Tmp)*dr_ij_drk;
	
	FvdW[atomk].z -= dEvdW_drij;
	FvdW[atomi].z += dEvdW_drij;
	
	if (fabs(inside_acos_ij) < 1.0 && Interlayer == 1){
	  n=0;
	  
	  atomk =Normal_atom[atomi*3 + n];
	  
	  dir = 0;
	 	  
	  d_Norm_length = Pow_normal_i*(Normal[atomi].not_norm_y*(-Normal[atomi].vector_2z) + Normal[atomi].not_norm_z*Normal[atomi].vector_2y);
	  dN[0].r[0]= -Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomi].vector_2z)/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= Normal[atomi].vector_2y/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
	  
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].x -= dEvdW_drij;
	  
	  dir = 1;
	 	  
	  d_Norm_length = Pow_normal_i*(Normal[atomi].not_norm_x*Normal[atomi].vector_2z + Normal[atomi].not_norm_z*(-Normal[atomi].vector_2x));
	  dN[0].r[0]= Normal[atomi].vector_2z/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= -Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomi].vector_2x)/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].y -= dEvdW_drij;
	  
	  dir = 2;
	 	  
	  d_Norm_length = Pow_normal_i*(Normal[atomi].not_norm_x*(-Normal[atomi].vector_2y) + Normal[atomi].not_norm_y*(Normal[atomi].vector_2x) );
	  dN[0].r[0]= (-Normal[atomi].vector_2y)/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (Normal[atomi].vector_2x)/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= -Normal[atomi].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].z -= dEvdW_drij;
	  
	  n=1;
	  
	  atomk =Normal_atom[atomi*3 + n];
	  
	  dir = 0;
	 	  
	  d_Norm_length = Pow_normal_i*( Normal[atomi].not_norm_y*(-Normal[atomi].vector_1z + Normal[atomi].vector_2z) + Normal[atomi].not_norm_z*(-Normal[atomi].vector_2y + Normal[atomi].vector_1y));
	  dN[0].r[0]= -Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomi].vector_1z + Normal[atomi].vector_2z)/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomi].vector_2y + Normal[atomi].vector_1y)/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
	  
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].x -= dEvdW_drij;
	  
	  dir = 1;
	 	  
	  //dN[0].r[0] = -Normal[atomi].vector_2z + Normal[atomi].vector_1z;
	  //dN[0].r[1] = 0.0;
	  //dN[0].r[2] = -Normal[atomi].vector_1x + Normal[atomi].vector_2x;
	  d_Norm_length = Pow_normal_i*(Normal[atomi].not_norm_x*(-Normal[atomi].vector_2z + Normal[atomi].vector_1z) +  Normal[atomi].not_norm_z*(-Normal[atomi].vector_1x + Normal[atomi].vector_2x));
	  dN[0].r[0]= (-Normal[atomi].vector_2z + Normal[atomi].vector_1z)/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= -Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomi].vector_1x + Normal[atomi].vector_2x)/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].y -= dEvdW_drij;
	  
	  dir = 2;
	  dr_ij_drk = 0.0;
	  
	  //dN[0].r[0] =-Normal[atomi].vector_1y + Normal[atomi].vector_2y;
	  //dN[0].r[1] =-Normal[atomi].vector_2x + Normal[atomi].vector_1x;
	  //dN[0].r[2] = 0.0;
	  d_Norm_length = Pow_normal_i*(Normal[atomi].not_norm_x*(-Normal[atomi].vector_1y + Normal[atomi].vector_2y) + Normal[atomi].not_norm_y*(-Normal[atomi].vector_2x + Normal[atomi].vector_1x));
	  dN[0].r[0]= (-Normal[atomi].vector_1y + Normal[atomi].vector_2y)/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomi].vector_2x + Normal[atomi].vector_1x)/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= -Normal[atomi].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].z -= dEvdW_drij;
	  
	  n=2;
	  
	  atomk =Normal_atom[atomi*3 + n];
	  
	  dir = 0;
	  dr_ij_drk = 0.0;
	  //dN[0].r[0] = 0.0;
	  //dN[0].r[1] =Normal[atomi].vector_1z;
	  //dN[0].r[2] =-Normal[atomi].vector_1y;
	  d_Norm_length = Pow_normal_i*( Normal[atomi].not_norm_y*(Normal[atomi].vector_1z) + Normal[atomi].not_norm_z*(-Normal[atomi].vector_1y));
	  dN[0].r[0]= -Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (Normal[atomi].vector_1z)/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomi].vector_1y)/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
	  
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].x -= dEvdW_drij;
	  
	  dir = 1;
	  dr_ij_drk = 0.0;
	  
	  //dN[0].r[0] =-Normal[atomi].vector_1z;
	  //dN[0].r[1] = 0.0;
	  //dN[0].r[2] =Normal[atomi].vector_1x;
	  d_Norm_length = Pow_normal_i*(Normal[atomi].not_norm_x*(-Normal[atomi].vector_1z) + Normal[atomi].not_norm_z*(Normal[atomi].vector_1x));
	  dN[0].r[0]= (-Normal[atomi].vector_1z)/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= -Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (Normal[atomi].vector_1x)/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].y -= dEvdW_drij;
	  
	  dir = 2;
	  dr_ij_drk = 0.0;
	  
	  //calc_dNormal_k(atomi,atomk,dir,dN,particle,Normal,Normal_atom,L);
	  //dN[0].r[0] =Normal[atomi].vector_1y;
	  //dN[0].r[1] =-Normal[atomi].vector_1x;
	  //dN[0].r[2] = 0.0;
	  d_Norm_length = Pow_normal_i*(Normal[atomi].not_norm_x*(Normal[atomi].vector_1y) + Normal[atomi].not_norm_y*(-Normal[atomi].vector_1x));
	  dN[0].r[0]= (Normal[atomi].vector_1y)/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomi].vector_1x)/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
	  dN[0].r[2]= -Normal[atomi].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ij*dN[0].r[0] + ry_ij*dN[0].r[1] + rz_ij*dN[0].r[2])/r_ij ;
	  dP_ij_drk = dP_ij_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term4*dP_ij_drk );
	  FvdW[atomk].z -= dEvdW_drij;
	  
	}
	if (fabs(inside_acos_ji) < 1.0 && Interlayer == 1 ){
	  n=0;
	  
	  atomk =Normal_atom[atomj*3 + n];
	  
	  dir = 0;
	 	  
	  d_Norm_length = Pow_normal_j*(Normal[atomj].not_norm_y*(-Normal[atomj].vector_2z) + Normal[atomj].not_norm_z*Normal[atomj].vector_2y);
	  dN[0].r[0]= -Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomj].vector_2z)/Normal[atomj].Normal_length-Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= Normal[atomj].vector_2y/Normal[atomj].Normal_length-Normal[atomj].not_norm_z*d_Norm_length;
	  
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].x -= dEvdW_drij;
	  
	  dir = 1;
	 	  
	  d_Norm_length = Pow_normal_j*(Normal[atomj].not_norm_x*Normal[atomj].vector_2z + Normal[atomj].not_norm_z*(-Normal[atomj].vector_2x));
	  dN[0].r[0]= Normal[atomj].vector_2z/Normal[atomj].Normal_length-Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= -Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomj].vector_2x)/Normal[atomj].Normal_length-Normal[atomj].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].y -= dEvdW_drij;
	  
	  dir = 2;
	 	  
	  d_Norm_length = Pow_normal_j*(Normal[atomj].not_norm_x*(-Normal[atomj].vector_2y) + Normal[atomj].not_norm_y*(Normal[atomj].vector_2x) );
	  dN[0].r[0]= (-Normal[atomj].vector_2y)/Normal[atomj].Normal_length-Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (Normal[atomj].vector_2x)/Normal[atomj].Normal_length-Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= -Normal[atomj].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].z -= dEvdW_drij;
	  
	  n=1;
	  
	  atomk =Normal_atom[atomj*3 + n];
	  
	  dir = 0;
	 	  
	  d_Norm_length = Pow_normal_j*( Normal[atomj].not_norm_y*(-Normal[atomj].vector_1z + Normal[atomj].vector_2z) + Normal[atomj].not_norm_z*(-Normal[atomj].vector_2y + Normal[atomj].vector_1y));
	  dN[0].r[0]= -Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomj].vector_1z + Normal[atomj].vector_2z)/Normal[atomj].Normal_length-Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomj].vector_2y + Normal[atomj].vector_1y)/Normal[atomj].Normal_length-Normal[atomj].not_norm_z*d_Norm_length;
	  
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].x -= dEvdW_drij;
	  
	  dir = 1;
	 	  
	  //dN[0].r[0] = -Normal[atomj].vector_2z + Normal[atomj].vector_1z;
	  //dN[0].r[1] = 0.0;
	  //dN[0].r[2] = -Normal[atomj].vector_1x + Normal[atomj].vector_2x;
	  d_Norm_length = Pow_normal_j*(Normal[atomj].not_norm_x*(-Normal[atomj].vector_2z + Normal[atomj].vector_1z) +  Normal[atomj].not_norm_z*(-Normal[atomj].vector_1x + Normal[atomj].vector_2x));
	  dN[0].r[0]= (-Normal[atomj].vector_2z + Normal[atomj].vector_1z)/Normal[atomj].Normal_length-Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= -Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomj].vector_1x + Normal[atomj].vector_2x)/Normal[atomj].Normal_length-Normal[atomj].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].y -= dEvdW_drij;
	  
	  dir = 2;
	  dr_ij_drk = 0.0;
	  
	  //dN[0].r[0] =-Normal[atomj].vector_1y + Normal[atomj].vector_2y;
	  //dN[0].r[1] =-Normal[atomj].vector_2x + Normal[atomj].vector_1x;
	  //dN[0].r[2] = 0.0;
	  d_Norm_length = Pow_normal_j*(Normal[atomj].not_norm_x*(-Normal[atomj].vector_1y + Normal[atomj].vector_2y) + Normal[atomj].not_norm_y*(-Normal[atomj].vector_2x + Normal[atomj].vector_1x));
	  dN[0].r[0]= (-Normal[atomj].vector_1y + Normal[atomj].vector_2y)/Normal[atomj].Normal_length-Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomj].vector_2x + Normal[atomj].vector_1x)/Normal[atomj].Normal_length-Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= -Normal[atomj].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].z -= dEvdW_drij;
	  
	  n=2;
	  
	  atomk =Normal_atom[atomj*3 + n];
	  
	  dir = 0;
	  dr_ij_drk = 0.0;
	  //dN[0].r[0] = 0.0;
	  //dN[0].r[1] =Normal[atomj].vector_1z;
	  //dN[0].r[2] =-Normal[atomj].vector_1y;
	  d_Norm_length = Pow_normal_j*( Normal[atomj].not_norm_y*(Normal[atomj].vector_1z) + Normal[atomj].not_norm_z*(-Normal[atomj].vector_1y));
	  dN[0].r[0]= -Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (Normal[atomj].vector_1z)/Normal[atomj].Normal_length-Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (-Normal[atomj].vector_1y)/Normal[atomj].Normal_length-Normal[atomj].not_norm_z*d_Norm_length;
	  
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].x -= dEvdW_drij;
	  
	  dir = 1;
	  dr_ij_drk = 0.0;
	  
	  //dN[0].r[0] =-Normal[atomj].vector_1z;
	  //dN[0].r[1] = 0.0;
	  //dN[0].r[2] =Normal[atomj].vector_1x;
	  d_Norm_length = Pow_normal_j*(Normal[atomj].not_norm_x*(-Normal[atomj].vector_1z) + Normal[atomj].not_norm_z*(Normal[atomj].vector_1x));
	  dN[0].r[0]= (-Normal[atomj].vector_1z)/Normal[atomj].Normal_length-Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= -Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= (Normal[atomj].vector_1x)/Normal[atomj].Normal_length-Normal[atomj].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].y -= dEvdW_drij;
	  
	  dir = 2;
	  dr_ij_drk = 0.0;
	  
	  //calc_dNormal_k(atomj,atomk,dir,dN,particle,Normal,Normal_atom,L);
	  //dN[0].r[0] =Normal[atomj].vector_1y;
	  //dN[0].r[1] =-Normal[atomj].vector_1x;
	  //dN[0].r[2] = 0.0;
	  d_Norm_length = Pow_normal_j*(Normal[atomj].not_norm_x*(Normal[atomj].vector_1y) + Normal[atomj].not_norm_y*(-Normal[atomj].vector_1x));
	  dN[0].r[0]= (Normal[atomj].vector_1y)/Normal[atomj].Normal_length-Normal[atomj].not_norm_x*d_Norm_length;
	  dN[0].r[1]= (-Normal[atomj].vector_1x)/Normal[atomj].Normal_length-Normal[atomj].not_norm_y*d_Norm_length;
	  dN[0].r[2]= -Normal[atomj].not_norm_z*d_Norm_length;
	  dinside_acos = (rx_ji*dN[0].r[0] + ry_ji*dN[0].r[1] + rz_ji*dN[0].r[2])/r_ij ;
	  dP_ji_drk = dP_ji_term*dinside_acos;
	  
	  dEvdW_drij = Tap*( term5*dP_ji_drk );
	  FvdW[atomk].z -= dEvdW_drij;
	  
	}
      }
    }
  }
  ERep_=ERep;
}
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
void spline_init()
{
  int i,j,k;

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      PCCf[i][j] = 0.0;
      PCCdfdx[i][j] = 0.0;
      PCCdfdy[i][j] = 0.0;
      PCHf[i][j] = 0.0;
      PCHdfdx[i][j] = 0.0;
      PCHdfdy[i][j] = 0.0;
    }
  }

  PCCf[0][2] = -0.00050;
  PCCf[0][3] = 0.0161253646;
  PCCf[1][1] = -0.010960;
  PCCf[1][2] = 0.00632624824;
  PCCf[2][0] = -0.0276030;
  PCCf[2][1] = 0.00317953083;

  PCHf[0][1] = 0.209336733;
  PCHf[0][2] = -0.0644496154;
  PCHf[0][3] = -0.303927546;
  PCHf[1][0] = 0.010;
  PCHf[1][1] = -0.125123401;
  PCHf[1][2] = -0.298905246;
  PCHf[2][0] = -0.122042146;
  PCHf[2][1] = -0.300529172;
  PCHf[3][0] = -0.307584705;

  for (i = 0; i < 5; i++) {
    for (j = 0; j < 5; j++) {
      for (k = 0; k < 10; k++) {
        piCCf[i][j][k] = 0.0;
        piCCdfdx[i][j][k] = 0.0;
        piCCdfdy[i][j][k] = 0.0;
        piCCdfdz[i][j][k] = 0.0;
        piCHf[i][j][k] = 0.0;
        piCHdfdx[i][j][k] = 0.0;
        piCHdfdy[i][j][k] = 0.0;
        piCHdfdz[i][j][k] = 0.0;
        piHHf[i][j][k] = 0.0;
        piHHdfdx[i][j][k] = 0.0;
        piHHdfdy[i][j][k] = 0.0;
        piHHdfdz[i][j][k] = 0.0;
        Tf[i][j][k] = 0.0;
        Tdfdx[i][j][k] = 0.0;
        Tdfdy[i][j][k] = 0.0;
        Tdfdz[i][j][k] = 0.0;
      }
    }
  }

  for (i = 3; i < 10; i++) piCCf[0][0][i] = 0.0049586079;
  piCCf[1][0][1] = 0.021693495;
  piCCf[0][1][1] = 0.021693495;
  for (i = 2; i < 10; i++) piCCf[1][0][i] = 0.0049586079;
  for (i = 2; i < 10; i++) piCCf[0][1][i] = 0.0049586079;
  piCCf[1][1][1] = 0.05250;
  piCCf[1][1][2] = -0.002088750;
  for (i = 3; i < 10; i++) piCCf[1][1][i] = -0.00804280;
  piCCf[2][0][1] = 0.024698831850;
  piCCf[0][2][1] = 0.024698831850;
  piCCf[2][0][2] = -0.00597133450;
  piCCf[0][2][2] = -0.00597133450;
  for (i = 3; i < 10; i++) piCCf[2][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][2][i] = 0.0049586079;
  piCCf[2][1][1] = 0.00482478490;
  piCCf[1][2][1] = 0.00482478490;
  piCCf[2][1][2] = 0.0150;
  piCCf[1][2][2] = 0.0150;
  piCCf[2][1][3] = -0.010;
  piCCf[1][2][3] = -0.010;
  piCCf[2][1][4] = -0.01168893870;
  piCCf[1][2][4] = -0.01168893870;
  piCCf[2][1][5] = -0.013377877400;
  piCCf[1][2][5] = -0.013377877400;
  piCCf[2][1][6] = -0.015066816000;
  piCCf[1][2][6] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[2][1][i] = -0.015066816000;
  for (i = 7; i < 10; i++) piCCf[1][2][i] = -0.015066816000;
  piCCf[2][2][1] = 0.0472247850;
  piCCf[2][2][2] = 0.0110;
  piCCf[2][2][3] = 0.0198529350;
  piCCf[2][2][4] = 0.01654411250;
  piCCf[2][2][5] = 0.013235290;
  piCCf[2][2][6] = 0.00992646749999 ;
  piCCf[2][2][7] = 0.006617644999;
  piCCf[2][2][8] = 0.00330882250;
  piCCf[3][0][1] = -0.05989946750;
  piCCf[0][3][1] = -0.05989946750;
  piCCf[3][0][2] = -0.05989946750;
  piCCf[0][3][2] = -0.05989946750;
  for (i = 3; i < 10; i++) piCCf[3][0][i] = 0.0049586079;
  for (i = 3; i < 10; i++) piCCf[0][3][i] = 0.0049586079;
  piCCf[3][1][2] = -0.0624183760;
  piCCf[1][3][2] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[3][1][i] = -0.0624183760;
  for (i = 3; i < 10; i++) piCCf[1][3][i] = -0.0624183760;
  piCCf[3][2][1] = -0.02235469150;
  piCCf[2][3][1] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[3][2][i] = -0.02235469150;
  for (i = 2; i < 10; i++) piCCf[2][3][i] = -0.02235469150;

  piCCdfdx[2][1][1] = -0.026250;
  piCCdfdx[2][1][5] = -0.0271880;
  piCCdfdx[2][1][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdx[2][1][i] = -0.0271880;
  piCCdfdx[1][3][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdx[2][3][i] = 0.031209;

  piCCdfdy[1][2][1] = -0.026250;
  piCCdfdy[1][2][5] = -0.0271880;
  piCCdfdy[1][2][6] = -0.0271880;
  for (i = 7; i < 10; i++) piCCdfdy[1][2][i] = -0.0271880;
  piCCdfdy[3][1][2] = 0.0187723882;
  for (i = 2; i < 10; i++) piCCdfdy[3][2][i] = 0.031209;

  piCCdfdz[1][1][2] = -0.0302715;
  piCCdfdz[2][1][4] = -0.0100220;
  piCCdfdz[1][2][4] = -0.0100220;
  piCCdfdz[2][1][5] = -0.0100220;
  piCCdfdz[1][2][5] = -0.0100220;
  for (i = 4; i < 9; i++) piCCdfdz[2][2][i] = -0.0033090;

  //  make top end of piCC flat instead of zero
  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCCf[i][j][k] = piCCf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){ // also enforces some symmetry
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCCf[i][j][k] = piCCf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCCf[4][4][k] = piCCf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCCf[i][j][k] = piCCf[i][j][k-1];
      }
  }

  piCHf[1][1][1] = -0.050;
  piCHf[1][1][2] = -0.050;
  piCHf[1][1][3] = -0.30;
  for (i = 4; i < 10; i++) piCHf[1][1][i] = -0.050;
  for (i = 5; i < 10; i++) piCHf[2][0][i] = -0.004523893758064;
  for (i = 5; i < 10; i++) piCHf[0][2][i] = -0.004523893758064;
  piCHf[2][1][2] = -0.250;
  piCHf[1][2][2] = -0.250;
  piCHf[2][1][3] = -0.250;
  piCHf[1][2][3] = -0.250;
  piCHf[3][1][1] = -0.10;
  piCHf[1][3][1] = -0.10;
  piCHf[3][1][2] = -0.125;
  piCHf[1][3][2] = -0.125;
  piCHf[3][1][3] = -0.125;
  piCHf[1][3][3] = -0.125;
  for (i = 4; i < 10; i++) piCHf[3][1][i] = -0.10;
  for (i = 4; i < 10; i++) piCHf[1][3][i] = -0.10;

  // make top end of piCH flat instead of zero
 // also enforces some symmetry

  i = 4;
  for (j = 0; j < 4; j++){
      for (k = 1; k < 11; k++){
          piCHf[i][j][k] = piCHf[i-1][j][k];
      }
  }
  for (i = 0; i < 4; i++){
      for (j = i+1; j < 5; j++){
          for (k = 1; k < 11; k++){
              piCHf[i][j][k] = piCHf[j][i][k];
          }
      }
  }
  for (k = 1; k < 11; k++) piCHf[4][4][k] = piCHf[3][4][k];
  k = 10;
  for (i = 0; i < 5; i++){
      for (j = 0; j < 5; j++){
      piCHf[i][j][k] = piCHf[i][j][k-1];
      }
  }

  piHHf[1][1][1] = 0.124915958;

  Tf[2][2][1] = -0.035140;
  for (i = 2; i < 10; i++) Tf[2][2][i] = -0.0040480;
}
double Calc_Potential(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec,double *rijShiled_1, FILE* logfile,double *angle,Normal_struct *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double &E_tot_,int *NBond_once,int *Bonded_once,int print_data,int counter,FILE* analysis)
{
  double E_tot;
  int atomi;
  double E_TORSION,E_KC,E_REBO;
  double E_TORSION_,E_KC_,E_REBO_;
  double ERep,EDisp,RI,Interlayer_dis;
  dr_ij *dN = new dr_ij[1];
  
  E_TORSION = 0.0; 
  E_REBO    = 0.0;
  E_KC      = 0.0;
  E_tot     = 0.0;
  
  //E_REBO = calc_dE_REBO_pot(angle,BondNum,Natoms,L,BondParams,Bonded,particle,Fc,NBond,E_REBO_,NBond_once,Bonded_once);
  
  if (Coulomb_flag==1 || Coulomb_flag==2)
    //E_TORSION = TORSION_pot(particle,BondNum,Natoms,L,Bonded,NBond,E_TORSION_,NBond_once,Bonded_once);
  
  if (Coulomb_flag==2)
    E_KC = calc_F_Rep_pot(Natoms,particle,rijShiled_1,BondParams,L,Normal,Normal_atom,Interlayer,Neighb,NNeighb,dN,E_KC_,ERep,EDisp,RI,Interlayer_dis);
  
  delete [] dN;
  //cerr<<E_REBO + E_TORSION + E_KC<<endl;
  E_tot = E_REBO + E_TORSION + E_KC;
  
  //if(print_data==0)
  //printf("Etot=%.8f REBO= %.8f TORSION= %.8f Rep= %.8f Disp= %.8f RI= %.4f Interlayer_dis=%.4f\n",E_tot,E_REBO,E_TORSION,ERep,EDisp,RI,Interlayer_dis);
  
  if(print_data==0)
    fprintf(analysis,"%i %.8f %.8f %.8f %.8f %.8f %.4f %.4f\n",counter,E_tot,E_REBO,E_TORSION,ERep,EDisp,RI,Interlayer_dis);
  
  fflush(analysis);
  
  return(E_tot);
  
  
}
double calc_dE_REBO_pot(double *angle,int *BondNum,int Natoms,point L,BondParamsStruct *BondParams,int *Bonded,atom *particle,double *Fc,int *NBond,double &E_REBO_,int *NBond_once,int *Bonded_once)
{
  double E_REBO;
  E_REBO=0.0;
#pragma omp parallel
  {
    int index,atomk,typek,typej,typei,index2;
    double g,Theta_ijk,b_ij,r_ij,r_jk,Tersoff_lambda1,Tersoff_lambda2,Tersoff_lambda3,Tersoff_c,Tersoff_d,Tersoff_h,Tersoff_n,A,B,beta;
    double Xsi,ETersoff,dTheta_ijk_drn,dg_drn,dXsi_drn,Fc_ij,Fc_jk,Tersoff_R,Tersoff_D;
    double dr_ij_drn,dr_jk_drn,dFc_ij_drn,dFc_jk_drn,Fa,Fr,dFa_drn,dFr_drn,db_ij_drn,n,dE_Tersoff_ij_drn;
    double Tersoff_sqr_c,Tersoff_sqr_d,sqr_Tersoff_h_Cos_Theta;
    //i=BondNum[atomj*Natoms + atomi];
    double inside_acos,cos_ijk;
    //cerr<<"min[0][0]2 ="<<rcmin[0][0]<<endl;
    //exit(0);
    
    int i,j,k,m,ii,itype,jtype,atomi,atomj;
    // tagint itag,jtag;
    int itag ,jtag;
    double delx,dely,delz,evdwl,fpair,xtmp,ytmp,ztmp;
    double rsq,rij,wij;
    double Qij,Aij,alphaij,VR,pre,dVRdi,VA,term,bij,dVAdi,dVA;
    double dwij,del[3];
    //double E_REBO;
    
    int vflag_atom=0;
    evdwl = 0.0;
    
    // two-body interactions from REBO neighbor list, skip half of them
#pragma omp for reduction(+:E_REBO)
    for (atomi = 0; atomi < Natoms; atomi++) {
      itag=atomi;
      itype = particle[atomi].type;
      xtmp = particle[atomi].r[0];
      ytmp = particle[atomi].r[1];
      ztmp = particle[atomi].r[2];
      //REBO_neighs = REBO_firstneigh[i];
      
      for (j = 0; j < NBond_once[atomi]; j++) {
	atomj = Bonded_once[atomi*MaxNBond + j];

	jtype = particle[atomj].type;

	delx = R_PBC(particle,atomi,atomj,0,L.x);
	dely = R_PBC(particle,atomi,atomj,1,L.y);	
	delz = R_PBC(particle,atomi,atomj,2,L.z);
	
	rsq = delx*delx + dely*dely + delz*delz;
	rij = sqrt(rsq);
	
	wij = Sp(rij,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
	
	if (wij <= TOL) continue;
	Qij = Q[itype][jtype];
	
	Aij = A_ij[itype][jtype];
	
	alphaij = alpha[itype][jtype];
	
	VR = wij*(1.0+(Qij/rij)) * Aij*exp(-alphaij*rij);

	pre = wij*Aij * exp(-alphaij*rij);
		
	VA = dVA = 0.0;
	
	for (m = 0; m < 3; m++) {
	  term = -wij * BIJc[itype][jtype][m] * exp(-Beta[itype][jtype][m]*rij);
	  VA += term;
	}
	
	del[0] = delx;
	del[1] = dely;
	del[2] = delz;
	
	bij = bondorder_pot(atomi,atomj,del,rij,VA,vflag_atom,particle,Bonded,BondNum,NBond,L);
	//bij=1;
	
	//#pragma omp atomic
	E_REBO+=VR + bij*VA;
     	
      }
    }
  }
  //cerr<<"E_REBO_03"<<E_REBO<<endl;
  return(E_REBO);
}


double bondorder_pot(int i, int j, double rij[3], double rijmag, double VA, int vflag_atom,atom *particle,int *Bonded,int *BondNum,int *NBond,point L){

  int atomi,atomj,k,n,l,atomk,atoml,atomn,atom1,atom2,atom3,atom4;
  int itype,jtype,ktype,ltype,ntype;
  double rik[3],rjl[3],rkn[3],rji[3],rki[3],rlj[3],rknmag,dNki,dwjl,bij;
  double NijC,NijH,NjiC,NjiH,wik,dwik,dwkn,wjl;
  double rikmag,rjlmag,cosjik,cosijl,g,tmp2,tmp3;
  double Etmp,pij,tmp,wij,dwij,NconjtmpI,NconjtmpJ,Nki,Nlj,dS;
  double lamdajik,lamdaijl,dgdc,dgdN,pji,Nijconj,piRC;
  double dcosjikdri[3],dcosijldri[3],dcosjikdrk[3];
  double dN2[2],dN3[3];
  double dcosjikdrj[3],dcosijldrj[3],dcosijldrl[3];
  double Tij;
  double r32[3],r32mag,cos321,r43[3],r13[3];
  double dNlj;
  double om1234,rln[3];
  double rlnmag,dwln,r23[3],r23mag,r21[3],r21mag;
  double w21,dw21,r34[3],r34mag,cos234,w34,dw34;
  double cross321[3],cross234[3],prefactor,SpN;
  double fcijpc,fcikpc,fcjlpc,fcjkpc,fcilpc;
  double dt2dik[3],dt2djl[3],dt2dij[3],aa,aaa1,aaa2,at2,cw,cwnum,cwnom;
  double sin321,sin234,rr,rijrik,rijrjl,rjk2,rik2,ril2,rjl2;
  double dctik,dctjk,dctjl,dctij,dctji,dctil,rik2i,rjl2i,sink2i,sinl2i;
  double rjk[3],ril[3],dt1dik,dt1djk,dt1djl,dt1dil,dt1dij;
  double F23[3],F12[3],F34[3],F31[3],F24[3],fi[3],fj[3],fk[3],fl[3];
  double f1[3],f2[3],f3[3],f4[4];
  double dcut321,PijS,PjiS;
  double rij2,tspjik,dtsjik,tspijl,dtsijl,costmp;
  int *REBO_neighs,*REBO_neighs_i,*REBO_neighs_j,*REBO_neighs_k,*REBO_neighs_l;

  //double **x = atom->x;
  //int *type = atom->type;

  atomi = i;
  atomj = j;
  itype = particle[atomi].type;
  jtype = particle[atomj].type;

  //index      = itype*NAtomParams + jtype;
  wij = Sp(rijmag,rcmin[itype][jtype],rcmax[itype][jtype],dwij);
  //printf("jtype=%i itype=%i wij=%.16f\n",jtype,itype,wij);  
  //  wij = Sp(rijmag,rcmin[index],rcmax[index],dwij);
  NijC = nC[i]-(wij*kronecker(jtype,0));
  NijH = nH[i]-(wij*kronecker(jtype,1));
  NjiC = nC[j]-(wij*kronecker(itype,0));
  NjiH = nH[j]-(wij*kronecker(itype,1));

  //printf("jtype=%i itype=%i wij=%.16f nC[i]=%.16f nH[i]=%.16f nC[j]=%.16f nH[j]=%.16f \n",jtype,itype,wij,nC[i],nH[i],nC[j],nH[j]);
  
  bij = 0.0;
  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  dgdc = 0.0;
  dgdN = 0.0;
  NconjtmpI = 0.0;
  NconjtmpJ = 0.0;
  Etmp = 0.0;
  
  //REBO_neighs = REBO_firstneigh[i];
  
    for(k=0 ; k < NBond[atomi] ; k++){ // Go over all j neighbors k != i.
      //atomk = REBO_neighs[k];
      atomk = Bonded[atomi*MaxNBond + k];
    if (atomk != atomj) {
      ktype = particle[atomk].type;
      //ktype = particle[atomk].type;
      rik[0] = R_PBC(particle,atomi,atomk,0,L.x);
      rik[1] = R_PBC(particle,atomi,atomk,1,L.y);
      rik[2] = R_PBC(particle,atomi,atomk,2,L.z);
      
      rikmag = sqrt((rik[0]*rik[0])+(rik[1]*rik[1])+(rik[2]*rik[2]));
      lamdajik = 4.0*kronecker(itype,1) *
        ((rho[ktype][1]-rikmag)-(rho[jtype][1]-rijmag));
      wik = Sp(rikmag,rcmin[itype][ktype],rcmax[itype][ktype],dS);
      Nki = nC[atomk]-(wik*kronecker(itype,0))+nH[atomk] -
        (wik*kronecker(itype,1));
      cosjik = ((rij[0]*rik[0])+(rij[1]*rik[1])+(rij[2]*rik[2])) /
        (rijmag*rikmag);
      cosjik = MIN(cosjik,1.0);
      cosjik = MAX(cosjik,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosjik,(NijC+NijH),itype,&dgdc,&dgdN);
      Etmp = Etmp+(wik*g*exp(lamdajik));
      tmp3 = tmp3+(wik*dgdN*exp(lamdajik));
      NconjtmpI = NconjtmpI+(kronecker(ktype,0)*wik*Sp(Nki,Nmin,Nmax,dS)); //parameter Nmin Nmax
    }
  }

  PijS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  PijS = PijSpline(NijC,NijH,itype,jtype,dN2);
  pij = pow(1.0+Etmp+PijS,-0.5);
  tmp = -0.5*cube(pij);
  //printf("pij=%.16f Etmp=%.16f PijS=%.16f \n",pij,Etmp,PijS);
  // pij forces

  tmp = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  Etmp = 0.0;

  //REBO_neighs = REBO_firstneigh[j];
  
  for (l = 0; l < NBond[atomj]; l++) {
    atoml = Bonded[atomj*MaxNBond + l];
    if (atoml != atomi) {
      ltype = particle[atoml].type;
      rjl[0] = R_PBC(particle,atomj,atoml,0,L.x);
      rjl[1] = R_PBC(particle,atomj,atoml,1,L.y);
      rjl[2] = R_PBC(particle,atomj,atoml,2,L.z);
      rjlmag = sqrt((rjl[0]*rjl[0])+(rjl[1]*rjl[1])+(rjl[2]*rjl[2]));
      lamdaijl = 4.0*kronecker(jtype,1) *
        ((rho[ltype][1]-rjlmag)-(rho[itype][1]-rijmag));
      wjl = Sp(rjlmag,rcmin[jtype][ltype],rcmax[jtype][ltype],dS);
      Nlj = nC[atoml]-(wjl*kronecker(jtype,0)) +
        nH[atoml]-(wjl*kronecker(jtype,1));
      cosijl = -1.0*((rij[0]*rjl[0])+(rij[1]*rjl[1])+(rij[2]*rjl[2])) /
        (rijmag*rjlmag);
      cosijl = MIN(cosijl,1.0);
      cosijl = MAX(cosijl,-1.0);

      // evaluate splines g and derivatives dg

      g = gSpline(cosijl,NjiC+NjiH,jtype,&dgdc,&dgdN);
      Etmp = Etmp+(wjl*g*exp(lamdaijl));
      //printf("Etmp=%.16f wjl=%.16f g=%.16f lamdaijl=%.16f\n",Etmp,wjl,g,lamdaijl);
      tmp3 = tmp3+(wjl*dgdN*exp(lamdaijl));
      NconjtmpJ = NconjtmpJ+(kronecker(ltype,0)*wjl*Sp(Nlj,Nmin,Nmax,dS));
    }
  }

  PjiS = 0.0;
  dN2[0] = 0.0;
  dN2[1] = 0.0;
  
  PjiS = PijSpline(NjiC,NjiH,jtype,itype,dN2);
  
  pji = pow(1.0+Etmp+PjiS,-0.5);
  tmp = -0.5*cube(pji);
   
  //REBO_neighs = REBO_firstneigh[j];
 
  // evaluate Nij conj

  Nijconj = 1.0+(NconjtmpI*NconjtmpI)+(NconjtmpJ*NconjtmpJ);
  piRC = piRCSpline(NijC+NijH,NjiC+NjiH,Nijconj,itype,jtype,dN3);

  // piRC forces
  
  //REBO_neighs_i = REBO_firstneigh[i];

  Tij = 0.0;
  dN3[0] = 0.0;
  dN3[1] = 0.0;
  dN3[2] = 0.0;
  if (itype == 0 && jtype == 0)
    Tij=TijSpline((NijC+NijH),(NjiC+NjiH),Nijconj,dN3);
  Etmp = 0.0;

  if (fabs(Tij) > TOL) {
    atom2 = atomi;
    atom3 = atomj;
    r32[0] = R_PBC(particle,atom3,atom2,0,L.x);
    r32[1] = R_PBC(particle,atom3,atom2,1,L.y);
    r32[2] = R_PBC(particle,atom3,atom2,2,L.z);
    r32mag = sqrt((r32[0]*r32[0])+(r32[1]*r32[1])+(r32[2]*r32[2]));
    r23[0] = -r32[0];
    r23[1] = -r32[1];
    r23[2] = -r32[2];
    r23mag = r32mag;
    //REBO_neighs_i = REBO_firstneigh[i];
    for (k = 0; k < NBond[atomi]; k++) {
      atomk = Bonded[atomi*MaxNBond + k];
      atom1 = atomk;
      ktype = particle[atomk].type;
      if (atomk != atomj) {
	r21[0] = R_PBC(particle,atom2,atom1,0,L.x);
	r21[1] = R_PBC(particle,atom2,atom1,1,L.y);
	r21[2] = R_PBC(particle,atom2,atom1,2,L.z);
        r21mag = sqrt(r21[0]*r21[0] + r21[1]*r21[1] + r21[2]*r21[2]);
        cos321 = -1.0*((r21[0]*r32[0])+(r21[1]*r32[1])+(r21[2]*r32[2])) /
          (r21mag*r32mag);
        cos321 = MIN(cos321,1.0);
        cos321 = MAX(cos321,-1.0);
        Sp2(cos321,thmin,thmax,dcut321);
        sin321 = sqrt(1.0 - cos321*cos321);
        sink2i = 1.0/(sin321*sin321);
        rik2i = 1.0/(r21mag*r21mag);
        if (sin321 != 0.0) {
          rr = (r23mag*r23mag)-(r21mag*r21mag);
          rjk[0] = r21[0]-r23[0];
          rjk[1] = r21[1]-r23[1];
          rjk[2] = r21[2]-r23[2];
	  
          rjk2 = (rjk[0]*rjk[0])+(rjk[1]*rjk[1])+(rjk[2]*rjk[2]);
          rijrik = 2.0*r23mag*r21mag;
          rik2 = r21mag*r21mag;
          dctik = (-rr+rjk2)/(rijrik*rik2);
          dctij = (rr+rjk2)/(rijrik*r23mag*r23mag);
          dctjk = -2.0/rijrik;
          w21 = Sp(r21mag,rcmin[itype][ktype],rcmaxp[itype][ktype],dw21);
          rijmag = r32mag;
          rikmag = r21mag;
          rij2 = r32mag*r32mag;
          rik2 = r21mag*r21mag;
          costmp = 0.5*(rij2+rik2-rjk2)/rijmag/rikmag;
          tspjik = Sp2(costmp,thmin,thmax,dtsjik);
          dtsjik = -dtsjik;

          //REBO_neighs_j = REBO_firstneigh[j];
          for (l = 0; l < NBond[atomj]; l++) {
            //atoml = REBO_neighs_j[l];
	    atoml = Bonded[atomj*MaxNBond + l];
            atom4 = atoml;
            ltype = particle[atoml].type;
            if (!(atoml == atomi || atoml == atomk)) {
	      r34[0] = R_PBC(particle,atom3,atom4,0,L.x);
	      r34[1] = R_PBC(particle,atom3,atom4,1,L.y);
	      r34[2] = R_PBC(particle,atom3,atom4,2,L.z);
              r34mag = sqrt((r34[0]*r34[0])+(r34[1]*r34[1])+(r34[2]*r34[2]));
              cos234 = (r32[0]*r34[0] + r32[1]*r34[1] + r32[2]*r34[2]) /
                (r32mag*r34mag);
              cos234 = MIN(cos234,1.0);
              cos234 = MAX(cos234,-1.0);
              sin234 = sqrt(1.0 - cos234*cos234);
              sinl2i = 1.0/(sin234*sin234);
              rjl2i = 1.0/(r34mag*r34mag);

              if (sin234 != 0.0) {
                w34 = Sp(r34mag,rcmin[jtype][ltype],rcmaxp[jtype][ltype],dw34);
                rr = (r23mag*r23mag)-(r34mag*r34mag);
                //ril[0] = r23[0]+r34[0];
                //ril[1] = r23[1]+r34[1];
                //ril[2] = r23[2]+r34[2];
		ril[0] = R_PBC(particle,atomi,atoml,0,L.x);
		ril[1] = R_PBC(particle,atomi,atoml,1,L.y);
		ril[2] = R_PBC(particle,atomi,atoml,2,L.z);
                ril2 = (ril[0]*ril[0])+(ril[1]*ril[1])+(ril[2]*ril[2]);
                rijrjl = 2.0*r23mag*r34mag;
                rjl2 = r34mag*r34mag;
                dctjl = (-rr+ril2)/(rijrjl*rjl2);
                dctji = (rr+ril2)/(rijrjl*r23mag*r23mag);
                dctil = -2.0/rijrjl;
                rjlmag = r34mag;
                rjl2 = r34mag*r34mag;
                costmp = 0.5*(rij2+rjl2-ril2)/rijmag/rjlmag;
                tspijl = Sp2(costmp,thmin,thmax,dtsijl);
                dtsijl = -dtsijl;
                prefactor = VA*Tij;

                cross321[0] = (r32[1]*r21[2])-(r32[2]*r21[1]);
                cross321[1] = (r32[2]*r21[0])-(r32[0]*r21[2]);
                cross321[2] = (r32[0]*r21[1])-(r32[1]*r21[0]);
                cross234[0] = (r23[1]*r34[2])-(r23[2]*r34[1]);
                cross234[1] = (r23[2]*r34[0])-(r23[0]*r34[2]);
                cross234[2] = (r23[0]*r34[1])-(r23[1]*r34[0]);

                cwnum = (cross321[0]*cross234[0]) +
                  (cross321[1]*cross234[1]) + (cross321[2]*cross234[2]);
                cwnom = r21mag*r34mag*r23mag*r23mag*sin321*sin234;
                om1234 = cwnum/cwnom;
                cw = om1234;
                Etmp += ((1.0-square(om1234))*w21*w34) *
                  (1.0-tspjik)*(1.0-tspijl);
		
	      }
            }
          }
	  
        }
      }
    }
  }
  bij = (0.5*(pij+pji))+piRC+(Tij*Etmp);
  
  return bij;
  
  
}
double TORSION_pot(atom *particle,int *BondNum,int Natoms,point L,int *Bonded,int *NBond,double &E_TORSION_,int *NBond_once,int *Bonded_once)
{
  double E_TORSION;
  E_TORSION=0;
#pragma omp parallel
  {
    int i,j,k,l,ii,inum;
    //tagint itag,jtag;
    int itag,jtag;
    double evdwl,fpair,xtmp,ytmp,ztmp;
    double cos321;
    double w21,dw21,cos234,w34,dw34;
    double cross321[3],cross321mag,cross234[3],cross234mag;
    double w23,dw23,cw2,ekijl,Ec;
    double cw,cwnum,cwnom;
    double rij,rij2,rik,rjl,tspjik,dtsjik,tspijl,dtsijl,costmp,fcpc;
    double sin321,sin234,rjk2,rik2,ril2,rjl2;
    double rjk,ril;
    double Vtors;
    double dndij[3],tmpvec[3],dndik[3],dndjl[3];
    double dcidij,dcidik,dcidjk,dcjdji,dcjdjl,dcjdil;
    double dsidij,dsidik,dsidjk,dsjdji,dsjdjl,dsjdil;
    double dxidij,dxidik,dxidjk,dxjdji,dxjdjl,dxjdil;
    double ddndij,ddndik,ddndjk,ddndjl,ddndil,dcwddn,dcwdn,dvpdcw,Ftmp[3];
    double del32[3],rsq,r32,del23[3],del21[3],r21;
    double deljk[3],del34[3],delil[3],delkl[3],r23,r34;
    double fi[3],fj[3],fk[3],fl[3];
    int itype,jtype,ktype,ltype,kk,ll,jj;
    int *ilist,*REBO_neighs_i,*REBO_neighs_j;
    
#pragma omp for reduction(+:E_TORSION)
    for (ii = 0; ii < Natoms; ii++) {
      i=ii;
      
      itag = ii;
      itype = particle[i].type;
      if (itype != 0) continue;
      xtmp = particle[i].r[0];
      ytmp = particle[i].r[1];
      ztmp = particle[i].r[2];
      
      for (jj = 0; jj < NBond_once[ii]; jj++) {
	
	j = Bonded_once[i*MaxNBond + jj];
	jtag = j;
	
	jtype = particle[j].type;
	if (jtype != 0) continue;
	
	del32[0] = R_PBC(particle,j,i,0,L.x);
	del32[1]  = R_PBC(particle,j,i,1,L.y);
	del32[2]  = R_PBC(particle,j,i,2,L.z);
	rsq = del32[0]*del32[0] + del32[1]*del32[1] + del32[2]*del32[2];
	r32 = sqrt(rsq);
	del23[0] = -del32[0];
	del23[1] = -del32[1];
	del23[2] = -del32[2];
	r23 = r32;
	w23 = Sp(r23,rcmin[itype][jtype],rcmax[itype][jtype],dw23);
	//cerr<<"hey3"<<endl;
	for (kk = 0; kk < NBond[i]; kk++) {
	  //k = REBO_neighs_i[kk];
	  k = Bonded[i*MaxNBond + kk];
	  ktype = particle[k].type;
	  if (k == j) continue;
	  del21[0] = R_PBC(particle,i,k,0,L.x);
	  del21[1] = R_PBC(particle,i,k,1,L.y);
	  del21[2] = R_PBC(particle,i,k,2,L.z);
	  rsq = del21[0]*del21[0] + del21[1]*del21[1] + del21[2]*del21[2];
	  r21 = sqrt(rsq);
	  cos321 = - ((del21[0]*del32[0]) + (del21[1]*del32[1]) +
		      (del21[2]*del32[2])) / (r21*r32);
	  cos321 = MIN(cos321,1.0);
	  cos321 = MAX(cos321,-1.0);
	  sin321 = sqrt(1.0 - cos321*cos321);
	  if (sin321 < TOL) continue;
	  
	  //deljk[0] = del21[0]-del23[0];
	  //deljk[1] = del21[1]-del23[1];
	  //deljk[2] = del21[2]-del23[2];
	  
	  deljk[0] = R_PBC(particle,j,k,0,L.x);
	  deljk[1] = R_PBC(particle,j,k,1,L.y);
	  deljk[2] = R_PBC(particle,j,k,2,L.z);
	  
	  rjk2 = deljk[0]*deljk[0] + deljk[1]*deljk[1] + deljk[2]*deljk[2];
	  rjk=sqrt(rjk2);
	  rik2 = r21*r21;
	  w21 = Sp(r21,rcmin[itype][ktype],rcmax[itype][ktype],dw21);
	  
	  rij = r32;
	  rik = r21;
	  rij2 = r32*r32;
	  rik2 = r21*r21;
	  costmp = 0.5*(rij2+rik2-rjk2)/rij/rik;
	  tspjik = Sp2(costmp,thmin,thmax,dtsjik);
	  dtsjik = -dtsjik;
	  //REBO_neighs_j = REBO_firstneigh[j];
	  for (ll = 0; ll <NBond[j]; ll++) {
	    //l = REBO_neighs_j[ll];
	    l = Bonded[j*MaxNBond + ll];
	    ltype = particle[l].type;
	    if (l == i || l == k) continue;
	    del34[0] = R_PBC(particle,j,l,0,L.x);
	    del34[1] = R_PBC(particle,j,l,1,L.y);
	    del34[2] = R_PBC(particle,j,l,2,L.z);
	    rsq = del34[0]*del34[0] + del34[1]*del34[1] + del34[2]*del34[2];
	    r34 = sqrt(rsq);
	    cos234 = (del32[0]*del34[0] + del32[1]*del34[1] +
		      del32[2]*del34[2]) / (r32*r34);
	    cos234 = MIN(cos234,1.0);
	    cos234 = MAX(cos234,-1.0);
	    sin234 = sqrt(1.0 - cos234*cos234);
	    if (sin234 < TOL) continue;
	    w34 = Sp(r34,rcmin[jtype][ltype],rcmax[jtype][ltype],dw34);
	    
	    delil[0] = R_PBC(particle,i,l,0,L.x);
	    delil[1] = R_PBC(particle,i,l,1,L.y);
	    delil[2] = R_PBC(particle,i,l,2,L.z);

	    //delil[0] = del23[0] + del34[0];
	    //delil[1] = del23[1] + del34[1];
	    //delil[2] = del23[2] + del34[2];
	    
	    ril2 = delil[0]*delil[0] + delil[1]*delil[1] + delil[2]*delil[2];
	    ril=sqrt(ril2);
	    rjl2 = r34*r34;
	    
	    rjl = r34;
	    rjl2 = r34*r34;
	    costmp = 0.5*(rij2+rjl2-ril2)/rij/rjl;
	    tspijl = Sp2(costmp,thmin,thmax,dtsijl);
	    dtsijl = -dtsijl; //need minus sign
	    cross321[0] = (del32[1]*del21[2])-(del32[2]*del21[1]);
	    cross321[1] = (del32[2]*del21[0])-(del32[0]*del21[2]);
	    cross321[2] = (del32[0]*del21[1])-(del32[1]*del21[0]);
	    cross321mag = sqrt(cross321[0]*cross321[0]+
			       cross321[1]*cross321[1]+
			       cross321[2]*cross321[2]);
	    cross234[0] = (del23[1]*del34[2])-(del23[2]*del34[1]);
	    cross234[1] = (del23[2]*del34[0])-(del23[0]*del34[2]);
	    cross234[2] = (del23[0]*del34[1])-(del23[1]*del34[0]);
	    cross234mag = sqrt(cross234[0]*cross234[0]+
			       cross234[1]*cross234[1]+
			       cross234[2]*cross234[2]);
	    cwnum = (cross321[0]*cross234[0]) +
	      (cross321[1]*cross234[1])+(cross321[2]*cross234[2]);
	    cwnom = r21*r34*r32*r32*sin321*sin234;
	    cw = cwnum/cwnom;
	    
	    cw2 = (.5*(1.0-cw));
	    ekijl = epsilonT[ktype][ltype];
	    Ec = 256.0*ekijl/405.0;
	    Vtors = (Ec*(powint(cw2,5)))-(ekijl/10.0);
	    
	    //#pragma omp atomic
	    E_TORSION += Vtors*w21*w23*w34*(1.0-tspjik)*(1.0-tspijl);
	    
	  }
	}
      }
    }
  }
  return(E_TORSION);
}
double calc_F_Rep_pot(int Natoms, atom *particle, double *rijShiled_1, BondParamsStruct *BondParams, point L,Normal_struct *Normal,int *Normal_atom,int Interlayer,int *Neighb,int *NNeighb,dr_ij *dN,double &ERep_,double &ERepulsion,double &EDisp,double &RI,double &Interlayer_dis)
{
  
  double ERep=0.0;
  double SCC,SCC1,SCC2, Pavg;
  int Nlayer1,Nlayer2,atomi;
  
  Nlayer1=0;
  Nlayer2=0;
  
  RI=EDisp=Interlayer_dis=0;
  
  for(atomi=0; atomi < Natoms; atomi++){
    if(particle[atomi].layer==0)Nlayer1++;
    if(particle[atomi].layer==1)Nlayer2++;
  }

  //cerr<<"Nlayer1= "<<Nlayer1<<endl;
  //cerr<<"Nlayer2= "<<Nlayer2<<endl;
  
  SCC=SCC1=SCC2=0.0;
  
  if(Nlayer1 > Nlayer2){
    SCC2=Nlayer2*CircleOverlap(0,rC,rC);
    SCC1=0.5*Nlayer2*CircleOverlap(0,rC,rC);
  }
  else{
    SCC2=Nlayer1*CircleOverlap(0,rC,rC);
    SCC1=0.5*Nlayer1*CircleOverlap(0,rC,rC);
  }
  
  EDisp = ERepulsion = 0.0;
  //#pragma omp parallel
  //{
    int atomk, typek, typei,atomj,typej, index,j;
    double alpha, rvdW, Epsilon, term1, dEvdW_drij;
    double dTap, r_ij, Tap,Epsilon_C_vdW_exp,term3,term4,term5,term7;
    double dP_ij_drk,dP_ji_drk,R,Pij,Pji,gamma,dr_ij_drk,exp_term1;
    double exp_Pij_gamma, exp_Pji_gamma,Pij_gamma,Pji_gamma,reff,r_ik,r_jk;
    double rx_ij,ry_ij,rz_ij,rx_ji,ry_ji,rz_ji,inside_acos_ij,inside_acos_ji,Pij_term1,Pji_term1,Pij_term,Pji_term,dr_ji_drk;
    double angle_ij, sqrt_inside_ij, sin_angle_ij,cos_angle_ij,angle_ji,sqrt_inside_ji,sin_angle_ji,cos_angle_ji,dinside_acos,dangle;
    double C0_kc,C2_kc,C4_kc,A_kc,fji,fij,term6,dfji,dfij,A_term;
    double SCC_tmp,Interlayer_dis_Tmp;
    int n,dir;
    dP_ij_drk = 0.0;
    dP_ji_drk = 0.0;
    
    C0_kc      = 15.71/1000;
    C2_kc      = 12.29/1000;
    C4_kc      = 4.933/1000;
    alpha   = 3.629;
    Epsilon = 3.03/1000;
    R       = 3.34;
    gamma   = 0.578;
    A_kc    = 10.238/1000;
    printf("%i \n\n",Natoms);
    //#pragma omp for schedule(dynamic,15) reduction(+:ERep)
    for(atomi=0 ; atomi < Natoms ; atomi++){
      typei = particle[atomi].type;
      SCC_tmp=0.0;
      for (j=0 ; j < NNeighb[atomi]; j++){
	atomk = atomj = Neighb[atomi*MaxNeighb + j] ;	  
	typek = typej = particle[atomj].type;	
	//cerr<<"hey2"<<endl;
	index = typej*NAtomParams + typei;
	//r_ij  = R_ij(particle,atomi,atomj, L);
	r_ij = rijShiled_1[atomi*MaxNeighb + j];
	Tap   = calc_Tap(r_ij);
	
	reff    = BondParams[index].r_eff;	    
	
	Pij_term1=0.0;
	Pji_term1=0.0;
	dN[0].r[0]=0;
	dN[0].r[1]=0;
	dN[0].r[2]=0;
	rx_ij = particle[atomi].r[0] - particle[atomj].r[0];
	ry_ij = particle[atomi].r[1] - particle[atomj].r[1];
	rz_ij = particle[atomi].r[2] - particle[atomj].r[2];
	rx_ij -= L.x * rint(rx_ij / (L.x + TINY));
	ry_ij -= L.y * rint(ry_ij / (L.y + TINY));
	rz_ij -= L.z * rint(rz_ij / (L.z + TINY));
	
	rx_ji = -rx_ij;
	ry_ji = -ry_ij;
	rz_ji = -rz_ij;
	
	inside_acos_ij = (rx_ij*Normal[atomi].x + ry_ij*Normal[atomi].y + rz_ij*Normal[atomi].z) / (r_ij);
	inside_acos_ji = (rx_ji*Normal[atomj].x + ry_ji*Normal[atomj].y + rz_ji*Normal[atomj].z) / (r_ij);
	Pij_term1 = inside_acos_ij/r_ij;
	Pji_term1 = inside_acos_ji/r_ij;
	
	angle_ij  = acos(inside_acos_ij);
	sqrt_inside_ij=sqrt(1.0-sqr(inside_acos_ij));
	sin_angle_ij=sin(angle_ij);
	cos_angle_ij=cos(angle_ij);
	angle_ji  = acos(inside_acos_ji);
	sqrt_inside_ji=sqrt(1.0-sqr(inside_acos_ji));
	sin_angle_ji=sin(angle_ji);
	cos_angle_ji=cos(angle_ji);
	
	if (fabs(inside_acos_ij) < 1.0){
	  Pij=r_ij*sin(angle_ij);
	}
	else Pij=0;
	
	if (fabs(inside_acos_ji) < 1.0){
	  Pji=r_ij*sin(angle_ji);
	}
	else Pji=0;
	
	Pavg=(Pij+Pji)/2.0;
	//#pragma omp critical
	SCC += CircleOverlap(Pavg,rC,rC);
	SCC_tmp += CircleOverlap(Pavg,rC,rC);
	
	if(r_ij < 4.0 ){
	  Interlayer_dis_Tmp = (sqrt(r_ij*r_ij -  Pij*Pij ) + sqrt(r_ij*r_ij -  Pji*Pji ))/2;
	  //if((atomi==552 || atomj==552)){
	  if((atomi==0 || atomj==0))
	    Interlayer_dis     = (sqrt(r_ij*r_ij -  Pij*Pij ) + sqrt(r_ij*r_ij -  Pji*Pji ))/2;
	  //cerr<<"Interlayer_dis= "<<Interlayer_dis<<endl;
	}
	
	
	Pij_gamma = Pij/gamma; //new
	Pji_gamma = Pji/gamma; //new
	exp_Pij_gamma = exp(-pow((Pij_gamma),2.0)); //new
	exp_Pji_gamma = exp(-pow((Pji_gamma),2.0)); //new
	exp_term1 = exp(-alpha*(r_ij-R));
	
	fij=0.0;
	fji=0.0;
	
	fij = C0_kc+ C2_kc*pow((Pij_gamma),2.0)+ C4_kc*pow((Pij_gamma),4.0);
	fji = C0_kc+ C2_kc*pow((Pji_gamma),2.0)+ C4_kc*pow((Pji_gamma),4.0);
	
	fij*=exp_Pij_gamma;
	fji*=exp_Pji_gamma;
	
	A_term=A_kc*pow((r_ij/R),-6);
	
	Epsilon_C_vdW_exp = Epsilon + fij + fji;
	//
	ERep += Tap*(exp_term1*(Epsilon_C_vdW_exp) - A_term);
	//#pragma omp atomic
	ERepulsion += Tap*(exp_term1*(Epsilon_C_vdW_exp));
	//#pragma omp atomic
	EDisp += Tap*(-A_term);
	
	
	//printf("atomi= %i atomj= %i ERep= %.16f \n",atomi,atomj,Tap*(exp_term1*(Epsilon_C_vdW_exp) - A_term));
      }
      //if(particle[atomi].layer==1)
      //printf("C %.16f %.16f %.16f %.16f \n",particle[atomi].r[0],particle[atomi].r[1],particle[atomi].r[2],(SCC_tmp)/(0.5*CircleOverlap(0,rC,rC))-1.0 );
      particle[atomi].RI    = (SCC_tmp)/(0.5*CircleOverlap(0,rC,rC))-1.0 ;
      particle[atomi].Inter = 2.0*(Interlayer_dis_Tmp-3.2)/(3.7-3.2)-1;
      
    }
    // }
  RI=(SCC-SCC1)/(SCC2-SCC1);
  return(ERep);
}
