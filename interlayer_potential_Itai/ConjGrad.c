#include "declarations.h"

// *************************************************************************************

inline void SHIFT( double &a , double &b , double &c , double &d ) {(a)=(b); (b)=(c); (c)=(d);}
inline void MOV3(double &a, double &b, double &c, double &d, double &e1, double &f) {(a)=(d); (b)=(e1); (c)=(f);}

// *************************************************************************************
// This function calculates the value of the potential at a point
// at distance x times the direction vector from the starting point.

void f1D(double &V1D, double &x, atom *temp_particle, point *direction,int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams,int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *angle, double *Bmat, double *qvec, double *Avec, double *rijShiled_1, FILE* logfile,point *Normal,int *Normal_atom,int Interlayer, double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double Xmin,double Ymin,double Zmin,int *NList,int *List,point *R0)
{
  int atomi;
  
#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    temp_particle[atomi].r[0]    = particle[atomi].r[0] + x * direction[atomi].x;
    temp_particle[atomi].r[1]    = particle[atomi].r[1] + x * direction[atomi].y;
    temp_particle[atomi].r[2]    = particle[atomi].r[2] + x * direction[atomi].z;

    temp_particle[atomi].r[0] -= L.x * floor((temp_particle[atomi].r[0] - Xmin) / (L.x + TINY));
    temp_particle[atomi].r[1] -= L.y * floor((temp_particle[atomi].r[1] - Ymin) / (L.y + TINY));
    temp_particle[atomi].r[2] -= L.z * floor((temp_particle[atomi].r[2] - Zmin) / (L.z + TINY));
    
    temp_particle[atomi].type    = particle[atomi].type;
    temp_particle[atomi].layer   = particle[atomi].layer;
    temp_particle[atomi].Mass    = particle[atomi].Mass;
    temp_particle[atomi].vx      = particle[atomi].vx = 0.0;
    temp_particle[atomi].vy      = particle[atomi].vy = 0.0;
    temp_particle[atomi].vz      = particle[atomi].vz = 0.0;
  }

  calc_force_arrays(Natoms,NBond,Bonded,temp_particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
  
  if (Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  V1D = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,temp_particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
}

// *************************************************************************************
// This function calculates the value of the force and the potential at a point
// at distance x times the direction vector from the starting point.

void df1D(double &dV, double &V1D, double &x, atom *temp_particle, point *direction, point *_dV, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,point *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double Xmin,double Ymin,double Zmin,int *NList,int *List,point *R0)
{
  int atomi;
  
  dV = 0.0;
  
#pragma omp for 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    temp_particle[atomi].r[0]    = particle[atomi].r[0] + x * direction[atomi].x;
    temp_particle[atomi].r[1]    = particle[atomi].r[1] + x * direction[atomi].y;
    temp_particle[atomi].r[2]    = particle[atomi].r[2] + x * direction[atomi].z;
    
    temp_particle[atomi].r[0] -= L.x * floor((temp_particle[atomi].r[0] - Xmin) / (L.x + TINY));
    temp_particle[atomi].r[1] -= L.y * floor((temp_particle[atomi].r[1] - Ymin) / (L.y + TINY));
    temp_particle[atomi].r[2] -= L.z * floor((temp_particle[atomi].r[2] - Zmin) / (L.z + TINY));
    
    temp_particle[atomi].type    = particle[atomi].type;
    temp_particle[atomi].layer   = particle[atomi].layer;
    temp_particle[atomi].Mass    = particle[atomi].Mass;
    temp_particle[atomi].vx      = particle[atomi].vx = 0.0;
    temp_particle[atomi].vy      = particle[atomi].vy = 0.0;
    temp_particle[atomi].vz      = particle[atomi].vz = 0.0;
  }
    
  calc_force_arrays(Natoms,NBond,Bonded,temp_particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
  
  if (Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  V1D = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,temp_particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
  
  Calc_Force(_dV,Natoms,NBond,Bonded,temp_particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag);
  
  //#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++)
    //#pragma omp atomic
    dV -= (temp_particle[atomi].Fx * direction[atomi].x + temp_particle[atomi].Fy * direction[atomi].y + temp_particle[atomi].Fz * direction[atomi].z);
  
}

// *************************************************************************************

double dbrent(double &a1, double &b1, double &c1, double tol, double &xmin, atom *temp_particle, point *direction, point *_dV, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams,int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,int &found,point *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double Xmin,double Ymin,double Zmin,int *NList,int *List,point *R0)
{
  int iter, ok1, ok2;
  double a,b,d,e1=0.0,olde,x,w,v,xm,d1,d2;
  double fx,fw,fv,fu,dx,dw,dv,du;
  double tol1,tol2;
  double u,u1,u2;
  double temp;
  
  a = (a1 < c1 ? a1 : c1);
  b = (a1 > c1 ? a1 : c1);

  x=w=v=b1;

  df1D(dx,fx,x,temp_particle,direction,_dV,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
  
  fw=fv=fx;
  dw=dv=dx;

  for(iter=0 ; iter < ITMAXdbrent ; iter++){
    xm = 0.5 * (a + b);
    tol1 = tol * fabs(x) + ZEPS;
    tol2 = 2.0 * tol1;

    if(fabs(x-xm) <= ( tol2 - 0.5*(b-a) )){
      found=1;
      xmin=x;
      return(fx);
    }
    
    if(fabs(e1) > tol1){
      d1=2.0*(b-a);
      d2=d1;

      if(dw != dx) d1=(w-x)*dx/(dx-dw);
      if(dv != dx) d2=(v-x)*dx/(dx-dv);

      u1=x+d1;
      u2=x+d2;

      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;

      olde=e1;
      e1=d;
      
      if(ok1 || ok2){
	if(ok1 && ok2) d = (fabs(d1) < fabs(d2) ? d1 : d2 );
	else if(ok1)   d = d1;
	else           d = d2;
	
	if(fabs(d) <= fabs(0.5*olde)){
	  u = x + d;
	  if(((u-a) < tol2) || ((b-u) < tol2)) d = SIGN(tol1,xm-x);
	}
	else d = 0.5 * ( e1 = ( dx >= 0.0 ? a-x : b-x ) );
      }
      else d = 0.5 * ( e1 = ( dx >= 0.0 ? a-x : b-x ) );
    }
    else d = 0.5 * ( e1 = ( dx >= 0.0 ? a-x : b-x ) );
    
    if(fabs(d) >= tol1){
      u = x + d;
      
      f1D(fu,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
    }
    else{
      u = x + SIGN(tol1,d);

      f1D(fu,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
      if(fu > fx){
	found=1;
	xmin = x;
	return fx;
      }
    }
    
    df1D(du,temp,u,temp_particle,direction,_dV,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
    
    if(fu <= fx){
      if(u >= x) a = x; else b = x;
      MOV3(v,fv,dv,w,fw,dw);
      MOV3(w,fw,dw,x,fx,dx);
      MOV3(x,fx,dx,u,fu,du);
    }
    else{
      if(u < x) a = u; else b = u;
      if(fu <= fw || w == x){
	MOV3(v,fv,dv,w,fw,dw);
	MOV3(w,fw,dw,u,fu,du);
      }
      else if((fu < fv) || (v == x) || (v == w)) MOV3(v,fv,dv,u,fu,du);
    }
  }
  cerr<<"Exceeded maximum number of iterations in routine dbrent\n";
  exit(0);
}

// *************************************************************************************

// This function initially brackets the closest minimum

void min_bracket(double &a, double &b, double &c, double &f_a, double &f_b, double &f_c, atom *temp_particle, point *direction,int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,point *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double Xmin,double Ymin,double Zmin,int *NList,int *List,point *R0)
{
  double temp, u, r, q, ulim, f_u;
  
  f1D(f_a,a,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
  f1D(f_b,b,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
  if(f_b > f_a){
    SHIFT(temp,a,b,temp);
    SHIFT(temp,f_b,f_a,temp);
  }
  
  c = b + GOLD * (b - a);

  f1D(f_c,c,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
  
  while(f_b > f_c){
    r = (b - a) * (f_b - f_c);
    q = (b - c) * (f_b - f_a);
    u = b - ((b-c) * q - (b-a) * r)/(2.0 * SIGN(MAX(fabs(q-r) , TINY) , q-r));
    ulim = b + GLIMIT * (c - b);

    if((b-u)*(u-c) > 0.0){
      f1D(f_u,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
      if(f_u < f_c){
	a = b;
	b = u;
	f_a = f_b;
	f_b = f_u;
	return;
      }
      else if(f_u > f_b){
	c = u;
	f_c = f_u;
	return;
      }
      
      u = c + GOLD * (c - b);
      f1D(f_u,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
    }
    
    else if((c-u)*(u-ulim) > 0.0){
      f1D(f_u,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
      
      if(f_u < f_c){
	temp = c+GOLD*(c-b);
	SHIFT(b,c,u,temp);
	
	f1D(temp,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
	
	SHIFT(f_b,f_c,f_u,temp);
      }
    }
    
    else if((u-ulim)*(ulim-c) >= 0.0){
      u = ulim;
      
      f1D(f_u,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
    }
    
    else{
      u = c + GOLD * (c-b);
      
      f1D(f_u,u,temp_particle,direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Bmat,qvec,Avec,rijShiled_1,logfile,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);  
    }
    
    SHIFT(a,b,c,u);
    SHIFT(f_a,f_b,f_c,f_u);
  }
}

// *************************************************************************************

// This function finds a minimum of the potential on a line in the conj_direction
// that passes through the present coordinates of the atoms.

#define TOL 2.0e-4

void optimize_1D(atom *temp_particle, point *conj_direction, double &min_V, point *_dV, int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double *angle,point *Normal,int *Normal_atom,int Interlayer,double *Fc,int *Neighb,int *NNeighb,int Coulomb_flag,double Xmin,double Ymin,double Zmin,int *NList,int *List,point *R0)
{
  int j=0, found=0, atomj;
  double xx,xmin,fx,fb,fa,bx,ax;
  
  do{
    ax = 0.0;
    xx = ++j;
    
    min_bracket(ax,xx,bx,fa,fx,fb,temp_particle,conj_direction,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
    
    min_V = dbrent(ax,xx,bx,TOL,xmin,temp_particle,conj_direction,_dV,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,found,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
    
  }while(!found);
  
#pragma omp parallel for firstprivate(xmin) 
  for(atomj=0 ; atomj < Natoms ; atomj++){
    //#pragma omp critical
    conj_direction[atomj].x *= xmin;
    //#pragma omp critical
    conj_direction[atomj].y *= xmin;
    //#pragma omp critical
    conj_direction[atomj].z *= xmin;
    //#pragma omp critical
    particle[atomj].r[0] += conj_direction[atomj].x;
    //#pragma omp critical
    particle[atomj].r[1] += conj_direction[atomj].y;
    //#pragma omp critical
    particle[atomj].r[2] += conj_direction[atomj].z;
    
    particle[atomj].r[0] -= L.x * floor((particle[atomj].r[0] - Xmin) / (L.x + TINY));
    particle[atomj].r[1] -= L.y * floor((particle[atomj].r[1] - Ymin) / (L.y + TINY));
    particle[atomj].r[2] -= L.z * floor((particle[atomj].r[2] - Zmin) / (L.z + TINY));
    
  }
}

#undef TOL

// *************************************************************************************
void conj_grad(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams,double *angle, double OptTol, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,int NConserve,int NTraj,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,point *Normal,int *Normal_atom,int Interlayer,double *Fc,int Periodic_flag,int *Neighb,int *NNeighb,int Coulomb_flag,int *NList,int *List,point *R0)
{  
  char optimization_file_name[50];
  int atomi, frame_counter, file_counter, precentage, iter;
  double V, V_init, gama, gama_mone, gama_mechane, min_V;
  
  frame_counter = 0;
  file_counter  = 1;
  //precentage  = 0;
  
  point *_grad_V        = new point[Natoms];
  point *_grad_V_new    = new point[Natoms];
  point *conj_direction = new point[Natoms];
  point *_dV            = new point[Natoms];
  atom  *temp_particle  = new atom[Natoms];
  
  // Initialization ************************************************************************************
  
  if(Periodic_flag == 1)    
    relax_cell(Natoms,NBond,AtomParams,BondParams,Bonded,Bmat,qvec,Avec,rij,rijShiled_1,L,particle,Normal,Fc,logfile,angle,BondNum,Neighb,NNeighb,Coulomb_flag,Interlayer,Normal_atom,NList,List,R0);
  
  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
  
  if(Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  V = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
  
  Calc_Force(_grad_V,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag);
  
  V_init = V;
  
#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    _grad_V[atomi].x = particle[atomi].Fx;
    _grad_V[atomi].y = particle[atomi].Fy;
    _grad_V[atomi].z = particle[atomi].Fz;
    
    conj_direction[atomi].x = _grad_V[atomi].x;
    conj_direction[atomi].y = _grad_V[atomi].y;
    conj_direction[atomi].z = _grad_V[atomi].z;
  }
  
  FILE* energy_file = fopen("E.dat" , "w");
  
  if(!energy_file){
    cerr<<"cannot open optimization energy file for writing\n";
    exit(0);
  }
  
  FILE* path_file = fopen("Traj-1.xyz" , "w");
  
  if(!path_file){
    cerr<<"cannot open optimization path file for writing\n";
    exit(0);
  }
  
  fprintf(energy_file,"%i %.16f\n",0,V_init);
  fflush(energy_file);
  
  Print_Traj(path_file,Natoms,0,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
  frame_counter ++;
  
  // Starting optimization *****************************************************************************
  
  //cerr<<"\n\nBeginning conjugate gradient algorithem\n\n";
  //cerr<<"0%                                            100%\n";
  
  for(iter=1 ; iter < ITMAXCG  ; iter++){
    
    if(Periodic_flag == 1)    
      relax_cell(Natoms,NBond,AtomParams,BondParams,Bonded,Bmat,qvec,Avec,rij,rijShiled_1,L,particle,Normal,Fc,logfile,angle,BondNum,Neighb,NNeighb,Coulomb_flag,Interlayer,Normal_atom,NList,List,R0);
    
    optimize_1D(temp_particle,conj_direction,min_V,_dV,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
    
    if(!(iter%NConserve)){
      fprintf(energy_file,"%i %.16f\n",iter,min_V);
      fflush(energy_file);
    }
    
    if(!(iter%NTraj)){
      if (!(frame_counter % Frame_limit)){
	fclose(path_file);
	file_counter ++;
	sprintf(optimization_file_name , "Traj-%d.xyz" , file_counter);
	path_file = fopen(optimization_file_name , "w");
	
	if(!path_file){
	  cerr<<"cannot open optimization file for writing\n";
	  exit(0);
	}
      }
      
      Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      frame_counter ++;
    }
    
    if( (2.0 * fabs( min_V - V )) <= (OptTol * ( fabs(min_V) + fabs(V) + EPS )) ){
      fprintf(energy_file,"%i %.16f\n",iter,min_V);
      fflush(energy_file);
      
      Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      Print_Traj_gv(Natoms,particle);
      frame_counter ++;
      
      delete [] _grad_V;
      delete [] _grad_V_new;
      delete [] conj_direction;
      delete [] temp_particle;
      delete [] _dV;
      
      cerr<<"Completed optimization after: "<<iter<<" iterations.\n\n";
      cerr<<"The initial value of the potential was: "<<V_init<<" and the final is: "<<min_V<<" (1).\n";
      fclose(path_file);
      fclose(energy_file);
      return;
    }

    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
    
    if(Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    V = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
    
    Calc_Force(_grad_V_new,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag);
    
    gama_mone = gama_mechane = 0.0;
    
    
#pragma omp parallel for reduction(+:gama_mechane) reduction(+:gama_mone)
    for(atomi=0 ; atomi < Natoms ; atomi++){
      _grad_V_new[atomi].x = particle[atomi].Fx;
      _grad_V_new[atomi].y = particle[atomi].Fy;
      _grad_V_new[atomi].z = particle[atomi].Fz;
      
      //#pragma omp atomic
      gama_mechane += (_grad_V[atomi].x * _grad_V[atomi].x + _grad_V[atomi].y * _grad_V[atomi].y + _grad_V[atomi].z * _grad_V[atomi].z);
      //#pragma omp atomic
      gama_mone += (
		    ( ( _grad_V_new[atomi].x - _grad_V[atomi].x ) * _grad_V_new[atomi].x ) +
		    ( ( _grad_V_new[atomi].y - _grad_V[atomi].y ) * _grad_V_new[atomi].y ) +
		    ( ( _grad_V_new[atomi].z - _grad_V[atomi].z ) * _grad_V_new[atomi].z )
		    );
    }
    
    if(gama_mechane == 0.0){
      fprintf(energy_file,"%i %.16f\n",iter,min_V);
      fflush(energy_file);
      
      Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      Print_Traj_gv(Natoms,particle);
      frame_counter ++;
      
      delete [] _grad_V;
      delete [] _grad_V_new;
      delete [] conj_direction;
      delete [] temp_particle;
      delete [] _dV;
      
      cerr<<"Completed optimization after: "<<iter<<" iterations.\n\n";
      cerr<<"The initial value of the potential was: "<<V_init<<" and the final is: "<<min_V<<" (2).\n";
      fclose(path_file);
      fclose(energy_file);
      return;
    }
    
    gama = gama_mone / gama_mechane ;
    
#pragma omp parallel for firstprivate(gama) 
    for(atomi=0 ; atomi < Natoms ; atomi++){
      _grad_V[atomi].x = _grad_V_new[atomi].x;
      _grad_V[atomi].y = _grad_V_new[atomi].y;
      _grad_V[atomi].z = _grad_V_new[atomi].z;
      //#pragma omp critical
      conj_direction[atomi].x *= gama;
      //#pragma omp critical
      conj_direction[atomi].y *= gama;
      //#pragma omp critical
      conj_direction[atomi].z *= gama;
      //#pragma omp critical
      conj_direction[atomi].x += _grad_V_new[atomi].x;
      //#pragma omp critical
      conj_direction[atomi].y += _grad_V_new[atomi].y;
      //#pragma omp critical
      conj_direction[atomi].z += _grad_V_new[atomi].z;
    }
    
    /*
      if((double(iter)/ITMAXCG)*100 > precentage){
      cerr<<'*';
      precentage += 2;
      }
    */
  }
  
  fprintf(energy_file,"%i %.16f\n",iter,min_V);
  fflush(energy_file);
  
  Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
  Print_Traj_gv(Natoms,particle);
  frame_counter ++;
  
  delete [] _grad_V;
  delete [] _grad_V_new;
  delete [] conj_direction;
  delete [] temp_particle;
  delete [] _dV;
    
  cerr<<"Exceeded Max. number of iterations for optimization\n";
  fclose(path_file);
  fclose(energy_file);
  
  exit(0);
  }

// *************************************************************************************

void steepest_descent(int Natoms, int *NBond, int *Bonded, atom *particle, BondParamsStruct *BondParams, int *BondNum, double *rij, point L, AtomParamsStruct *AtomParams,double *angle, int NTraj, int NConserve, double OptTol, double *Bmat, double *qvec, double *Avec, double *dqvec, double *rijShiled_1, FILE* logfile,double Xmin,double Ymin,double Zmin,double Xmax,double Ymax,double Zmax,point *Normal,int *Normal_atom,int Interlayer,double *Fc,int Periodic_flag,int *Neighb,int *NNeighb,int Coulomb_flag,int *NList,int *List,point *R0)
{
  char optimization_file_name[50];
  int atomi, iter, frame_counter, file_counter, precentage;
  double V, V_init, _grad_V_2, min_V;

  atom *temp_particle = new atom[Natoms];
  point *_grad_V      = new point[Natoms];
  point *_dV          = new point[Natoms];
  
  // Initialization ************************************************************************************
  
  frame_counter = 0;
  file_counter  = 1;
  //precentage  = 0;
  
  if(Periodic_flag == 1)    
    relax_cell(Natoms,NBond,AtomParams,BondParams,Bonded,Bmat,qvec,Avec,rij,rijShiled_1,L,particle,Normal,Fc,logfile,angle,BondNum,Neighb,NNeighb,Coulomb_flag,Interlayer,Normal_atom,NList,List,R0);
  
  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
  
  if(Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  
  V = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
  
  Calc_Force(_grad_V,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag);
  
  
#pragma omp parallel for 
  for(atomi=0 ; atomi < Natoms ; atomi++){
    _grad_V[atomi].x = particle[atomi].Fx;
    _grad_V[atomi].y = particle[atomi].Fy;
    _grad_V[atomi].z = particle[atomi].Fz;
  }
  
  V_init = V;
  
  FILE* energy_file = fopen ("E.dat" , "w");

  if(!energy_file){
    cerr<<"cannot open optimization energy file for writing\n";
    exit(0);
  }
  
  FILE* path_file = fopen("Traj-1.xyz" , "w");
  
  if(! path_file){
    cout<<"cannot open optimization path file for writing\n";
    exit(0);
  }
  
  fprintf(energy_file,"%i %.16f\n",0,V_init);
  fflush(energy_file);
  
  Print_Traj(path_file,Natoms,0,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
  frame_counter ++;
  
// Starting optimization *****************************************************************************

  //cerr<<"\n\nBeginning steepest descent algorithem\n\n";
  //cerr<<"0%                                            100%\n";
  
  for(iter=1 ; iter < ITMAXSD  ; iter++){
    
    if(Periodic_flag == 1)    
      relax_cell(Natoms,NBond,AtomParams,BondParams,Bonded,Bmat,qvec,Avec,rij,rijShiled_1,L,particle,Normal,Fc,logfile,angle,BondNum,Neighb,NNeighb,Coulomb_flag,Interlayer,Normal_atom,NList,List,R0);
    
    optimize_1D(temp_particle,_grad_V,min_V,_dV,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag,Xmin,Ymin,Zmin,NList,List,R0);
    
    if(!(iter%NConserve)){
      fprintf(energy_file , "%i %.16f\n" , iter , min_V);
      fflush(energy_file);
    }
    
    if(!(iter%NTraj)){
      if(!(frame_counter % Frame_limit)){
	fclose(path_file);
	file_counter ++;
	sprintf ( optimization_file_name , "Traj-%d.xyz" , file_counter );
	path_file = fopen ( optimization_file_name , "w" );
	
	if(!path_file){
	  cerr<<"cannot open optimization file for writing\n";
	  exit(0);
	}
      }

      Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      frame_counter ++;
    }
    
    if(( 2.0 * fabs( min_V - V ) ) <= (OptTol * ( fabs(min_V) + fabs(V) + EPS ))){
      fprintf(energy_file,"%i %.16f\n",iter,min_V);
      fflush(energy_file);
      
      Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      Print_Traj_gv(Natoms,particle);
      frame_counter ++;
      
      delete [] _grad_V;
      delete [] temp_particle;
      delete [] _dV;
      
      cerr<<"Completed optimization after: "<<iter<<" iterations.\n\n";
      cerr<<"The initial value of the potential was: "<<V_init<<" and the final is: "<<min_V<<" .\n";
      fclose(path_file);
      fclose(energy_file);
      return;
    }
    
    //calc_force(_grad_V, V, atoms, num_of_atoms, num_of_neighbors , list_of_neighbors, fix , num_of_edge_atoms);

    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
    
    if(Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    
    V = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);    
    
    Calc_Force(_grad_V,Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,Bmat,qvec,Avec,dqvec,rijShiled_1,logfile,angle,Normal,Normal_atom,Interlayer,Fc,Neighb,NNeighb,Coulomb_flag);    
    
#pragma omp parallel for 
    for(atomi=0 ; atomi < Natoms ; atomi++){
      _grad_V[atomi].x = particle[atomi].Fx;
      _grad_V[atomi].y = particle[atomi].Fy;
      _grad_V[atomi].z = particle[atomi].Fz;
    }
    
    //#pragma omp parallel for reduction(+:_grad_V_2) 
    for(atomi=0 ; atomi < Natoms ; atomi++) 
      //#pragma omp critical
      _grad_V_2 += (_grad_V[atomi].x * _grad_V[atomi].x + _grad_V[atomi].y * _grad_V[atomi].y + _grad_V[atomi].z * _grad_V[atomi].z);
    
    
    if(_grad_V_2 == 0.0){
      fprintf(energy_file,"%i %.16f\n",iter,min_V);
      fflush(energy_file);
      
      Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
      frame_counter ++;
      
      delete [] _grad_V;
      delete [] temp_particle;
      delete [] _dV;
      
      cerr<<"Completed optimization after: "<<iter<<" iterations.\n\n";
      cerr<<"The initial value of the potential was: "<<V_init<<" and the final is: "<<min_V<<" .\n";
      fclose(path_file);
      fclose(energy_file);
      return;
    }
    /* 
    if((double(iter)/ITMAXSD)*100 > precentage){
    cerr<<'*';
    precentage += 2;
    }
    */
  }
  
  fprintf(energy_file,"%i %.16f\n",iter,min_V);
  fflush(energy_file);

  Print_Traj(path_file,Natoms,iter,particle,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,L);
  frame_counter ++;
    
  cerr<<"Exceeded Max. number of iterations for optimization\n";
  fclose(path_file);
  fclose(energy_file);
  
  exit(0);
}


void scale_coords(int dir,int Natoms,point L,atom *particle, double factor){
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

void relax_cell(int Natoms,int *NBond,AtomParamsStruct *AtomParams,BondParamsStruct *BondParams,int *Bonded,double *Bmat,double *qvec,double *Avec,double *rij,double *rijShiled_1,point L, atom *particle,point *Normal,double *Fc,FILE* logfile,double *angle,int *BondNum,int *Neighb,int *NNeighb,int Coulomb_flag,int Interlayer,int *Normal_atom,int *NList,int *List,point *R0){
  double origEnergy;
  double newEnergy;
  double scaleup = 1.001;
  double scaledown = 1./scaleup;
  ///static const char message[]="\rrelax box size %c: %.3f    ";
  int nchange=0, maxchange=100,dir;
  //try to increase the volume
  //printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);
  
  calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
  if(Interlayer == 1)
    calc_Normal(Natoms,Normal_atom,particle,Normal,L);
  origEnergy=calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
  
  for(dir=0;dir<3;dir++){
    scale_coords(dir,Natoms,L,particle,scaleup);
    //newEnergy=E_interlayer()+E_intralayer();
    calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
    if(Interlayer == 1)
      calc_Normal(Natoms,Normal_atom,particle,Normal,L);
    newEnergy=calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
    if(newEnergy < origEnergy){ //keep increasing the volume
      origEnergy = newEnergy;
      while(nchange<maxchange){
	//printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);
	scale_coords(dir,Natoms,L,particle,scaleup);
	calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
	if(Interlayer == 1)
	  calc_Normal(Natoms,Normal_atom,particle,Normal,L);
	newEnergy = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
	if(newEnergy < origEnergy){
	  origEnergy = newEnergy;
	}
	else{ //step back and exit
	  //printf(message,dir,pow(scaleup,--nchange)); fflush(stdout);
	  scale_coords(dir,Natoms,L,particle,scaledown);
	  break;
	}
      }
    }
    else { //try to reduce the volume
      origEnergy = newEnergy;
      while(nchange>-maxchange){
	//printf(message,dir,pow(scaleup,--nchange)); fflush(stdout);
	scale_coords(dir,Natoms,L,particle,scaledown);
	calc_force_arrays(Natoms,NBond,Bonded,particle,BondParams,BondNum,rij,L,AtomParams,angle,Fc,Neighb,NNeighb,rijShiled_1,NList,List,R0);
	if(Interlayer == 1)
	  calc_Normal(Natoms,Normal_atom,particle,Normal,L);
	newEnergy = calc_Potential(Natoms,NBond,Bonded,AtomParams,BondParams,particle,Bmat,qvec,Avec,rijShiled_1,logfile,rij,angle,BondNum,Normal,Fc,Neighb,NNeighb,L,Coulomb_flag);
	if(newEnergy < origEnergy){
	  origEnergy = newEnergy;
	}
	else{ //step back and exit
	  //printf(message,dir,pow(scaleup,++nchange)); fflush(stdout);
	  
	  scale_coords(dir,Natoms,L,particle,scaleup);
	  break;
	}
      }
    }
  }
  //printf("\n"); fflush(stdout);
  return;
}
 
