#include "declarations.h"

/***********************************************************/
/* Calculate and fill the following arrays:                */
/*                                                         */
/* rij[Natoms*Natoms] - interatom distance array.          */
/*                                                         */
/* NBond[Natoms]      - number of atoms bonded to each     */
/*                      atom - two atoms are consedered to */
/*			be bonded if BO_ij != 0 namely if  */
/*                      their distance is less than 3A     */
/*			(BOCutoff).                        */
/*                                                         */
/* Continue explanation ...                                */
/***********************************************************/


void calc_angles(double *angle, int Natoms, int *NBond, int *Bonded, double *rij, atom *particle,point L)
{
    // This routine calculates and fills the angles array.  For each
  // atom j all the neighboring atom pairs (i and k) are considered
  // and the angle ijk is calculated, where j is the central atom. The
  // calculation is done using the scalar product formula and the acos
  // function.
  //#pragma omp parallel
  //{
  int i, k, atomi, atomk, atomj, index;
  double DummyAngle,inside_acos,r_ij,r_jk;
  double rx,ry,rz;
  point rji, rjk;
  
  DummyAngle = -1.0;
  
#pragma omp parallel for shared(DummyAngle) private(i,k,atomi,atomk,index,inside_acos,rji,rjk,rx,ry,rz,r_jk,r_ij)
  for(atomj=0 ; atomj < Natoms ; atomj++){ // Go over all central atoms j.
    
    for(i=0 ; i < NBond[atomj] ; i++){ // Go over all j neighbors i.

      atomi = Bonded[atomj*MaxNBond + i];
      
      //r_ij = R_ij(particle,atomi,atomj,L);
      r_ij   = rij[atomj*MaxNBond + i];
      rji.x = particle[atomi].r[0] - particle[atomj].r[0];
      rji.y = particle[atomi].r[1] - particle[atomj].r[1];
      rji.z = particle[atomi].r[2] - particle[atomj].r[2];

      rji.x -= L.x * rint(rji.x / (L.x + TINY));
      rji.y -= L.y * rint(rji.y / (L.y + TINY));
      rji.z -= L.z * rint(rji.z / (L.z + TINY));
      
      for(k=0 ; k < i ; k++){ // Go over all j neighbors k != i.
	index = atomj * MaxNBond * MaxNBond + i * MaxNBond + k;
	atomk = Bonded[atomj*MaxNBond + k];
	
	//r_jk = R_ij(particle,atomk,atomj,L);
	r_jk   = rij[atomj*MaxNBond + k];
	rjk.x = particle[atomk].r[0] - particle[atomj].r[0];
	rjk.y = particle[atomk].r[1] - particle[atomj].r[1];
	rjk.z = particle[atomk].r[2] - particle[atomj].r[2];
	
	rjk.x -= L.x * rint(rjk.x / (L.x + TINY));
	rjk.y -= L.y * rint(rjk.y / (L.y + TINY));
	rjk.z -= L.z * rint(rjk.z / (L.z + TINY));
	
	inside_acos = ((rji.x*rjk.x + rji.y*rjk.y + rji.z*rjk.z) / (r_ij* r_jk));
	
	if (fabs(inside_acos) < 1.0){
	  angle[index] = acos((rji.x*rjk.x + rji.y*rjk.y + rji.z*rjk.z) / (r_ij* r_jk));
	}
	else if (inside_acos >= 1.0) angle[index] = 0.0;
	else angle[index] = PIE;
      }
      
      k = i;
      index = atomj * MaxNBond * MaxNBond + i * MaxNBond + k;
      angle[index] = DummyAngle;
      
      for(k=i+1 ; k < NBond[atomj] ; k++){ // Go over all j neighbors k != i.
	index = atomj * MaxNBond * MaxNBond + i * MaxNBond + k;
	atomk = Bonded[atomj*MaxNBond + k];
	
	//r_jk = R_ij(particle,atomk,atomj,L);
	r_jk   = rij[atomj*MaxNBond + k];
	rjk.x = particle[atomk].r[0] - particle[atomj].r[0];
	rjk.y = particle[atomk].r[1] - particle[atomj].r[1];
	rjk.z = particle[atomk].r[2] - particle[atomj].r[2];

	rjk.x -= L.x * rint(rjk.x / (L.x + TINY));
	rjk.y -= L.y * rint(rjk.y / (L.y + TINY));
	rjk.z -= L.z * rint(rjk.z / (L.z + TINY));
	
	inside_acos = ((rji.x*rjk.x + rji.y*rjk.y + rji.z*rjk.z) / (r_ij* r_jk));
	
	if (fabs(inside_acos) < 1.0){
	  angle[index] = acos((rji.x*rjk.x + rji.y*rjk.y + rji.z*rjk.z) / (r_ij* r_jk));
	}
	else if (inside_acos >= 1.0) angle[index] = 0.0;
	else angle[index] = PIE;
	
      }
      
      for(k=NBond[atomj] ; k < MaxNBond ; k++){
	index = atomj * MaxNBond * MaxNBond + i * MaxNBond + k;
	angle[index] = DummyAngle;
      }
    }
    
    for(i=NBond[atomj] ; i < MaxNBond ; i++){
      for(k=0 ; k < MaxNBond ; k++){
	index = atomj * MaxNBond * MaxNBond + i * MaxNBond + k;
	angle[index] = DummyAngle;
      }
    }
    // }
  }
}

double calc_val_angle(int atomi,int atomj, int atomk, int Natoms, int *NBond, int *Bonded,double r_jk, double r_ij, atom *particle,point L,dr_ij &rji,dr_ij &rjk,double inside_acos)
{
  // This routine calculates and fills the angles array.  For each
  // atom j all the neighboring atom pairs (i and k) are considered
  // and the angle ijk is calculated, where j is the central atom. The
  // calculation is done using the scalar product formula and the acos
  // function.
  //#pragma omp parallel
  //{

  //cerr<<"dir= " <<0<<"rji.r[dir]_01= "<<rji.r[0]<<endl;
  //cerr<<"dir= " <<1<<"rji.r[dir]_01= "<<rji.r[1]<<endl;
  //cerr<<"dir= " <<2<<"rji.r[dir]_01= "<<rji.r[2]<<endl;
  
  //double inside_acos;
  double rx,ry,rz,angle;
  
  angle=0;
  
  //point rji, rj;k
  
  rji.r[0] = particle[atomi].r[0] - particle[atomj].r[0];
  rji.r[1] = particle[atomi].r[1] - particle[atomj].r[1];
  rji.r[2] = particle[atomi].r[2] - particle[atomj].r[2];
  
  rji.r[0] -= L.x * rint(rji.r[0] / (L.x + TINY));
  rji.r[1] -= L.y * rint(rji.r[1] / (L.y + TINY));
  rji.r[2] -= L.z * rint(rji.r[2] / (L.z + TINY));
  
  //cerr<<"dir= " <<0<<"rji.r[dir]_01= "<<rji.r[0]<<endl;
  //cerr<<"dir= " <<1<<"rji.r[dir]_01= "<<rji.r[1]<<endl;
  //cerr<<"dir= " <<2<<"rji.r[dir]_01= "<<rji.r[2]<<endl;
  
  rjk.r[0] = particle[atomk].r[0] - particle[atomj].r[0];
  rjk.r[1] = particle[atomk].r[1] - particle[atomj].r[1];
  rjk.r[2] = particle[atomk].r[2] - particle[atomj].r[2];
  
  rjk.r[0] -= L.x * rint(rjk.r[0] / (L.x + TINY));
  rjk.r[1] -= L.y * rint(rjk.r[1] / (L.y + TINY));
  rjk.r[2] -= L.z * rint(rjk.r[2] / (L.z + TINY));
  
  inside_acos = ((rji.r[0]*rjk.r[0] + rji.r[1]*rjk.r[1] + rji.r[2]*rjk.r[2]) / (r_ij* r_jk));
  
  if (fabs(inside_acos) < 1.0){
    angle = acos((rji.r[0]*rjk.r[0] + rji.r[1]*rjk.r[1] + rji.r[2]*rjk.r[2]) / (r_ij* r_jk));
  }
  else if (inside_acos >= 1.0) angle = 0.0;
  else angle = PIE;

  return(angle);
}
/*********************************************************************************************************/

/*********************************************************************************************************/

double calc_dihedral(int atomi,int atomj, int atomk, int atoml, atom *particle,point L)
{
  // This routine calculates and fills the omega array.  For each
  // two neighboring atoms jk all the neighboring atoms (i and l where i is a neighbor of j and l is a neighbor of k)
  // are considered and the dihedral angle ijkl is calculated, where jk is the central atom pair. The
  // calculation is done using the atan funtion.
  
  double a,b,size_jk,Omega;
  
  point rji,rjk,rkl,Normal_ij_jk,Normal_jk_kl;
  
  rji.x = particle[atomj].r[0] - particle[atomi].r[0];
  rji.y = particle[atomj].r[1] - particle[atomi].r[1];
  rji.z = particle[atomj].r[2] - particle[atomi].r[2];
  
  rji.x -= L.x * rint(rji.x / (L.x + TINY));
  rji.y -= L.y * rint(rji.y / (L.y + TINY));
  rji.z -= L.z * rint(rji.z / (L.z + TINY));
  
  rjk.x = particle[atomk].r[0] - particle[atomj].r[0];
  rjk.y = particle[atomk].r[1] - particle[atomj].r[1];
  rjk.z = particle[atomk].r[2] - particle[atomj].r[2];
  
  rjk.x -= L.x * rint(rjk.x / (L.x + TINY));
  rjk.y -= L.y * rint(rjk.y / (L.y + TINY));
  rjk.z -= L.z * rint(rjk.z / (L.z + TINY));
  
  rkl.x = particle[atoml].r[0] - particle[atomk].r[0];
  rkl.y = particle[atoml].r[1] - particle[atomk].r[1];
  rkl.z = particle[atoml].r[2] - particle[atomk].r[2];
  
  rkl.x -= L.x * rint(rkl.x / (L.x + TINY));
  rkl.y -= L.y * rint(rkl.y / (L.y + TINY));
  rkl.z -= L.z * rint(rkl.z / (L.z + TINY));
  
  Normal_ij_jk.x =  rji.y * rjk.z - rji.z * rjk.y;
  Normal_ij_jk.y =-(rji.x * rjk.z - rji.z * rjk.x);
  Normal_ij_jk.z =  rji.x * rjk.y - rji.y * rjk.x;
  
  Normal_jk_kl.x =  rjk.y * rkl.z - rjk.z * rkl.y;
  Normal_jk_kl.y =-(rjk.x * rkl.z - rjk.z * rkl.x);
  Normal_jk_kl.z =  rjk.x * rkl.y - rjk.y * rkl.x;
  
  size_jk = sqrt(sqr(rjk.x) + sqr(rjk.y) + sqr(rjk.z));
  a       = size_jk*(rji.x*Normal_jk_kl.x + rji.y*Normal_jk_kl.y + rji.z*Normal_jk_kl.z);
  b       = Normal_ij_jk.x*Normal_jk_kl.x + Normal_ij_jk.y*Normal_jk_kl.y + Normal_ij_jk.z*Normal_jk_kl.z;
  Omega = atan2(a,b);
  return(Omega);
}
void calc_Normal(int Natoms,int *Normal_atom,atom *particle,Normal_struct *Normal,point L)
{
  //#pragma omp parallel
  //{
  int atomi,atomj,atomk,atoml;
  double vector_1x,vector_1y,vector_1z;
  double vector_2x,vector_2y,vector_2z;
  double Normal_length;
  
  //#pragma omp for
#pragma omp parallel for private(atomi,atomk,atomj,vector_1x,vector_1y,vector_1z,vector_2x,vector_2y,vector_2z,Normal_length)
  for(atomi=0; atomi < Natoms; atomi++){
    Normal[atomi].x = 0.0;
    Normal[atomi].y = 0.0;
    Normal[atomi].z = 0.0;
    if (particle[atomi].type ==1){
      Normal[atomi].x = 0.0;
      Normal[atomi].y = 0.0;
      Normal[atomi].z = 0.0;
    }
    else{
      atomj=Normal_atom[atomi*3 + 0];
      atomk=Normal_atom[atomi*3 + 1];
      atoml=Normal_atom[atomi*3 + 2];
      
      Normal[atomi].vector_1x = particle[atomj].r[0] - particle[atomk].r[0];
      Normal[atomi].vector_1y = particle[atomj].r[1] - particle[atomk].r[1];
      Normal[atomi].vector_1z = particle[atomj].r[2] - particle[atomk].r[2];
      
      Normal[atomi].vector_1x -= L.x * rint(Normal[atomi].vector_1x / (L.x + TINY));
      Normal[atomi].vector_1y -= L.y * rint(Normal[atomi].vector_1y / (L.y + TINY));
      Normal[atomi].vector_1z -= L.z * rint(Normal[atomi].vector_1z / (L.z + TINY));
      
      Normal[atomi].vector_2x = particle[atoml].r[0] - particle[atomk].r[0];
      Normal[atomi].vector_2y = particle[atoml].r[1] - particle[atomk].r[1];
      Normal[atomi].vector_2z = particle[atoml].r[2] - particle[atomk].r[2];
      
      Normal[atomi].vector_2x -= L.x * rint(Normal[atomi].vector_2x / (L.x + TINY));
      Normal[atomi].vector_2y -= L.y * rint(Normal[atomi].vector_2y / (L.y + TINY));
      Normal[atomi].vector_2z -= L.z * rint(Normal[atomi].vector_2z / (L.z + TINY));
      
      Normal[atomi].not_norm_x = Normal[atomi].vector_1y*Normal[atomi].vector_2z - Normal[atomi].vector_2y*Normal[atomi].vector_1z;
      Normal[atomi].not_norm_y = Normal[atomi].vector_1z*Normal[atomi].vector_2x - Normal[atomi].vector_1x*Normal[atomi].vector_2z;
      Normal[atomi].not_norm_z = Normal[atomi].vector_1x*Normal[atomi].vector_2y - Normal[atomi].vector_1y*Normal[atomi].vector_2x;
      
      Normal[atomi].Normal_length   = sqrt(sqr(Normal[atomi].not_norm_x)+sqr(Normal[atomi].not_norm_y)+sqr(Normal[atomi].not_norm_z));
      
      Normal[atomi].x = Normal[atomi].not_norm_x/(Normal[atomi].Normal_length);
      Normal[atomi].y = Normal[atomi].not_norm_y/(Normal[atomi].Normal_length);
      Normal[atomi].z = Normal[atomi].not_norm_z/(Normal[atomi].Normal_length);
      /*
      if(atomi==640){
	cerr<<"atomj= "<<atomj<<endl;
	cerr<<"atoml= "<<atoml<<endl;
	cerr<<"atomk= "<<atomk<<endl;
	
	cerr<<"Normal[atomi].Normal_length="<<Normal[atomi].Normal_length<<endl;
	cerr<<"Normal[atomi].not_norm_x="<<Normal[atomi].not_norm_x<<endl;
	cerr<<"Normal[atomi].not_norm_y="<<Normal[atomi].not_norm_y<<endl;
	cerr<<"Normal[atomi].not_norm_z="<<Normal[atomi].not_norm_z<<endl;
	cerr<<"Normal[atomi].vector_1y="<<Normal[atomi].vector_1y<<endl;
	cerr<<"Normal[atomi].vector_2x="<<Normal[atomi].vector_2x<<endl;
	cerr<<"=Normal[atomi].vector_2y"<<Normal[atomi].vector_2y<<endl;
	cerr<<"=Normal[atomi].vector_1x"<<Normal[atomi].vector_1x<<endl;
	cerr<<"Normal[atomi].x="<<Normal[atomi].x<<endl;
	cerr<<"Normal[atomi].y="<<Normal[atomi].y<<endl;
	cerr<<"Normal[atomi].z="<<Normal[atomi].z<<endl;
      }
      */
      //Normal[atomi].x = 0.0;
      //Normal[atomi].y = 0.0;
      //Normal[atomi].z = 1.0;
    }
  }
  //}
}
void calc_dNormal_k(int &atomi,int &atomn,int &dir,dr_ij *dN,atom *particle,Normal_struct *Normal,int *Normal_atom,point &L)
{
  double vector_1x,vector_1y,vector_1z;
  double vector_2x,vector_2y,vector_2z;
  int atomj,atomk,atoml;
  double Normal_length,d_Norm_length;
  double Normalx,Normaly,Normalz;
  
  atomj=Normal_atom[atomi*3 + 0];
  atomk=Normal_atom[atomi*3 + 1];
  atoml=Normal_atom[atomi*3 + 2];
  
  if(particle[atomi].type==1){
    dN[0].r[0] = 0.0;
    dN[0].r[1] = 0.0;
    dN[0].r[2] = 0.0;
  }
  else {
    if (atomn == atomj){
      if (dir == 0){
	dN[0].r[0] = 0.0;
	dN[0].r[1] = -Normal[atomi].vector_2z;
	dN[0].r[2] = Normal[atomi].vector_2y;
      }
      else if(dir == 1){
	dN[0].r[0] = Normal[atomi].vector_2z;
	dN[0].r[1] = 0.0;
	dN[0].r[2] =-Normal[atomi].vector_2x;
      }
      else if(dir == 2){
	dN[0].r[0] =-Normal[atomi].vector_2y;
	dN[0].r[1] = Normal[atomi].vector_2x;
	dN[0].r[2] = 0.0;
      }
    }
    else if (atomn == atomk){
      if (dir == 0){
	dN[0].r[0] = 0.0;
	dN[0].r[1] =-Normal[atomi].vector_1z + Normal[atomi].vector_2z;
	dN[0].r[2] =-Normal[atomi].vector_2y + Normal[atomi].vector_1y;
      }
      else if(dir == 1){
	dN[0].r[0] = -Normal[atomi].vector_2z + Normal[atomi].vector_1z;
	dN[0].r[1] = 0.0;
	dN[0].r[2] = -Normal[atomi].vector_1x + Normal[atomi].vector_2x;
      }
      else if(dir == 2){
	dN[0].r[0] =-Normal[atomi].vector_1y + Normal[atomi].vector_2y;
	dN[0].r[1] =-Normal[atomi].vector_2x + Normal[atomi].vector_1x;
	dN[0].r[2] = 0.0;
      }
    }
    else if (atomn == atoml){
      if (dir == 0){
	dN[0].r[0] = 0.0;
	dN[0].r[1] =Normal[atomi].vector_1z;
	dN[0].r[2] =-Normal[atomi].vector_1y;
      }
      else if(dir == 1){
	dN[0].r[0] =-Normal[atomi].vector_1z;
	dN[0].r[1] = 0.0;
	dN[0].r[2] =Normal[atomi].vector_1x;
      }
      else if(dir == 2){
	dN[0].r[0] =Normal[atomi].vector_1y;
	dN[0].r[1] =-Normal[atomi].vector_1x;
	dN[0].r[2] = 0.0;
      }
    }
    else{
      dN[0].r[0] =0.0;
      dN[0].r[1] =0.0;
      dN[0].r[2] =0.0;
    }
    d_Norm_length = pow(Normal[atomi].Normal_length,-3.0)*(Normal[atomi].not_norm_x*dN[0].r[0] + Normal[atomi].not_norm_y*dN[0].r[1] + Normal[atomi].not_norm_z*dN[0].r[2]);
    dN[0].r[0]= dN[0].r[0]/Normal[atomi].Normal_length-Normal[atomi].not_norm_x*d_Norm_length;
    dN[0].r[1]= dN[0].r[1]/Normal[atomi].Normal_length-Normal[atomi].not_norm_y*d_Norm_length;
    dN[0].r[2]= dN[0].r[2]/Normal[atomi].Normal_length-Normal[atomi].not_norm_z*d_Norm_length;
  }
  
}
void set_Normal_atom(double *rij,int Natoms,Normal_struct *Normal,atom *particle,int *Normal_atom,point L)
{
 //This function fills the Normal_atom array of the 3 nearest neighbors surrounding an atom which we use to calculate its normal. 
  double rx,ry,rz;
  int n,atomj,atomi;
  double r_ij,R;

  n=0;
  R=Normal_cuttoff;
  
  for(atomi=0; atomi < Natoms; atomi++){
    //R=Normal_cuttoff;
    //n=0;
    for(atomj = 0;atomj < Natoms;atomj++){
      
      if (atomi == atomj || (particle[atomi].layer!=particle[atomj].layer))continue;
      
      //r_ij=R_ij(particle,atomi,atomj, L);
      rx = particle[atomi].r[0] - particle[atomj].r[0];
      ry = particle[atomi].r[1] - particle[atomj].r[1];
      rz = particle[atomi].r[2] - particle[atomj].r[2];
      
      rx -= L.x * rint(rx / (L.x + TINY));
      ry -= L.y * rint(ry / (L.y + TINY));
      rz -= L.z * rint(rz / (L.z + TINY));
      
      r_ij = sqr(rx) + sqr(ry) + sqr(rz);
      
      if(r_ij > R*R)continue;
      else
	{
	  Normal_atom[atomi*3+n] = atomj;
	  n++;
	}
      
      if(n==3)break;
      //else if(atomj==(Natoms-1)){
      //atomj=-1;
      //n=0;
      //R = R + 0.1;
      //}
    }
    
    if(n!=3){
      atomi-=1;
      R = R + 0.2;
    }
    else R=Normal_cuttoff;
    
    n=0;
    //cerr<<"atomi="<<atomi<<endl;
    //cerr<<"atomi type="<<particle[atomi].type<<endl;
    //cerr<<"n="<<n<<endl;
    
    //if (n != 3 && particle[atomi].type != 1){
    //cerr<<"Error in set_Normal_atom function with atom="<<atomi<<": each atom should have three nearest neighbors excluding hydorgens"<<" atom type="<<particle[atomi].type<<" lock="<<particle[atomi].lock<<" num of neighbors="<<n<<endl;
    //exit(0);
    //}
    
  }
}

double R_ij(atom *particle,int atomi,int atomj, point L){
  
  double rx,ry,rz;
  
  rx = particle[atomi].r[0] - particle[atomj].r[0];
  ry = particle[atomi].r[1] - particle[atomj].r[1];
  rz = particle[atomi].r[2] - particle[atomj].r[2];
  
  rx -= L.x * rint(rx / (L.x + TINY));
  ry -= L.y * rint(ry / (L.y + TINY));
  rz -= L.z * rint(rz / (L.z + TINY));
  
  return(sqrt(sqr(rx) + sqr(ry) + sqr(rz)));
}
double R_PBC(atom *particle,int atomi,int atomj,int dir,double L){
  double R;
  R = particle[atomi].r[dir] - particle[atomj].r[dir];
  
  R -= L * rint(R / (L + TINY));
  
  return(R);
}
double Fc_(double r_ij,double Tersoff_R,double Tersoff_D)
{
  double Fc;
  
  if(r_ij < Tersoff_R+Tersoff_D){
    if(r_ij < Tersoff_R-Tersoff_D){
      Fc = 1;
    }
    else{
      Fc = 1/2-1/2*sin(PIE*0.5*(r_ij - Tersoff_R)/Tersoff_D);
    }
  }
  else Fc=0;
  
  return(Fc);
}

void UpdateNeighbList(int Dim, atom *particle, int *List, int *NList, point *R0, point L)
{
  int i, j;
  double R, rx, ry, rz;
  
#pragma omp parallel for  
  for(i=0 ; i < Dim ;i++) NList[i] = 0;
#pragma omp parallel for 
  for(i=0 ; i < Dim*MaxList ;i++) List[i] = -1;
  
  #pragma omp parallel for private(j,i,rx,ry,rz,R)
  for(i=0 ; i < Dim ;i++){
    for(j=0 ; j < i ; j++){
      rx = particle[i].r[0] - particle[j].r[0];
      ry = particle[i].r[1] - particle[j].r[1];
      rz = particle[i].r[2] - particle[j].r[2];
      
      rx -= L.x * rint(rx / (L.x + TINY));
      ry -= L.y * rint(ry / (L.y + TINY));
      rz -= L.z * rint(rz / (L.z + TINY));
      
      R  = sqrt(sqr(rx) + sqr(ry) + sqr(rz));
      
      if(R < Rcut_List){
   	List[i*MaxList + NList[i]] = j;
	//List[j*MaxList + NList[j]] = i;
	//#pragma omp critical
	NList[i]++;
	//NList[j]++;
      }
    }
    
    R0[i].x = particle[i].r[0];
    R0[i].y = particle[i].r[1];
    R0[i].z = particle[i].r[2];
    
  }
  
#pragma omp parallel for  
  for(i=0 ; i < Dim ;i++){
    if(NList[i] > MaxList){
      cerr<<"Number of neighbors of atom in List "<<i<<" exceeded maximum number of allowed neighbors in list! ("<<MaxList<<") ending session.\n";
      exit(0);
    }
  }
}

void Check_List(atom *particle, point *R0,int Natoms, int *List, int *NList, point L){
  int UpdatedList, i;
  double dx1, dy1, dz1, r1, R,R1, R2;
  R1=0;
  R2=0;
  
#pragma omp parallel for private(R,dx1,dy1,dz1)
  for(i=0; i<Natoms; i++){
    //for(UpdatedList=0, i=0; ((i<Natoms) && !UpdatedList) ; i++){
    dx1 = particle[i].r[0] - R0[i].x;
    dy1 = particle[i].r[1] - R0[i].y;
    dz1 = particle[i].r[2] - R0[i].z;
    
    dx1 -= L.x * rint(dx1/(L.x+TINY));
    dy1 -= L.y * rint(dy1/(L.y+TINY));
    dz1 -= L.z * rint(dz1/(L.z+TINY));
    
    R = sqrt(sqr(dx1) + sqr(dy1) + sqr(dz1));
    
    if(R > R1){
      R2=R1;
      R1=R;
    }
    
    if ((R1 + R2) > (Rcut_List-non_bond_cut)){
      UpdatedList = 1;
    }
  }
  
  if ((R1 + R2) > (Rcut_List-non_bond_cut)){
    UpdateNeighbList(Natoms,particle,List,NList,R0,L);
  }
  
}


