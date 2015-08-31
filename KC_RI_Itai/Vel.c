#include "declarations.h"

//point vel_center(atom *particle, double MTot, int Natoms);

void vel(atom *particle, double T, int Natoms, long *idum)
{
  int i;
  double SigmaB, MTot;
  point Vcm; 

  MTot = 0.0;
  
  for(i=0 ; i < Natoms ; i++){
    SigmaB = sqrt(kB * T / particle[i].Mass);
    
    MTot += particle[i].Mass;
    
    particle[i].vx = SigmaB * randv(idum);
    particle[i].vy = SigmaB * randv(idum);
    particle[i].vz = SigmaB * randv(idum);
  }
  
  // subtract the  center of mass velocity
  
  Vcm = vel_center(particle,MTot,Natoms);
  
  for(i=0 ; i < Natoms ; i++){
    particle[i].vx -= Vcm.x;
    particle[i].vy -= Vcm.y;
    particle[i].vz -= Vcm.z;
  }
}

/****************************************************************************/

point vel_center(atom *particle, double MTot, int Natoms)
{
  int i;
  point Vcm;
  double MTot_1 = ( 1.0 / MTot );
  
  Vcm.x = Vcm.y = Vcm.z = 0.0;
    
  for(i=0 ; i < Natoms ; i++){
    Vcm.x += particle[i].Mass * particle[i].vx;
    Vcm.y += particle[i].Mass * particle[i].vy;
    Vcm.z += particle[i].Mass * particle[i].vz;
  }
  
  Vcm.x *= MTot_1;
  Vcm.y *= MTot_1;
  Vcm.z *= MTot_1;
  
  return(Vcm);
}
