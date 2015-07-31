#include "declarations.h"

void go_to_eol( FILE * file )
{
  char tmp='1';
  while( tmp != '\n' ) fscanf( file , "%c", &tmp );
}

void read_RUN_parameters(double &tini, double &tend, double &dt, double &T, int &NConserve, int &NTraj, double &OptTol, int &CalcType,int &Interlayer, int &Periodic_flag,int &nthreads,int &Coulomb_flag,int &Input_flag)
{
  double tConserve, tTraj;
  FILE* ParFile = fopen("run.par","r");
  
  fscanf(ParFile, "%lf", &tini);
  go_to_eol(ParFile);
  fscanf(ParFile, "%lf", &tend);
  go_to_eol(ParFile);
  fscanf(ParFile, "%lf", &dt);
  go_to_eol(ParFile);
  fscanf(ParFile, "%lf", &T);
  go_to_eol(ParFile);
  fscanf(ParFile, "%lf", &tConserve);
  go_to_eol(ParFile);
  fscanf(ParFile, "%lf", &tTraj);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &NConserve);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &NTraj);
  go_to_eol(ParFile);
  fscanf(ParFile, "%lf", &OptTol);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &CalcType);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &Interlayer);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &Periodic_flag);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &nthreads);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &Coulomb_flag);
  go_to_eol(ParFile);
  fscanf(ParFile, "%i",  &Input_flag);
  fclose(ParFile);
  
  
  //fclose(ParFile);
  if((CalcType != 0) && (CalcType != 1) && (CalcType != 2) && (CalcType != 3) && (CalcType != 4)){
    cerr<<"Unrecognized requested calculation type! Ending session.\n";
    exit(0);
  }
  if((Interlayer != 0) && (Interlayer != 1) && (Interlayer != 2)){
    cerr<<"Unrecognized requested Interlayer_flag! Ending session.\n";
    exit(0);
  }
  if((Periodic_flag != 0) && (Periodic_flag != 1)){
    cerr<<"Unrecognized requested Periodic_flag! Ending session.\n";
    exit(0);
  }
  if((Coulomb_flag != 0) && (Coulomb_flag != 1) && (Coulomb_flag != 2) && (Coulomb_flag != 3) && (Coulomb_flag != 4) && (Coulomb_flag != 5)){
    cerr<<"Unrecognized requested Coulomb_flag! Ending session.\n";
    exit(0);
  }
  if((Input_flag != 0) && (Input_flag != 1) && (Input_flag != 2) && (Input_flag != 3) && (Input_flag != 4)){
    cerr<<"Unrecognized requested Coulomb_flag! Ending session.\n";
    exit(0);
  }
  
  // Moving to atomic units
  
  // tini      *= (1e-15 * SEC2AKMA);
  //tend      *= (1e-15 * SEC2AKMA);
  //dt        *= (1e-15 * SEC2AKMA);
  //tConserve *= (1e-15 * SEC2AKMA);
  //tTraj     *= (1e-15 * SEC2AKMA);

  if(CalcType == 0){
    NConserve = int(tConserve / dt);
    NTraj     = int(tTraj     / dt);
  }
}
