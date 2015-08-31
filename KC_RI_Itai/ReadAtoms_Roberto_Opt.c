#include "declarations.h"

int GetNatoms(char *filename)
{
  int counter=0;
  char tmp;
  FILE *file;
  
  if((file = fopen(filename,"r")) == NULL){
    cerr<<"\n ERROR: cannot open file: "<<filename<<" . Ending session!"<<endl;
    exit(0);
  }
  
  while(!feof(file)){
    fscanf(file, "%c", &tmp);
    if((tmp == '\n') && !feof(file)) counter++;
    if(feof(file) && (tmp != '\n'))  counter++;
  }
  
  fclose(file);
  
  return(counter-2);
}

void ReadCoor(char *filename, int Natoms, atom* particle,int interlayer,point &L)
{
  int i,Nlayer01;
  FILE *file;
  int Tmp,Tmp2,Tmp3;
  char atom_type;
  double Tmp4,Tmp5,Tmp6;
    
  if((file = fopen(filename,"r")) == NULL){
    cerr<<"\n ERROR: cannot open file: "<<filename<<" . Ending session!"<<endl;
    exit(0);
  }
  
  fscanf(file,"%i %lg %lg %lg %i\n",&Tmp,&L.x,&L.y,&L.z,&Tmp2);
  fscanf(file,"%i %i\n",&Nlayer01,&Tmp3);
  
  for(i=0 ; i < Natoms ; i++){
    
    //fscanf(file,"%i %lg %lg %lg\n",&particle[i].type, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    //fscanf(file,"%i %i %lg %lg %lg\n",&particle[i].type,&particle[i].layer, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    
    fscanf(file,"%c %lg %lg %lg %lg %lg %lg\n",&atom_type, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2],&Tmp4,&Tmp5,&Tmp6);
    
    if(atom_type == 'N')particle[i].type=7;
    else if(atom_type == 'B')particle[i].type=5;
    else if(atom_type == 'C')particle[i].type=6;
    else{
      cerr<<"Currently works for B-Nsystems! Ending session.";
      exit(0);
    }
    
    if(i < Nlayer01)particle[i].layer=0;
    else particle[i].layer=1;
    
    if(particle[i].type == 1){
      particle[i].Mass = 1.0080;
      particle[i].type =1;
    }
    else if(particle[i].type == 5){
      //particle[i].type =0;
      particle[i].Mass    = 12.01;
    }
    else if(particle[i].type == 7){
      //particle[i].type =0;
      particle[i].Mass    = 12.01;
    }
    else if(particle[i].type == 6){
      particle[i].type =0;
      particle[i].Mass    = 12.01;
    }
     else if(particle[i].type == 0){
      particle[i].type =0;
      particle[i].Mass    = 12.01;
    }
    else{
      cerr<<"Currently works for B-N-H systems! Ending session.";
      exit(0);
    }
  }
  fclose(file);
}
