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
  
  return(counter-1);
}

void ReadCoor(char *filename, int Natoms, atom* particle,int interlayer,point &L)
{
  int i,Nlayer01;
  FILE *file;

  
  if((file = fopen(filename,"r")) == NULL){
    cerr<<"\n ERROR: cannot open file: "<<filename<<" . Ending session!"<<endl;
    exit(0);
  }
  
  fscanf(file,"%lg %lg %lg\n",&L.x,&L.y,&L.z);
  
  for(i=0 ; i < Natoms ; i++){
    
    //fscanf(file,"%i %lg %lg %lg\n",&particle[i].type, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    fscanf(file,"%i %i %lg %lg %lg\n",&particle[i].type,&particle[i].layer, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    
    //if(i < Nlayer01)particle[i].layer=0;
    //else particle[i].layer=1;
    
    if(particle[i].type == 1){
      particle[i].Mass = 1.0080;
      particle[i].type =1;
    }
    else if(particle[i].type == 5){
      particle[i].type =0;
      particle[i].Mass    = 12.01;
    }
    else if(particle[i].type == 7){
      particle[i].type =0;
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
