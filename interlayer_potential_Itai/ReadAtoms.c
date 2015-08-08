#include "declarations.h"

int GetNatoms(char *filename,int Input_flag)
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
  
  if(Input_flag==0 || Input_flag==1 || Input_flag==2)
    return(counter-1);
  else if (Input_flag==3 || Input_flag==4)
    return(counter-2);
  else{
    cerr<<"wrong Input_flag value in function GetNatoms"<<endl;
    exit(0);
  }
}

void ReadCoor(char *filename, int Natoms, atom* particle,int interlayer,point &L,int Input_flag)
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
  
  if(Input_flag==0 || Input_flag==1 || Input_flag==2)//My_MD
    fscanf(file,"%lg %lg %lg\n",&L.x,&L.y,&L.z);
  else if(Input_flag==3){ //JMOL
    fscanf(file,"%i\n\n",&Tmp);
  }
  else if (Input_flag==4){
    fscanf(file,"%i %lg %lg %lg %i\n",&Tmp,&L.x,&L.y,&L.z,&Tmp2);//Robertos Input scheme
    fscanf(file,"%i %i\n",&Nlayer01,&Tmp3);
  }
  else{
    cerr<<"wrong Input_flag value in function ReadCoor"<<endl;
    exit(0);
  }
  
  for(i=0 ; i < Natoms ; i++){
    particle[i].layer=0;
    particle[i].lock =0;
    
    //fscanf(file,"%i %lg %lg %lg\n",&particle[i].type, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    if(Input_flag==0)
      fscanf(file,"%i %lg %lg %lg\n",&particle[i].type, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    else if(Input_flag==1) 
      fscanf(file,"%i %i %lg %lg %lg\n",&particle[i].type,&particle[i].layer, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    else if(Input_flag==2) 
      fscanf(file,"%i %i %i %lg %lg %lg\n",&particle[i].type,&particle[i].layer,&particle[i].lock, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
    else if(Input_flag==3){
      fscanf(file,"%c %lg %lg %lg \n",&atom_type, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2]);
      if(atom_type == 'N')particle[i].type=7;
      else if(atom_type == 'B')particle[i].type=5;
      else if(atom_type == 'C')particle[i].type=6;
      else{
	cerr<<"Currently works for B-N systems! Ending session.";
	exit(0);
      }
    }
    else if(Input_flag==4){
      fscanf(file,"%c %lg %lg %lg %lg %lg %lg\n",&atom_type, &particle[i].r[0], &particle[i].r[1], &particle[i].r[2],&Tmp,&Tmp,&Tmp);
      if(atom_type == 'N')particle[i].type=7;
      else if(atom_type == 'B')particle[i].type=5;
      else if(atom_type == 'C')particle[i].type=6;
      else{
	cerr<<"Currently works for B-N systems! Ending session.";
	exit(0);
      }
      if(i < Nlayer01)particle[i].layer=0;
      else particle[i].layer=1;
    }
    
    if(particle[i].type == 1){
      particle[i].Mass = 1.0080;
      particle[i].type =1;
    }
    else if(particle[i].type == 5){
      //particle[i].type =0;
      particle[i].Mass    = 12.01;
    }
    else if(particle[i].type == 7){
      particle[i].Mass    = 12.01;
    }
    else if(particle[i].type == 6){
      particle[i].type =0;
      particle[i].Mass    = 12.01;
    }
    else if(particle[i].type == 0){
      particle[i].Mass    = 12.01;
    }
    else{
      cerr<<"Currently works for B-N-H systems! Ending session.";
      exit(0);
    }
    
  }
  fclose(file);
}
