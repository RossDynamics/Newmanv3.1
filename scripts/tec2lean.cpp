/* Utility to convert Tecplot file to lean velocity file */
#include <math.h>
#include <stdio.h>
#include <iostream>
using namespace std;
#include <string.h>
#include <stdlib.h>


int main(int argc, char *argv[]) {
  
  if(argc != 3) {
    cout << "Usage: tec2lean [tecplotfile.dat] [ND]" << endl;
    exit(1);
  }
  
  
  int ND = atoi(argv[2]);  
  FILE *ftec, *flean;
  char leanfile[200];
  char tecfile[200];
  
  strcpy(tecfile,argv[1]);
  
  if((ftec = fopen( tecfile,"r"))==NULL) {
		cout << "Error opening "<< tecfile << endl;
		exit(1);
	}

  sprintf(leanfile,"lean_%s",tecfile);
  
  if((flean = fopen( leanfile,"w"))==NULL) {
		cout << "Error opening "<< leanfile << endl;
		exit(1);
	}
  
  char buf[400];
  
  // Copy Header
  fgets(buf,400,ftec);
  fputs(buf,flean);
  cout << buf << endl;
  
  for(int d=1;d<=ND*3;++d)
    fprintf(flean,"#%d\n",d);
  
  double u=0.0;
  double v=0.0;
  double w=0.0;
  
  while( fgets(buf,400,ftec) != NULL) {
  
    if(!( strncmp(buf,"Z",1) )) {
      fputs(buf,flean);
      cout << buf << endl;
    }
    else {
      if(ND==2) {
        sscanf(buf,"%*f %*f %lf %lf\n",&u,&v);
        fprintf(flean,"%f %f\n",u,v);
      }
      else {
        sscanf(buf,"%*f %*f %*f %lf %lf %lf\n",&u,&v,&w);
        fprintf(flean,"%f %f %f\n",u,v,w);
      }
    }
  }
  
  
  fclose(ftec);
  fclose(flean);
  
  return 0;
}