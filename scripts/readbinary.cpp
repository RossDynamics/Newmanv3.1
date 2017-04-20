#include <iostream>
#include <fstream>
using namespace std;

int  main(int argc, char *argv[]) {
  
  int numdrifters;
  double frametime;
  double X[1000];
  
  ifstream fin(argv[1],ios::binary);
  
  for(int tt=0;tt<3;++tt) {
 
  
  fin.read((char *)&numdrifters, sizeof(int) );
  fin.read((char *)&frametime, sizeof(double) );
  
  for(int d=0;d<3;++d) {
    
    
    fin.read((char *)&X[0], numdrifters*sizeof(double));
    
    cout << "numdrifters = " << numdrifters << endl;
    cout << "frametime = " << frametime << endl;
    cout << "X[" << d << "] = " << endl;
    
    for(int ii=0;ii<numdrifters;++ii)
      cout << X[ii] << endl;
    
  }
  }
  fin.close();
  
  
}