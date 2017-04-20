// raw2tec.h
//
// Header file for the raw2tec module
//
//   Author: Philip Du Toit
//   Rev 1.0
//

#ifndef _RAW2TEC_H_
#define _RAW2TEC_H_

//#include <cstring>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include "tclap/CmdLine.h"
#include "TECXXX.h"

using namespace std;

ifstream fin;

int   SlideNumber;



double FrameTime;
double NumDrifters;
char *Path_Output = ".";
int Debug     = 0;
int VIsDouble = 0;
int DIsDouble = 1;
int JMax      = 1;
int KMax      = 1;
int ZoneType  = 0;      /* Ordered */
int StrandID  = 0;     /* StaticZone */
int ParentZn  = 0;      /* No Parent */
int IsBlock   = 1;      /* Block */
int ICellMax  = 0;
int JCellMax  = 0;
int KCellMax  = 0;
int NFConns   = 0;
int FNMode    = 0;
int ShrConn   = 0;
int sharing[7]={0,0,0,0,0,0,0};
int Time_Origin[6] = {0,0,0,0,0,0}; 
int Output_TRes = 1;
int Track_Storm;
int Track_Flag;
int  Info_Flag = 0;
int Date_Flag = 0;
int Verbose_Flag = 0;
int Track_Method;
int Atmos_Set = 0;
bool Atmos_Spherical = 0;
double Atmos_AltScale = 1.0;
double Atmos_Radius = 6378.1;
int Compute_Type = 0;
bool swap_endian = 0;
char ftlefilename[200];
string rawstring;
char rawfilename[200];
char outfilename[200];
char syscommand[200];
  double ftlemax[ND];
  double ftlemin[ND];
  int FTLE_Res[ND];
  double FTLE_Delta[ND];
  double Plot_Delta[ND];
  double Plot_Max[ND];
  double Plot_Min[ND];
  int Plot_Res[ND];
  int   LCS_NumFields;
  double *U;
vector<double> X;
vector<double> color;

void ConvertTimeToDate(double frametime,char *zonetitle);
bool FileExists(string strFilename) ;
void Format_FTLEData();
void Format_TraceData();
void Format_PlotData();
void ByteSwap(unsigned char * b, int n);



#endif


