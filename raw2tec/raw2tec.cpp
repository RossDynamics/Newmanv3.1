#define pout(a,b) cout << a << "\t\t\t " << b << endl;
#define poutd(a,b) cout << a << "[" << d << "]" << "\t\t\t " << b << endl;
#define BS(x) ByteSwap((unsigned char *) &x,sizeof(x))
#define LONGSTRING 200
#include "raw2tec.h"

using namespace TCLAP;

int main(int argc, char *argv[]) {
  
	// Define the command line object.
	CmdLine cmd("Convert newman binary output to Tecplot format.", ' ', "1.1");
  
  	// Define a switch and add it to the command line.
	SwitchArg infoSwitch("i","info","Display the data header info of a raw file.", false);
	cmd.add( infoSwitch );

  
  	// Define a switch and add it to the command line.
	SwitchArg endianSwitch("e","endian","Switch endian type of binary input.", false);
	cmd.add( endianSwitch );
  
    	// Define a switch and add it to the command line.
	SwitchArg verboseSwitch("v","verbose","Provide verbose output.", false);
	cmd.add( verboseSwitch );
  
  // Define a switch and add it to the command line.
	SwitchArg sphericalSwitch("s","spherical","Convert spherical coordinates to cartesian for plotting on the sphere.", false);
	cmd.add( sphericalSwitch );
  
  // Define a switch and add it to the command line.
	SwitchArg trackSwitch("t","track","Track the storm center.", false);
	cmd.add( trackSwitch );
  
   // Define a switch and add it to the command line.
	SwitchArg dateSwitch("d","date","Add date and time to the Zone header", false);
	cmd.add( dateSwitch );

  
	// Define a value argument and add it to the command line.
	ValueArg<double> altscaleValue("a","altscale","Factor by which to scale the atmosphere relative to the sphere radius.",false,1.0,"double");
	cmd.add( altscaleValue );

  UnlabeledValueArg<string> rawfile("rawfile", "Name of raw binary input file. For FTLE data, provide only the basefile name common to all files.",true,"rawfile.raw","filename",false,NULL);
	cmd.add( rawfile );
  
	// Parse the args.
	cmd.parse( argc, argv );
  
	// Get the value parsed by each arg. 
  Info_Flag = infoSwitch.getValue();
	Atmos_AltScale = altscaleValue.getValue();
  Atmos_Spherical = sphericalSwitch.getValue();
  Track_Flag = trackSwitch.getValue();
  Date_Flag = dateSwitch.getValue();
  Verbose_Flag = verboseSwitch.getValue();
	swap_endian = endianSwitch.getValue();
  
  rawstring = rawfile.getValue();
  strcpy(rawfilename,rawstring.c_str());

  if( FileExists(rawfilename) )  {
    fin.open(rawfilename,ios::binary);
    fin.read((char *)&Compute_Type, sizeof(int) );
    if(swap_endian) BS(Compute_Type);
    fin.close();
  }
  else {
    sprintf(ftlefilename,"%s0000.raw",rawfilename);
    if(FileExists(ftlefilename)) {
      Compute_Type=0;
    }
    else {
      cout << "ERROR: Could not open " << ftlefilename << endl;
      exit(1);
    }
  }
  
  switch(Compute_Type) {
    case 0:
      cout << "FTLE data ..." << endl;
      Format_FTLEData();
      break;
    case 1:
      cout << "Trace data ..." << endl;
      Format_TraceData();
      break;
    case 2:
      cout << "Plot data ..." << endl;
      Format_PlotData();
      break;
    default :
      cout << "ERROR: Invalid Compute type in " << rawfilename << endl;
      exit(1);
  }
  
  return 0;
}

void Format_FTLEData() {
  
  sprintf(outfilename,"%s.plt",rawfilename);
  
  int KMax=1;
  int ss=0;
  
  while(1){
    
    sprintf(ftlefilename,"%s%04d.raw",rawfilename,ss);
    
    if(!FileExists(ftlefilename)) 
      break;
    
    fin.open(ftlefilename,ios::binary);
    
    fin.read((char *)&Compute_Type, sizeof(int) );
    fin.read((char *)&Time_Origin, 6*sizeof(int) );
    fin.read((char *)&FrameTime, sizeof(double) );
    fin.read((char *)&Output_TRes, sizeof(int) );
    fin.read((char *)&Atmos_Set, sizeof(int) );
    fin.read((char *)&Atmos_Radius, sizeof(double) );
    fin.read((char *)&SlideNumber, sizeof(int) );
    fin.read((char *)&Track_Storm, sizeof(int) );
    fin.read((char *)&ftlemin[0], ND * sizeof(double) );
    fin.read((char *)&ftlemax[0], ND * sizeof(double) );
    fin.read((char *)&FTLE_Res[0], ND * sizeof(int) );
    fin.read((char *)&LCS_NumFields, sizeof(int) );
    
    if(swap_endian) {
      BS(Compute_Type);
      for(int d=0;d<6;++d)
        BS(Time_Origin[d]);
      BS(FrameTime);
      BS(Output_TRes);
      BS(Atmos_Set);
      BS(Atmos_Radius);
      BS(SlideNumber);
      BS(Track_Storm);
      for(int d=0;d<ND;++d) {
        BS(ftlemin[d]);
        BS(ftlemax[d]);
        BS(FTLE_Res[d]);
      }
      BS(LCS_NumFields);
    }
    
    if(Verbose_Flag || Info_Flag) {
      pout("Compute_Type",Compute_Type);
      cout << "Time_Origin  \t" << flush;
      for(int d=0;d<6;++d)
        cout << Time_Origin[d] << "\t " << flush;
      cout << endl;
      pout("FrameTime",FrameTime);
      pout("Output_TRes",Output_TRes);
      pout("Atmos_Set ",Atmos_Set );
      pout("Atmos_Radius",Atmos_Radius);
      pout("SlideNumber",SlideNumber);
      pout("Track_Storm ",Track_Storm );
      for(int d=0;d<ND;++d) {
        poutd("ftlemin",ftlemin[d]);
        poutd("ftlemax",ftlemax[d]);
        poutd("FTLE_Res",FTLE_Res[d]);
      }
      pout("LCS_NumFields",LCS_NumFields);
    }
    
    if(Info_Flag)
      exit(1);
    
    if(Track_Storm)
      Track_Flag ? Track_Method = 2 : Track_Method = 1;
    else
      Track_Method = 0;

    int FTLE_BlockSize = 1;
    for(int d=0;d<ND;++d)
      FTLE_BlockSize*=FTLE_Res[d];
    for(int d=0;d<ND;++d) {
      if(FTLE_Res[d] == 1 )
        FTLE_Delta[d] = ftlemax[d]-ftlemin[d];
      else
        FTLE_Delta[d] = (ftlemax[d]-ftlemin[d])/(FTLE_Res[d]-1); 
    }
    int IMax=FTLE_Res[0];
    int JMax=FTLE_Res[1];
    
    if(ss==0) {
      char varstring[100];
      char fieldstring[100];
      if(LCS_NumFields==1)
        strcpy(fieldstring,"FTLE");
      else
        strcpy(fieldstring,"FTLE Omega Eval smFTLE");
#if ND == 6
    KMax      = FTLE_Res[2];
    sprintf(varstring,"x y z xdot ydot zdot %s",fieldstring);
    TECINI110("FTLE Data", varstring, outfilename, Path_Output, &Debug, &VIsDouble);
      
#elif ND == 4 
      KMax      = FTLE_Res[2];
      sprintf(varstring,"x xdot H y %s",fieldstring);
      TECINI110("FTLE Data", varstring, outfilename, Path_Output, &Debug, &VIsDouble);
#elif ND == 3 
      KMax      = FTLE_Res[2];
      if(Atmos_Set) {
        if(Atmos_Spherical) {
          sprintf(varstring,"X Y Z %s",fieldstring);
          TECINI110("FTLE Data", varstring, outfilename, Path_Output, &Debug, &VIsDouble);
        }
        else {
          sprintf(varstring,"Longitude Latitude Altitude %s",fieldstring);
          TECINI110("FTLE Data",varstring, outfilename, Path_Output, &Debug, &VIsDouble);
        }
      }
      else {
        sprintf(varstring,"X Y Z %s",fieldstring);
        TECINI110("FTLE Data", varstring, outfilename, Path_Output, &Debug, &VIsDouble);
      }
#else
      if(Atmos_Set) {
        if(Atmos_Spherical) {
          sprintf(varstring,"X Y Z %s",fieldstring);
          TECINI110("FTLE Data", varstring, outfilename, Path_Output, &Debug, &VIsDouble);
        }
        else {
          sprintf(varstring,"Longitude Latitude %s",fieldstring);
          TECINI110("FTLE Data", varstring, outfilename, Path_Output, &Debug, &VIsDouble);
        }
      }
      else {
        sprintf(varstring,"X Y %s",fieldstring);
        TECINI110("FTLE Data", varstring, outfilename, Path_Output, &Debug, &VIsDouble);
      }
#endif 
    }
    
    if((ss>0) && (Track_Method!=1)) {
      for(int d=0;d<ND;++d)
        sharing[d]=1;
      if(Atmos_Set && Atmos_Spherical)
        sharing[2]=1;
    }
    
    char zonetitle[LONGSTRING];
    if(Date_Flag) {
      ConvertTimeToDate(FrameTime,zonetitle);
    }
    else {
      sprintf(zonetitle,"T = %f",FrameTime);
    }
    cout << zonetitle << endl;
    TECZNE110(zonetitle, &ZoneType, &IMax, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax, &FrameTime,
              &StrandID, &ParentZn, &IsBlock, &NFConns, &FNMode, NULL, NULL, sharing, &ShrConn);
    
    
   /* Track_method      Sharing     ftlemin
          0                 1         1
          1                 0         1
          2                 1         0
*/    
    if((ss==0) || (Track_Method==1)) {
      
      double *XYZvals = new (nothrow) double[ND*FTLE_BlockSize];
      
      if(XYZvals == NULL) {
        cout << "Error allocating memory for XYZvals in output" << endl;
        exit(1);
      }
      
      for(int d=0;d<ND;++d) {
        //set the denominator
        int den = 1;
        for(int dc=0;dc<d;++dc)
          den *= FTLE_Res[dc];
        
        if(Track_Method == 2) {
          double center_offset = 0;
          if (d<2) center_offset = 0.5 * (ftlemax[d] - ftlemin[d]);
          
          for(int ii=0;ii<FTLE_BlockSize;++ii)
            XYZvals[d*FTLE_BlockSize + ii] =  ((ii/den) % FTLE_Res[d])  * FTLE_Delta[d] - center_offset;
        }
        else {
          for(int ii=0;ii<FTLE_BlockSize;++ii) 
            XYZvals[d*FTLE_BlockSize + ii] = ftlemin[d] + ( (ii/den) % FTLE_Res[d]) * FTLE_Delta[d];
        }
      }
      
      
#if ND==2        
      if(Atmos_Set && Atmos_Spherical) {
        vector<double> Z(FTLE_BlockSize);
        for(int cc=0;cc<FTLE_BlockSize;++cc) {
          double theta = XYZvals[cc]*(M_PI/180.0);
          double phi = XYZvals[FTLE_BlockSize + cc]*(M_PI/180.0);
          XYZvals[cc]=Atmos_Radius*cos(theta)*cos(phi);
          XYZvals[FTLE_BlockSize + cc] = Atmos_Radius*sin(theta)*cos(phi);
          Z[cc] = Atmos_Radius*sin(phi);
        }
        TECDAT110(&FTLE_BlockSize,&XYZvals[0],&DIsDouble);
        TECDAT110(&FTLE_BlockSize,&XYZvals[FTLE_BlockSize],&DIsDouble);
        TECDAT110(&FTLE_BlockSize,&Z[0],&DIsDouble);
      }
      else {
        TECDAT110(&FTLE_BlockSize,&XYZvals[0],&DIsDouble);
        TECDAT110(&FTLE_BlockSize,&XYZvals[FTLE_BlockSize],&DIsDouble);
      }
#elif ND==3
      if(Atmos_Set && Atmos_Spherical) {
        for(int cc=0;cc<FTLE_BlockSize;++cc) {
          double theta = XYZvals[cc]*(M_PI/180.0);
          double phi = XYZvals[FTLE_BlockSize + cc]*(M_PI/180.0);
          double r = Atmos_Radius+Atmos_AltScale*XYZvals[2*FTLE_BlockSize + cc];
          XYZvals[cc] = r*cos(theta)*cos(phi);
          XYZvals[FTLE_BlockSize + cc] = r*sin(theta)*cos(phi);
          XYZvals[2*FTLE_BlockSize + cc] = r*sin(phi);
        }
      }
      TECDAT110(&FTLE_BlockSize,&XYZvals[0],&DIsDouble);
      TECDAT110(&FTLE_BlockSize,&XYZvals[FTLE_BlockSize],&DIsDouble); 
      TECDAT110(&FTLE_BlockSize,&XYZvals[2*FTLE_BlockSize],&DIsDouble); 
#elif ND==4
      TECDAT110(&FTLE_BlockSize,&XYZvals[0],&DIsDouble);
      TECDAT110(&FTLE_BlockSize,&XYZvals[FTLE_BlockSize],&DIsDouble); 
      TECDAT110(&FTLE_BlockSize,&XYZvals[2*FTLE_BlockSize],&DIsDouble); 
      TECDAT110(&FTLE_BlockSize,&XYZvals[3*FTLE_BlockSize],&DIsDouble); 
#elif ND==6
      TECDAT110(&FTLE_BlockSize,&XYZvals[0],&DIsDouble);
      TECDAT110(&FTLE_BlockSize,&XYZvals[1*FTLE_BlockSize],&DIsDouble); 
      TECDAT110(&FTLE_BlockSize,&XYZvals[2*FTLE_BlockSize],&DIsDouble); 
      TECDAT110(&FTLE_BlockSize,&XYZvals[3*FTLE_BlockSize],&DIsDouble); 
      TECDAT110(&FTLE_BlockSize,&XYZvals[4*FTLE_BlockSize],&DIsDouble); 
      TECDAT110(&FTLE_BlockSize,&XYZvals[5*FTLE_BlockSize],&DIsDouble); 
#endif
      
      
      delete []XYZvals;
    }
    
    double *field = new (nothrow) double[FTLE_BlockSize];
    if(field == NULL) {
      cout << "Error allocating memory for field in output" << endl;
      exit(1);
    }

    for(int nf=0;nf<LCS_NumFields;++nf) {
      
      fin.read((char *)&field[0],FTLE_BlockSize*sizeof(double));
      if(swap_endian) {
        for(register int ii=0;ii<FTLE_BlockSize;++ii) {
          BS(field[ii]);
        }
      }
      TECDAT110(&FTLE_BlockSize,field,&DIsDouble);
    }
    fin.close();
    
    delete []field;
    
    ++ss;
    
  }
  TECEND110();
  
  cout << "Output file is " << outfilename << endl;
}


void ConvertTimeToDate(double frametime,char *zonetitle) {

  char  *datestr;
  time_t timezerosec;
  time_t timenowsec;
  
  struct tm timezero;
 
  timezero.tm_year = Time_Origin[0]-1900;
  timezero.tm_mon = Time_Origin[1]-1;
  timezero.tm_mday = Time_Origin[2];
  timezero.tm_hour = Time_Origin[3];
  timezero.tm_min = Time_Origin[4];
  timezero.tm_sec = Time_Origin[5];
  timezero.tm_isdst = 1;

  timezerosec=mktime(&timezero);
  
  timenowsec=timezerosec+(time_t)(floor((frametime)*3600));
  datestr = strtok(ctime(&timenowsec),"\n");
  sprintf(zonetitle,"%s GMT",datestr);
  
}


void Format_TraceData() {
  
  fin.open(rawfilename,ios::binary);
  fin.read((char *)&Compute_Type, sizeof(int) );
  fin.read((char *)&Time_Origin,6*sizeof(int) );
  fin.read((char *)&Output_TRes, sizeof(int) );
  fin.read((char *)&Atmos_Set, sizeof(int) );
  fin.read((char *)&Atmos_Radius, sizeof(double) );
  
  if(swap_endian) {
    BS(Compute_Type);
     for(int d=0;d<6;++d)
        BS(Time_Origin[d]);
    BS(Output_TRes);
    BS(Atmos_Set);
    BS(Atmos_Radius);
  }
  
  /*cout << Compute_Type << endl;
  cout << Output_TRes << endl;
  cout << Atmos_Set << endl;
  cout << Atmos_Radius << endl;*/
  
  sprintf(outfilename,"%s.plt",rawfilename);
#if ND==2
  if(Atmos_Set) {
    if(Atmos_Spherical) 
      TECINI110("Drifter Trajectories","X Y Z Longitude Latitude color", outfilename, Path_Output, &Debug, &VIsDouble);
    else   
      TECINI110("Drifter Trajectories","Longitude Latitude color", outfilename, Path_Output, &Debug, &VIsDouble);
  }
  else {
    TECINI110("Drifter Trajectories","X Y color", outfilename, Path_Output, &Debug, &VIsDouble);
  }
#elif ND==3
  if(Atmos_Set) {
    if(Atmos_Spherical)
      TECINI110("Drifter Trajectories","X Y Z  Longitude Latitude Altitude color", outfilename, Path_Output, &Debug, &VIsDouble);
    else      
      TECINI110("Drifter Trajectories","Longitude Latitude Altitude color", outfilename, Path_Output, &Debug, &VIsDouble);
  }
  else {
    TECINI110("Drifter Trajectories","X Y Z color", outfilename, Path_Output, &Debug, &VIsDouble);
  }
#elif ND==4
  TECINI110("Drifter Trajectories","X Y Z W color", outfilename, Path_Output, &Debug, &VIsDouble);
#elif ND==6
  TECINI110("Drifter Trajectories","X Y Z Xd Yd Zd color", outfilename, Path_Output, &Debug, &VIsDouble);
#endif
  
  
  for(int tt=0;tt<Output_TRes;++tt) {
    
    int numdrifters;
    double frametime;  
    /* Read in header data */
    fin.read((char *)&numdrifters, sizeof(int) );
    fin.read((char *)&frametime, sizeof(double) );
    
    /* Do endian change if required */
    if(swap_endian) {
      BS(numdrifters);
      BS(frametime);
    }
    
    int blocksize = ND*numdrifters;
    
    X.resize(blocksize,0);
    color.resize(numdrifters,0);
    
    fin.read((char *)&X[0], blocksize*sizeof(double));
    fin.read((char *)&color[0], numdrifters*sizeof(double));
    
    if(swap_endian) {
      for(register int ii=0;ii<blocksize;++ii) {
        BS(X[ii]);
      }
      for(register int ii=0;ii<numdrifters;++ii) {
        BS(color[ii]);
      }
    }
    
    char zonetitle[LONGSTRING];
    if(Date_Flag) {
      ConvertTimeToDate(frametime,zonetitle);
    }
    else {
      sprintf(zonetitle,"T = %f",frametime);
    }
    cout << numdrifters << " drifters at " << zonetitle << endl;
    
    TECZNE110(zonetitle, &ZoneType, &numdrifters, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax, &frametime,
              &StrandID, &ParentZn, &IsBlock, &NFConns, &FNMode, NULL, NULL, sharing, &ShrConn);
    
#if ND==2        
    if(Atmos_Set && Atmos_Spherical) {
      vector<double> Xc(numdrifters);
      vector<double> Yc(numdrifters);
      vector<double> Zc(numdrifters);
      for(int cc=0;cc<numdrifters;++cc) {
        double theta = X[cc]*(M_PI/180.0);
        double phi = X[numdrifters + cc]*(M_PI/180.0);
        Xc[cc] = Atmos_Radius*cos(theta)*cos(phi);
        Yc[cc] = Atmos_Radius*sin(theta)*cos(phi);
        Zc[cc] = Atmos_Radius*sin(phi);
      }
      TECDAT110(&numdrifters,&Xc[0],&DIsDouble);
      TECDAT110(&numdrifters,&Yc[0],&DIsDouble);
      TECDAT110(&numdrifters,&Zc[0],&DIsDouble);
    }

    TECDAT110(&blocksize,&X[0],&DIsDouble);
    
#elif ND==3
    if(Atmos_Set && Atmos_Spherical) {
      vector<double> Xc(numdrifters);
      vector<double> Yc(numdrifters);
      vector<double> Zc(numdrifters);
      for(int cc=0;cc<numdrifters;++cc) {
        double theta  =   X[cc]*(M_PI/180.0);
        double phi    =   X[numdrifters + cc]*(M_PI/180.0);
        double r      =   Atmos_Radius + X[2*numdrifters + cc]; 
        Xc[cc] =  r*cos(theta)*cos(phi);
        Yc[cc] =  r*sin(theta)*cos(phi);
        Zc[cc] =  r*sin(phi);
      }
      TECDAT110(&numdrifters,&Xc[0],&DIsDouble);
      TECDAT110(&numdrifters,&Yc[0],&DIsDouble);
      TECDAT110(&numdrifters,&Zc[0],&DIsDouble);
    }
    
    TECDAT110(&blocksize,&X[0],&DIsDouble);  

#elif ND==4
    TECDAT110(&blocksize,&X[0],&DIsDouble);
#elif ND==6
    TECDAT110(&blocksize,&X[0],&DIsDouble);
#endif
    TECDAT110(&numdrifters,&color[0],&DIsDouble);
  }
  
  TECEND110();
  
  fin.close();
  
  cout << "Output file is " << outfilename << endl;
  
  
}

void Format_PlotData() {
  
  int   blocksize;
  
  fin.open(rawfilename,ios::binary);
  fin.read((char *)&Compute_Type, sizeof(int) );
  fin.read((char *)&Time_Origin, 6*sizeof(int) );
  fin.read((char *)&Output_TRes, sizeof(int) );
  fin.read((char *)&Atmos_Set, sizeof(int) );
  fin.read((char *)&Atmos_Radius, sizeof(double) );
  fin.read((char *)&blocksize, sizeof(int) );
  fin.read((char *)&Plot_Res[0], ND * sizeof(int));
  fin.read((char *)&Plot_Min[0], ND * sizeof(double));
  fin.read((char *)&Plot_Max[0], ND * sizeof(double));
  fin.read((char *)&Plot_Delta[0], ND * sizeof(double));
  
  if(swap_endian) {
    BS(Compute_Type);
     for(int d=0;d<6;++d)
        BS(Time_Origin[d]);
    BS(Output_TRes);
    BS(Atmos_Set);
    BS(Atmos_Radius);
    BS(blocksize);
    for(int d=0;d<ND;++d) {
      BS(Plot_Res[d]);
      BS(Plot_Min[d]);
      BS(Plot_Max[d]);
      BS(Plot_Delta[d]);
    }
  }
  /*
  cout << Compute_Type << endl;
  cout << Output_TRes << endl;
  cout << Atmos_Set << endl;
  cout << Atmos_Radius << endl;
  cout << blocksize << endl;
  for(int d=0;d<ND;++d) {
    cout << Plot_Res[d] << endl;
    cout << Plot_Min[d] << endl;
    cout << Plot_Max[d] << endl; 
    cout << Plot_Delta[d] << endl;
  }*/
  
  int IMax = Plot_Res[0];
  int JMax = Plot_Res[1];
  int KMax = 1;
  
  sprintf(outfilename,"%s.plt",rawfilename);
#if ND==2
  TECINI110("Velocity Plot","X Y U V", outfilename, Path_Output, &Debug, &VIsDouble);
#elif ND==3  
  TECINI110("Velocity Plot","X Y Z U V W", outfilename, Path_Output, &Debug, &VIsDouble);
  KMax = Plot_Res[2];
#endif  
  
  U = new double[ND*blocksize];
  
  for(int tt=0;tt<Output_TRes;++tt) {
    
    double frametime; 
    /* Read in header data */
    fin.read((char *)&frametime, sizeof(double) );
    
    /* Do endian change if required */
    if(swap_endian) BS(frametime);
    
    fin.read((char *)&U[0], ND*blocksize*sizeof(double));
    
    if(swap_endian) {
      for(register int ii=0;ii<(ND*blocksize);++ii) {
        BS(U[ii]);
      }
    }
    
    if(tt>0) {
      for(int d=0;d<ND;++d)
        sharing[d]=1;
    }
    
    char zonetitle[LONGSTRING];
    if(Date_Flag) {
      ConvertTimeToDate(frametime,zonetitle);
    }
    else {
      sprintf(zonetitle,"T = %f",frametime);
    }

    cout << blocksize << " plot points at " << zonetitle << endl;
    
    TECZNE110(zonetitle, &ZoneType, &IMax, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax, &frametime,
              &StrandID, &ParentZn, &IsBlock, &NFConns, &FNMode, NULL, NULL, sharing, &ShrConn);
    
    if(tt == 0) {
      double *XYZvals = new (nothrow) double[blocksize];
      if(XYZvals == NULL) {
        cout << "Error allocating memory for XYZvals in output" << endl;
        exit(1);
      }
      
      for(int d=0;d<ND;++d) {
        //set the denominator
        int den = 1;
        for(int dc=0;dc<d;++dc)
          den *= Plot_Res[dc];
        
        for(int ii=0;ii<blocksize;++ii) 
          XYZvals[ii] = Plot_Min[d] + ( (ii/den) % Plot_Res[d]) * Plot_Delta[d];
        
        TECDAT110(&blocksize,XYZvals,&DIsDouble);
      }
      delete []XYZvals;
    }
    
    int datasize = ND*blocksize;
    TECDAT110(&datasize,&U[0],&DIsDouble);
  }
  
  TECEND110();
  
  fin.close();
  
  cout << "Output file is " << outfilename << endl;
  
}

#include <sys/stat.h> 

bool FileExists(string strFilename) { 
  struct stat stFileInfo; 
  bool blnReturn; 
  int intStat; 
  
  // Attempt to get the file attributes 
  intStat = stat(strFilename.c_str(),&stFileInfo); 
  if(intStat == 0) { 
    // We were able to get the file attributes 
    // so the file obviously exists. 
    blnReturn = true; 
  } else { 
    // We were not able to get the file attributes. 
    // This may mean that we don't have permission to 
    // access the folder which contains this file. If you 
    // need to do that level of checking, lookup the 
    // return values of stat which will give you 
    // more details on why stat failed. 
    blnReturn = false; 
  } 
  
  return(blnReturn); 
}

void ByteSwap(unsigned char * b, int n)
{
  register int i = 0;
  register int j = n-1;
  while (i<j)
  {
    swap(b[i], b[j]);
    i++, j--;
  }
}
