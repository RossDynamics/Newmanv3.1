# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:00:06 2016

Based on raw2mat.m by Philip du Toit From Newmanv3.1

Reads Raw Binary files and converts them into a Matlab file.
Each Binary file corresponds to one time slice of data.

@author: Peter Nolan, Virginia Tech, Department of 
Biomedical Engineering And Mechanics
"""
from struct import *
import numpy as np
import os
import sys
import h5py
import hdf5storage

'!!!!!! PARAMETERS !!!!!!'

'number of spacial dimensions'
nd = 3

'Name of the sequence of files, example FTLErep for FTLErepXXXX.raw'
#rawfilename = 'FTLErep'
#rawfilename = 'Output/FTLEatt'
rawfilename = 'Output/drift.raw'
'encoding, l for little-endian, b for big-endian'
endianformat = 'l'


'!!!!!!! BEGIN CODE !!!!!!!'

Track_Flag = 0

'Check and set endian format'
if endianformat == 'l':
    ef = '<'
elif endianformat == 'b':
    ef = '>'
else:
    sys.exit('ERROR: Unrecognized endian format.')


'Find the first file and Set the Compute_Type'
if os.path.isfile(rawfilename):
    f = open(rawfilename,'rb')
    Compute_Type = np.array(unpack(ef+'i',f.read(4)))
    f.close
else:
    ftlefilename = '%s0000.raw' % (rawfilename)
    if os.path.isfile(ftlefilename): 
        Compute_Type =0
    else:
        print 'Could not open: %s or %s' % (rawfilename,ftlefilename)


if Compute_Type == 0:
    
    print 'Creating FTLE data'
    
    'Open first in sequence file'
    f = open(rawfilename+'0000.raw','rb')
    
    'Decode and Read the file header'
    Compute_Type = np.array(unpack(ef+'i',f.read(4)))
    Time_Origin = np.array(unpack(ef+'i'*6,f.read(4*6)))
    FrameTime = np.array(unpack(ef+'d',f.read(8)))
    Output_TRes = np.array(unpack(ef+'i',f.read(4)))
    Atmos_Set = np.array(unpack(ef+'i',f.read(4)))
    Atmos_Radius = np.array(unpack(ef+'d',f.read(8)))
    Slide_Number = np.array(unpack(ef+'i',f.read(4)))
    Track_Storm = np.array(unpack(ef+'i',f.read(4)))
    
    'ftlemin and ftleman are the min and max of the'
    'various spacial dimensions. Not the min and max'
    'of the FTLE field over the domain'
    ftlemin = np.array(unpack(ef+'d'*nd,f.read(8*nd)))
    ftlemax = np.array(unpack(ef+'d'*nd,f.read(8*nd)))
    FTLE_Res = np.array(unpack(ef+'i'*nd,f.read(4*nd)))
    LCS_NumFields = np.array(unpack(ef+'i',f.read(4)))
    
    'Close file'
    f.close()
    
    'Output file header for user'
    print 'Compute Type = ' + str(Compute_Type[0])
    print 'Time Origin = ' + str(Time_Origin)
    print 'Frame Time = ' + str(FrameTime[0])
    print 'Output TRes = ' + str(Output_TRes[0])
    print 'Atmos Set = ' + str(Atmos_Set[0])
    print 'Atmos Radius = ' + str(Atmos_Radius[0])
    print 'Slide Number = ' + str(Slide_Number[0])
    print 'Track Storm = ' + str(Track_Storm[0])
    print 'ftlemin = ' + str(ftlemin)
    print 'ftlemax = ' + str(ftlemax)
    print 'FTLE Res = ' + str(FTLE_Res)
    print 'LCS NumFields = ' + str(LCS_NumFields[0])
    
    #FTLE_Res=np.abs(FTLE_Res)
    'Set the total blocksize of the spacial FTLE domain, length*width*height...'
    FTLE_BlockSize = np.prod(FTLE_Res)
    
    'Set the blocksize for just the xy portion of the spacial domain, length*width'
    XYblock = np.prod(FTLE_Res[0:2])
    
    'Set the size of the non-xy portion of the spacial domain'
    'If the data is 2D then Zblock = 1, otherwise Zblock is equal'
    'to the product of the non-xy domain spacial resolutions'
    if FTLE_Res.size > 2:
        Zblock = np.prod(FTLE_Res[2:FTLE_Res.size])
    else:
        Zblock = 1
    
    print FTLE_Res[0]
    print FTLE_Res[1]

           
    'Initialize a 5 dimensional numpy array of zeros'
    'First index is the y domain, width'
    'Second index is the x domain, lenth'
    'Third index is the remainder of the spacial domain, height...'
    'Fourth index is the size of the LCS NumFields'
    'Fifth index is the size of the temporal domain, time'
    
    F = np.zeros((FTLE_Res[1],FTLE_Res[0],Zblock,LCS_NumFields[0],Output_TRes[0]))
    
    'Initialize ss as 0, where ss is the iteration of the following while loop'
    'ss is also index of the RAW file we want to read'
    ss = 0
    'Initialize an INFINITE loop to read all of the binary files in the sequence'
    while 1: 
        'Set the name of the file for this iteration to the'
        'rawfilename concatonated with the current'
        'iteration of the loop, with a minimum of four digits'
        ftlefilename = '%s%04d.raw' % (rawfilename, ss)
        
        'Check that the file actually exists, if it doesnt break the loop'
        if not(os.path.isfile(ftlefilename)):
            break
    
        'Output the filename for the user'
        print ftlefilename
        
        'Open the file'
        f = open(ftlefilename,'rb')
        
        'Read the header'
        Compute_Type = np.array(unpack(ef+'i',f.read(4)))
        Time_Origin = np.array(unpack(ef+'i'*6,f.read(4*6)))
        FrameTime = np.array(unpack(ef+'d',f.read(8)))
        Output_TRes = np.array(unpack(ef+'i',f.read(4)))
        Atmos_Set = np.array(unpack(ef+'i',f.read(4)))
        Atmos_Radius = np.array(unpack(ef+'d',f.read(8)))
        Slide_Number = np.array(unpack(ef+'i',f.read(4)))
        Track_Storm = np.array(unpack(ef+'i',f.read(4)))
        ftlemin = np.array(unpack(ef+'d'*nd,f.read(8*nd)))
        ftlemax = np.array(unpack(ef+'d'*nd,f.read(8*nd)))
        FTLE_Res = np.array(unpack(ef+'i'*nd,f.read(4*nd)))
        LCS_NumFields = np.array(unpack(ef+'i',f.read(4)))
        
        'loop through all of the LCS NumFields'
        for nf in range(0,LCS_NumFields):
            'loop through the Zblock'
            for nb in range(0,Zblock):
                'for each LCS NumField and each Zblock, read the xy data'
                'store the xy data in a vector'
                fdata = np.array(unpack(ef+'d'*XYblock,f.read(8*XYblock)))
                'reshape that vector into an x by y matrix using FORTRAN ordering'
                fdata = np.reshape(fdata,(FTLE_Res[0],FTLE_Res[1]),order='F')
                'save the transpose of that matrix to the F array'
                'for the current NumField, Zblock and Time entry'
                F[:,:,nb,nf,ss]=np.transpose(fdata)
        'Close file'                
        f.close()
        'Set ss for next iteration'
        ss = ss + 1
    
    'Prepare F to be saved in Matlab format'
    if nd==2:
        'If the data was only 2D then reshape F removing the Zblock'
        'using FORTRAN ordering, then squeeze out singular dimentions'
        F=np.squeeze(np.reshape(F,(FTLE_Res[1],FTLE_Res[0],LCS_NumFields[0],Output_TRes[0]),order='F'))
    else:
        'If the data was more than 2D then reshape F around the Zblock'
        'using FORTRAN ordering, then squeeze out singular dimentions'
        F=np.squeeze(np.reshape(F,(FTLE_Res[1],FTLE_Res[0],FTLE_Res[2:FTLE_Res.size],LCS_NumFields[0],Output_TRes[0]),order='F'))
    
    'Create Excel cells of all the X and Y coordinates for Matlab'
    X=[0]*nd
    for d in range(0,nd):
        X[d] = np.linspace(ftlemin[d],ftlemax[d],FTLE_Res[d])
        
    'Save Results'
    'hdf5storage allows for saving data to matlab v7.3 files'
    matfiledata = {}
    matfiledata[u'X']=X
    del X
    matfiledata[u'F']=F
    del F
    matfiledata[u'ND']=nd
    matfiledata[u'Output_TRes']=Output_TRes
    hdf5storage.write(matfiledata, '.', 'FTLEOutput.mat', matlab_compatible=True)
    
    'Notify user that results are saved'
    print 'Data is stored in FTLEOutput.mat'

elif Compute_Type == 1:
    print 'Notice: this functionality has not been tested!'
    print 'Creating Trace file'
    sys.stdout=file('output.txt','w')
    'Open the file'
    f = open(rawfilename,'rb')
        
    'Read the header'
    Compute_Type = np.array(unpack(ef+'i',f.read(4)))
    Time_Origin = np.array(unpack(ef+'i'*6,f.read(4*6)))
    Output_TRes = np.array(unpack(ef+'i',f.read(4)))
    Atmos_Set = np.array(unpack(ef+'i',f.read(4)))
    Atmos_Radius = np.array(unpack(ef+'d',f.read(8)))
    
    X=[0]*Output_TRes[0]
    
    'Loop through each time step'
    for tt in range(0,Output_TRes[0]):
        
        numdrifters = np.array(unpack(ef+'i',f.read(4)))
        FrameTime = np.array(unpack(ef+'d',f.read(8)))
        
        blocksize = nd*numdrifters[0]
        
        xdata = np.array(unpack(ef+'d'*blocksize,f.read(8*blocksize)))
        xdata = np.reshape(xdata,(numdrifters[0],nd),order='F')
        color = np.array(unpack(ef+'d'*numdrifters[0],f.read(8*numdrifters[0])))

        X[tt] = np.concatenate([xdata[:,0], xdata[:,1], xdata[:,2], color])
        print '%d drifters at time %f.' % (numdrifters, FrameTime)
    
    f.close()

    'Save Results'
    'hdf5storage allows for saving data to matlab v7.3 files'
    matfiledata = {}
    matfiledata[u'X']=X
    matfiledata[u'ND']=nd
    matfiledata[u'Output_TRes']=Output_TRes
    hdf5storage.write(matfiledata, '.', 'TraceOutput.mat', matlab_compatible=True)    
        
    print 'Data is stored in TraceOutput.mat'
    print 'Notice: this functionality has not been tested!'   
    
elif Compute_Type == 2:
    print 'Creating Velocity file : NOT YET SUPPORTED'

else:
    print 'ERROR: Incorrect Compute_Type = %d' % (Compute_Type)