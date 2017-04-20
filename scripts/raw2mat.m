clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%

jobid = 2;

if jobid == 1
    ND = 3;
    endianformat = 'b';
    rawfilename = '/Users/philipdutoit/Documents/Research/LCS/Code/Newmanv1/examples/titan3D/Output/titandrift.raw';
elseif jobid == 2
    ND = 3;
    endianformat = 'l';
    rawfilename = '/Users/pdutoit/Documents/Research/LCS/Examples/Sigmoid/Output/FTLErep';
elseif jobid == 5
    ND = 3;
    endianformat = 'l';
    rawfilename = '/Users/philipdutoit/Documents/Research/LCS/Examples/ROMS3D/Output/roms3ddrifterout.raw';
elseif jobid == 3
    ND = 2;
    endianformat = 'b';
    rawfilename = '/Users/philipdutoit/Documents/Research/LCS/Examples/Overflow/Output/overflowdrift.raw';
elseif jobid == 4
    ND = 4;
    endianformat = 'b';
    rawfilename = '/Users/philipdutoit/Documents/Research/LCS/Examples/Mercury/Output/drifterout.raw';
elseif jobid == 6
    ND = 6;
    endianformat = 'l';
    rawfilename = '/Users/philipdutoit/Desktop/Output/drifterout.raw';
elseif jobid == 7
    ND = 6;
    endianformat = 'l';
    rawfilename = '/Users/philipdutoit/Desktop/Output/FTLErep';
elseif jobid == 8
    ND = 2;
    endianformat = 'l';
    rawfilename = '/Users/philipdutoit/Desktop/Output/FTLErep';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Track_Flag = 0;

if exist(rawfilename,'file')

    fin = fopen(rawfilename,'rb',endianformat);
    Compute_Type = fread(fin,1,'int');
    fclose(fin);

else
    ftlefilename = sprintf('%s0000.raw',rawfilename);
    if exist(ftlefilename,'file')==2

        Compute_Type = 0;

    else
        fprintf('Could not open: \n%s\n%s\n',rawfilename,ftlefilename);

    end
end




switch Compute_Type
    case 0,
        fprintf('Creating FTLE data\n');

        ftlefilename = sprintf('%s0000.raw',rawfilename);
        fin = fopen(ftlefilename,'rb',endianformat);
        Compute_Type = fread(fin,1,'int');
        Time_Origin = fread(fin,6,'int')
        FrameTime = fread(fin,1,'double')
        Output_TRes = fread(fin,1,'int')
        Atmos_Set = fread(fin,1,'int')
        Atmos_Radius = fread(fin,1,'double')
        Slide_Number = fread(fin,1,'int')
        Track_Storm = fread(fin,1,'int')
        ftlemin = fread(fin,ND,'double')
        ftlemax = fread(fin,ND,'double')
        FTLE_Res = fread(fin,ND,'int')'
        LCS_NumFields = fread(fin,1,'int')
        fclose(fin);

        FTLE_BlockSize = prod(FTLE_Res);
        XYblock = prod(FTLE_Res(1:2));
        if numel(FTLE_Res > 2)
            Zblock = prod(FTLE_Res(3:end));
        else
            Zblock = 1;
        end

        F = zeros(FTLE_Res(2),FTLE_Res(1),Zblock,LCS_NumFields,Output_TRes);

        ss = 1;
        while 1

            ftlefilename = sprintf('%s%04d.raw',rawfilename,ss-1);

            if ~exist(ftlefilename,'file')
                break
            end

            fprintf('%s\n',ftlefilename);
            fin = fopen(ftlefilename,'rb',endianformat);

            Compute_Type = fread(fin,1,'int');
            Time_Origin = fread(fin,6,'int');
            FrameTime = fread(fin,1,'double');
            Output_TRes = fread(fin,1,'int');
            Atmos_Set = fread(fin,1,'int');
            Atmos_Radius = fread(fin,1,'double');
            Slide_Number = fread(fin,1,'int');
            Track_Storm = fread(fin,1,'int');
            ftlemin = fread(fin,ND,'double');
            ftlemax = fread(fin,ND,'double');
            FTLE_Res = fread(fin,ND,'int')';
            LCS_NumFields = fread(fin,1,'int');

            for nf=1:LCS_NumFields
                for nb=1:Zblock
                    fdata = fread(fin,XYblock,'double');
                    fdata = reshape(fdata,FTLE_Res(1),FTLE_Res(2));
                    F(:,:,nb,nf,ss) = fdata';
                end
            end
            fclose(fin);
            ss = ss+1;

        end

        if ND==2
            F=squeeze(reshape(F,[FTLE_Res(2),FTLE_Res(1),LCS_NumFields,Output_TRes]));
        else
            F=squeeze(reshape(F,[FTLE_Res(2),FTLE_Res(1),FTLE_Res(3:end),LCS_NumFields,Output_TRes]));
        end
        X = cell(ND,1);
        for d=1:ND
            X{d} = linspace(ftlemin(d),ftlemax(d),FTLE_Res(d));
        end

        save FTLEOutput X F ND Output_TRes;

        fprintf('\nData is stored in FTLEOutput.mat\n');
        
        

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1,

        fprintf('Creating Trace file\n');


        fin = fopen(rawfilename,'rb',endianformat);
        Compute_Type = fread(fin,1,'int')
        Time_Origin = fread(fin,6,'int');
        Output_TRes = fread(fin,1,'int')
        Atmos_Set = fread(fin,1,'int');
        Atmos_Radius = fread(fin,1,'double');


        X = cell(Output_TRes,1);
        for tt=1:Output_TRes

            numdrifters = fread(fin,1,'int');
            FrameTime = fread(fin,1,'double');

            blocksize = ND*numdrifters;

            xdata = fread(fin,blocksize,'double');
            xdata = reshape(xdata,numdrifters,ND);
            color = fread(fin,numdrifters,'double');

            X{tt} = [xdata color];
            fprintf('%d drifters at time %f\n',numdrifters,FrameTime);
        end


        fclose(fin);
        
        save TraceOutput X ND Output_TRes;
        
        fprintf('\nData is stored in TraceOutput.mat\n');

        
        
    case 2,
        fprintf('Creating Velocity file : NOT YET SUPPORTED\n');
        

    otherwise,
        fprintf('ERROR: Incorrect Compute_Type = %d\n',Compute_Type);
end



clear





