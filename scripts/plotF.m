clear all
clc
clf

load FTLEOutput


if ND==2
    figure(1)
    for tt=1:Output_TRes
        contourf(X{1},X{2},F(:,:,1))
        getframe;
    end
elseif ND==3
    figure(1)
    for tt=1:Output_TRes
        for zz=1:31
            contourf(X{1},X{2},F(:,:,zz,tt))
            getframe;
        end
    end

elseif ND==4
    figure(1)
    for tt=1:Output_TRes
        contourf(X{1},X{2},squeeze(F(:,:,3,tt)))
        getframe;
    end

elseif ND==6
    figure(1)
    for tt=1:Output_TRes
        contourf(X{1},X{2},squeeze(F(:,:,1,1,1,1,1,tt)))
        getframe;
    end
end