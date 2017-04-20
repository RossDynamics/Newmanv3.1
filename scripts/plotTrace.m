clear
clc

load TraceOutput


if ND == 2
    
    for tt=1:Output_TRes
        figure(1)
        plot(X{tt}(:,1),X{tt}(:,2),'.r');
        %axis([21 29 2 9])
        getframe;
    end

elseif ND == 3
    figure(1)
    for tt=1:Output_TRes

        plot3(X{tt}(:,1),X{tt}(:,2),X{tt}(:,3),'.r');
        %axis([21 29 2 9])
        getframe;
    end

elseif ND == 4
    figure(1)
    for tt=1:Output_TRes
        plot(X{tt}(:,1),X{tt}(:,2),'.r');
        hold all
    end
elseif ND == 6
    figure(1)
    for tt=1:Output_TRes
        plot(X{tt}(:,1),X{tt}(:,2),'.r');
        hold all
    end

end