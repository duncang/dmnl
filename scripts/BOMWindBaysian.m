




%read in wind data



%curve fit to the data, determine best order



%subtract polynomial approximation from truth, then analyse the  noise
%statistics.


%read data file



Filename = 'data.txt';


filename_obsAMicro = fopen(Filename,'r');

TimeVec = [4, 5, 6,10,11,15, 22];
%for 15.00

k =  1;
j = 1;


Direction(1,1) = 0;
Speed(1,1) = 0;



Hgt = 5;  %5 ft



for i = 1:8000

    obslineMicro = fgetl(filename_obsAMicro);


    for m = 1:length(TimeVec)

        Time = TimeVec(m);

        if (str2num(obslineMicro(22:23)) == Time & str2num(obslineMicro(74:79)) == Hgt)

            if obslineMicro(53:55) == ' '

                Direction(m,k) = Direction(m,k-1);
            else

                Direction(m,k) = str2num(obslineMicro(53:55));
                
                if m ==7
                k = k+1;
                end
            end

            if obslineMicro(46:49) == ' '

                Speed(m,j) = Speed(m,j-1);
            else

                Speed(m,j) = str2num(obslineMicro(46:49));
                
                
                if m ==7
                j = j+1;

                end
            end
        else    %if obslineMicro(

            %do nothing
            % Direction(i) =  Direction(i-1);
            %Speed(i) = Speed(i-1);
        end

        %end  %if (obslineMicro(22:23) ~= ' ' &

        %end %& obslineMicro(76:79) ~= ' ')
    end  % for m = 1:length(TimeVec)

end

fclose(filename_obsAMicro);


plot(Speed(1,:),'*');
hold;

plot(Speed(2,:),'r*');

plot(Speed(3,:),'g*');
plot(Speed(4,:),'y*');
plot(Speed(5,:),'k*');
plot(Speed(6,:),'c*');
plot(Speed(7,:),'m*');




%remove zeros
SpeedNoZero = Speed;
DirectionNoZero=Direction;

for m = 1:7
    
    for i = 1:length(SpeedNoZero)
        
        if SpeedNoZero(m,i) == 0 && DirectionNoZero(m,i) == 0 && i~=1    %if both are zero then its likely theres no data
            SpeedNoZero(m,i) = SpeedNoZero(m,i-1);
            
               DirectionNoZero(m,i) =  DirectionNoZero(m,i-1);
            
        end   
      
      
    end
    
end




plot(Speed(1,:),'-');
hold;

plot(Speed(2,:),'r-');

plot(Speed(3,:),'g-');
plot(Speed(4,:),'y-');
plot(Speed(5,:),'k-');
plot(Speed(6,:),'c-');
plot(Speed(7,:),'m-');



plot(Direction(1,:),'-');
hold;

plot(Direction(2,:),'r-');

plot(Direction(3,:),'g-');
plot(Direction(4,:),'y-');
plot(Direction(5,:),'k-');
plot(Direction(6,:),'c-');
plot(Direction(7,:),'m-');



startplot = 2;
endplot = 25;

plot(SpeedNoZero(1,startplot:endplot),'-');
hold;
%plot(SpeedNoZero(2,:),'r-');

plot(SpeedNoZero(3,startplot:endplot),'g-');
plot(SpeedNoZero(4,startplot:endplot),'y-');
%plot(SpeedNoZero(5,:),'k-');
%plot(SpeedNoZero(6,:),'c-');
plot(SpeedNoZero(7,startplot:endplot),'m-');





plot(DirectionNoZero(1,startplot:endplot),'-');
hold;

%plot(DirectionNoZero(2,:),'r-');

plot(DirectionNoZero(3,startplot:endplot),'g-');
plot(DirectionNoZero(4,startplot:endplot),'y-');
%plot(DirectionNoZero(5,:),'k-');
%plot(DirectionNoZero(6,:),'c-');
plot(DirectionNoZero(7,startplot:endplot),'m-');



%find zeros

find(SpeedNoZero(6,:) == 0)

find(DirectionNoZero(4,:) ==0)



SpeedPlot(1,:) = SpeedNoZero(1,startplot:endplot);
SpeedPlot(2,:) = SpeedNoZero(3,startplot:endplot);
SpeedPlot(3,:) = SpeedNoZero(4,startplot:endplot);
SpeedPlot(4,:) = SpeedNoZero(7,startplot:endplot);



DirectionPlot(1,:) = DirectionNoZero(1,startplot:endplot);
DirectionPlot(2,:) = DirectionNoZero(3,startplot:endplot);
DirectionPlot(3,:) = DirectionNoZero(4,startplot:endplot);
DirectionPlot(4,:) = DirectionNoZero(7,startplot:endplot);




plot(SpeedPlot(1,:),'-');
hold;
%plot(SpeedNoZero(2,:),'r-');

plot(SpeedPlot(2,:),'g-');
plot(SpeedPlot(3,:),'y-');
%plot(SpeedNoZero(5,:),'k-');
%plot(SpeedNoZero(6,:),'c-');
plot(SpeedPlot(4,:),'m-');





plot(DirectionPlot(1,:),'-');
hold;

%plot(DirectionNoZero(2,:),'r-');

plot(DirectionPlot(2,:),'g-');
plot(DirectionPlot(3,:),'y-');
%plot(DirectionNoZero(5,:),'k-');
%plot(DirectionNoZero(6,:),'c-');
plot(DirectionPlot(4,:),'m-');


[CurveFit, BestFit] = LeastSquaresBestFitWind(SpeedPlot, 10, 0);


[CurveFit, BestFit] = LeastSquaresBestFitWind(DirectionPlot, 10, 0);




plot(CurveFit,'r','LineWidth',4) 




[CurveFit, Param,limit] = LeastSquaresBestFitWindFourier(DirectionPlot);



[CurveFit, BestFit,  Limit] = LeastSquaresBestFitWindFourierNOrder(SpeedPlot(:,:), 3, 0);

