




%Constants


c = 2.99792458e8;


L1_f = 1575.42e6; %Hz
L2_f = 1227.6e6; %Hz

gamma = (L1_f/L2_f)^2;  % unitless

L1_Wavelength = Speedoflight/L1_f; %Metres



%This is for the plotting and analysis of results from the Repeater Test done in September 2005




%compare positions


PosTypeRoof = T1PosSess1Roof;
PosTypeRoom = T1PosSess1Room;


Size = size(PosTypeRoof);


%Size(1) = 600 used this for test 1 session 2 cos there was a lot of zeros



%Convert to Lat Long HHgth

for i = 1:Size(1)


[LatRoof(i) LonRoof(i) HgtRoof(i)] = ECEF2LLH([PosTypeRoof(i,1),PosTypeRoof(i,2),PosTypeRoof(i,3)]);

end


LatRoof = LatRoof*180/pi;
LonRoof = LonRoof*180/pi;


for i = 1:Size(1)


[LatRoom(i) LonRoom(i) HgtRoom(i)] = ECEF2LLH([PosTypeRoom(i,1),PosTypeRoom(i,2),PosTypeRoom(i,3)]);

end



LatRoom = LatRoom*180/pi;
LonRoom = LonRoom*180/pi;


%plot the height difference between the two calculated positions, is easier to understand it then.

HgtDiff = HgtRoom(100:Size(1)) - HgtRoof(100:Size(1));


plot(HgtDiff); title 'Room - Roof Height (m)'; xlabel('Time (Secs)');


pause;

meanheightDiff = mean( HgtRoom(100:Size(1)) - HgtRoof(100:Size(1)))


%plot the values
%plot X

%did it this way to scale the y axis of the subplots to the same values


plot1 = PosTypeRoof(100:Size(1),1)- mean(PosTypeRoof(100:Size(1),1));
plot2 = PosTypeRoom(100:Size(1),1)- mean(PosTypeRoom(100:Size(1),1));


plot3 = PosTypeRoof(100:Size(1),2)- mean(PosTypeRoof(100:Size(1),2));
plot4 = PosTypeRoom(100:Size(1),2)- mean(PosTypeRoom(100:Size(1),2));


plot5 = PosTypeRoof(100:Size(1),3)- mean(PosTypeRoof(100:Size(1),3));
plot6 = PosTypeRoom(100:Size(1),3)- mean(PosTypeRoom(100:Size(1),3));


min1 = min(plot1);
min2 = min(plot2);
max1 = max(plot1);
max2 = max(plot2);

axesvaluemin = min([min1 min2]) - 5;
axesvaluemax = max([max1 max2]) + 5;

subplot(2,1,1); plot(plot1); title 'Norm X Pos Roof (m)'; axis([1 Size(1)+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(plot2); title 'Norm X Pos Room (m)'; xlabel('Time (Secs)'); axis([1 Size(1)+10 axesvaluemin axesvaluemax ]);

pause;


min1 = min(plot3);
min2 = min(plot4);
max1 = max(plot3);
max2 = max(plot4);

axesvaluemin = min([min1 min2]) - 5;
axesvaluemax = max([max1 max2]) + 5;

subplot(2,1,1); plot(plot3); title 'Norm Y Pos Roof (m)'; axis([1 Size(1)+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(plot4); title 'Norm Y Pos Room (m)'; xlabel('Time (Secs)');axis([1 Size(1)+10 axesvaluemin axesvaluemax ]);

pause;

min1 = min(plot5);
min2 = min(plot6);
max1 = max(plot5);
max2 = max(plot6);

axesvaluemin = min([min1 min2]) - 5;
axesvaluemax = max([max1 max2]) + 5;

subplot(2,1,1); plot(plot5); title 'Norm Z Pos Roof (m)'; axis([1 Size(1)+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(plot6); title 'Norm Z Pos Room (m)'; xlabel('Time (Secs)');axis([1 Size(1)+10 axesvaluemin axesvaluemax ]);

pause;
subplot(2,1,1); plot(PosTypeRoof(100:Size(1),4)/c); title 'dt Pos Roof (Secs)'; 
subplot(2,1,2); plot(PosTypeRoom(100:Size(1),4)/c); title 'dt Pos Room (Secs)'; xlabel('Time (Secs)');

pause;


clf;
%plot the differences

XPosDiff = PosTypeRoom(100:Size(1),1) - PosTypeRoof(100:Size(1),1);
YPosDiff = PosTypeRoom(100:Size(1),2) - PosTypeRoof(100:Size(1),2);
ZPosDiff = PosTypeRoom(100:Size(1),3) - PosTypeRoof(100:Size(1),3);
TPosDiff = PosTypeRoom(100:Size(1),4)/c - PosTypeRoof(100:Size(1),4)/c;


plot(XPosDiff); title 'X Pos Room - X Pos Roof (m)'; xlabel('Time (Secs)');
pause;
plot(YPosDiff); title 'Y Pos Room - Y Pos Roof (m)'; xlabel('Time (Secs)');
pause;
plot(ZPosDiff); title 'Z Pos Room - Z Pos Roof (m)'; xlabel('Time (Secs)');
pause;
plot(TPosDiff); title 'dt Pos Room - dt Pos Roof (Secs)'; xlabel('Time (Secs)');
pause;


  
    
%to show the first 150 seconds with test 1 session 1 - can see the uncorrected clock;

 

    
subplot(2,1,1); plot(PosTypeRoof(1:150,4)/c); title 'dt Pos Roof (Secs)'; 
subplot(2,1,2); plot(PosTypeRoom(1:150,4)/c); title 'dt Pos Room (Secs)'; xlabel('Time (Secs)');
    




meanXRoof = meanNZ(PosTypeRoof(100:Size(1),1)')
meanYRoof = meanNZ(PosTypeRoof(100:Size(1),2)')
meanZRoof = meanNZ(PosTypeRoof(100:Size(1),3)')
meanTRoof = meanNZ(PosTypeRoof(100:Size(1),4)'/c)

meanXRoom = meanNZ(PosTypeRoom(100:Size(1),1)')
meanYRoom = meanNZ(PosTypeRoom(100:Size(1),2)')
meanZRoom = meanNZ(PosTypeRoom(100:Size(1),3)')
meanTRoom = meanNZ(PosTypeRoom(100:Size(1),4)'/c)

meanXDiff = mean(XPosDiff)
meanYDiff = mean(YPosDiff)
meanZDiff = mean(ZPosDiff)
meanTDiff = mean(TPosDiff)







for i = 100:Size(1)
    NormDiffXYZ(i) = sqrt( (PosTypeRoom(i,1) - PosTypeRoof(i,1))^2 + (PosTypeRoom(i,2) - PosTypeRoof(i,2))^2 + (PosTypeRoom(i,3) - PosTypeRoof(i,3))^2 );
end

%room - roof as a range

plot(NormDiffXYZ); title 'Range Room-Roof XYZ(m)'; xlabel('Time (Secs)');
    

    


NormDiffXYZ = sqrt(meanXDiff^2 + meanYDiff^2 +meanZDiff^2)

%including the time component:
NormDiffXYZT = sqrt(meanXDiff^2 + meanYDiff^2 +meanZDiff^2 +(meanTDiff*c)^2)










%compare with Ashtechs own position solution for the room  (this is for test 2 session 3)


Xashcompare = data(100:Size(1),3) - PosTypeRoom(100:Size(1),1);


Yashcompare = data(100:Size(1),4) - PosTypeRoom(100:Size(1),2);

Zashcompare = data(100:Size(1),5) - PosTypeRoom(100:Size(1),3);


plot(Xashcompare); title 'X Pos Room Ash - X Pos Room GARDSim (m)'; xlabel('Time (Secs)');
pause;


plot(Yashcompare); title 'Y Pos Room Ash - Y Pos Room GARDSim (m)'; xlabel('Time (Secs)');
pause;

plot(Zashcompare); title 'Z Pos Room Ash - Z Pos Room GARDSim (m)'; xlabel('Time (Secs)');

pause;




%--------------------------------------------------

%Compare pseudoranges

%Test 1 Session1
 

%subtract the differences from each other.

% for Sat = 1:32
% difference = (PRDiffNoClk(Sat-(Sat+1),100:Size(1)))
% end




PRTypeRoof =  T1PRSess1Roof;
PRTypeRoom =  T1PRSess1Room;

CPTypeRoof =  T1CPSess1Roof;
CPTypeRoom =  T1CPSess1Room;











%----------------------------------------
%Multipath analysis: Code - Carrier. 


%Code minus carrier for Room:
clear CMCminusNRoom CMCminusNRoof SignalRoom SignalRoof BestFitRoom BestFitRoof;


PRRoof = PRTypeRoof;
CPRoof = CPTypeRoof*L1_Wavelength;  %don't forget to convert this to metres like i did at first!!


for Sat = 1:32
    for i = 1:Size(1)
        
        if (PRRoof(Sat,i) ~= 0)
            
            CMCRoof(Sat,i) = PRRoof(Sat,i) - CPRoof(Sat,i) ;
        else
            CMCRoof(Sat,i) = 0;
        end
    end
end





SignalRoof = CMCRoof(29,5100:6200) ;
OrderRoof = 1;
FindBestFitRoof = 1;


[CurveFitRoof, BestFitRoof] = LeastSquaresBestFit(SignalRoof, OrderRoof, FindBestFitRoof);


CMCminusNRoof = SignalRoof - CurveFitRoof; %this is the multipath + 2*iono etc






PRRoom = PRTypeRoom;
CPRoom = CPTypeRoom*L1_Wavelength;  %don't forget to convert this to metres like i did at first!!


for Sat = 1:32
    for i = 1:Size(1)
        
        if (PRRoom(Sat,i) ~= 0)
            
            CMCRoom(Sat,i) = PRRoom(Sat,i) - CPRoom(Sat,i) ;
        else
            CMCRoom(Sat,i) = 0;
        end
    end
end






SignalRoom = CMCRoom(29,5100:6200) ;
OrderRoom = 1;
FindBestFitRoom = 1;


[CurveFitRoom, BestFitRoom] = LeastSquaresBestFit(SignalRoom, OrderRoom, FindBestFitRoom);



CMCminusNRoom = SignalRoom - CurveFitRoom; %this is the multipath + 2*iono etc



plot(CMCminusNRoom(:)); title('Code - Carrier (m) SV 29'); xlabel('Time (Secs)');
hold
plot(CMCminusNRoof(:),'r'); title('Code - Carrier (m) SV 29'); xlabel('Time (Secs)');
 

pause;




%takes too long for it to do the above calculations so do it for juts one or two.



%-----------------------------------------------------------


%analyse the pseudoranges


%SV 26

% StartTime = 1;
% 
% EndTime = Size(1);

for Sat = 1:32
    for i = 1:Size(1)
      
    PRDiff(Sat,i) = PRTypeRoom(Sat,i) - PRTypeRoof(Sat,i);
    
    end
end


%Subtract the receiver clock bias from the measurements



 
for Sat = 1:32
    for i = 1:Size(1)
        
        if (PRTypeRoom(Sat,i) ~= 0 & PRTypeRoof(Sat,i) ~=0)
            
            PRDiffNoClk(Sat,i) = (PRTypeRoom(Sat,i)-PosTypeRoom(i,4)) - (PRTypeRoof(Sat,i)-PosTypeRoof(i,4));
        else
            PRDiffNoClk(Sat,i) = 0;
        end
    end
end



plot(PosTypeRoom(5475:5480,4)

PosTypeRoof(5475:5480,4));


   
subplot(2,1,1); plot(PosTypeRoof(5475:6003,4)); title 'dt Pos Roof (Secs)'; 


plot(diff(PosTypeRoom(5450:6003,4))); title 'dt dot Pos Room (Secs/sec)'; xlabel('Time (Secs)');






%do plotting now


for Sat = 1:32
   
    if max(PRDiffNoClk(Sat,1:Size(1))) ~= 0
          
       %plot(PRDiffNoClk(Sat,100:Size(1))); title (['Room-Roof Pseudorange no Rx clock (m) SV ', num2str(Sat)]); xlabel('Time (Secs)');
       plot(PRDiffNoClk(Sat,100:Size(1))); title (['Room-Roof Pseudorange no Rx clock (m) SV ', num2str(Sat)]); xlabel('Time (Secs)');
       %hold on;
       
      
  
pause;
end
end


pause;






%diff(PRTypeRoom(9,:)


%look at quality of the PR itself
% 
% SignalRoof = (PRTypeRoom(9,:));
% OrderRoof = 1;
% FindBestFitRoof = 1;
% 
% 
% [CurveFitRoof, BestFitRoof] = LeastSquaresBestFit(SignalRoof, OrderRoof, FindBestFitRoof);
% 
% 
% uu = SignalRoof-CurveFitRoof;
% 
% plot(uu); title (['Pseudorange residual SV 9']); xlabel('Time (Secs)');
% 
% 
% 
% 
% 
% SignalRoof = diff(PRTypeRoom(3,:));
% OrderRoof = 1;
% FindBestFitRoof = 1;
% 
% 
% [CurveFitRoof, BestFitRoof] = LeastSquaresBestFit(SignalRoof, OrderRoof, FindBestFitRoof);
% 
% 
% uu = SignalRoof-CurveFitRoof;
% 
% plot(uu); title (['Pseudorange residual SV 3']); xlabel('Time (Secs)');
% 
% 
% 



%---------------------------------------------------------------
%Test effect of adding a constant bias to the PR measurements (for BiasAddingTestingNoBias.mat workspace)



difference = PosWithoutPRbias-PosWithbias;



plot1 = PosWithoutPRbias(1:3999,1)- mean(PosWithoutPRbias(1:3999,1));
plot2 = PosWithbias(1:3999,1)- mean(PosWithbias(1:3999,1));


plot3 = PosWithoutPRbias(1:3999,2)- mean(PosWithoutPRbias(1:3999,2));
plot4 = PosWithbias(1:3999,2)- mean(PosWithbias(1:3999,2));


plot5 = PosWithoutPRbias(1:3999,3)- mean(PosWithoutPRbias(1:3999,3));
plot6 = PosWithbias(1:3999,3)- mean(PosWithbias(1:3999,3));


min1 = min(plot1);
min2 = min(plot2);
max1 = max(plot1);
max2 = max(plot2);

axesvaluemin = min([min1 min2]) - 5;
axesvaluemax = max([max1 max2]) + 5;

subplot(2,1,1); plot(plot1); title 'Norm X Pos No Bias (m)'; axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(plot2); title 'Norm X Pos Bias (m)'; xlabel('Time (Secs)'); axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;


min1 = min(plot3);
min2 = min(plot4);
max1 = max(plot3);
max2 = max(plot4);

axesvaluemin = min([min1 min2]) - 5;
axesvaluemax = max([max1 max2]) + 5;

subplot(2,1,1); plot(plot3); title 'Norm Y Pos No Bias (m)'; axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(plot4); title 'Norm Y Pos Bias (m)'; xlabel('Time (Secs)');axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;

min1 = min(plot5);
min2 = min(plot6);
max1 = max(plot5);
max2 = max(plot6);

axesvaluemin = min([min1 min2]) - 5;
axesvaluemax = max([max1 max2]) + 5;

subplot(2,1,1); plot(plot5); title 'Norm Z Pos No Bias (m)'; axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(plot6); title 'Norm Z Pos Bias (m)'; xlabel('Time (Secs)');axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
subplot(2,1,1); plot(PosWithoutPRbias(1:3999,4)); title 'dt Pos No Bias (m)'; 
subplot(2,1,2); plot(PosWithbias(1:3999,4)); title 'dt Pos Bias (m)'; xlabel('Time (Secs)');

pause;


clf;
%plot the differences

XPosDiff = PosWithoutPRbias(1:3999,1) - PosWithbias(1:3999,1);
YPosDiff = PosWithoutPRbias(1:3999,2) - PosWithbias(1:3999,2);
ZPosDiff = PosWithoutPRbias(1:3999,3) - PosWithbias(1:3999,3);
TPosDiff = PosWithoutPRbias(1:3999,4) - PosWithbias(1:3999,4);


plot(XPosDiff); title 'X Pos No Bias - X Pos Bias(m)'; xlabel('Time (Secs)');
pause;
plot(YPosDiff); title 'Y Pos No Bias - Y Pos Bias(m)'; xlabel('Time (Secs)');
pause;
plot(ZPosDiff); title 'Z Pos No Bias - Z Pos Bias(m)'; xlabel('Time (Secs)');
pause;
plot(TPosDiff); title 'dt Pos No Bias - dt Pos Bias(m)'; xlabel('Time (Secs)');
pause;