



DataFileName = 'data\NG05JA05.SP3';

[NumberSVs, ValidData, NumberEpochs, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data] = readsp3(DataFileName);




data = SV_Y_Data(1:90,1);





%X is the values of the polynomail of order BestFitSim


%each point is 15 minutes, so split up into 15*60 seconds = 900 secs


%for j = 1
% datafixed(1:900) = data(1);
% 
% 
% for j = 2:10
% for i =  901*(j-1) : 900*j + (j-1)
%     
%     datafixed(i) = data(j);
% 
% end
% end


%best fit is 11th order polynomial

%approximate between each poitn with 11th order poly


[CurveFitObs, BestFitSim,Xest] = LeastSquaresBestFit(data', 5, 1);





%now fit a curve through the points at the 15 minute marks
%can't use least squares best fit because want to pass through the points exactly.



data = SV_Y_Data(1:90,1);
x = 1:90;
y = data;
xx = 1:0.1:90;
yy = spline(x,y,xx);
















