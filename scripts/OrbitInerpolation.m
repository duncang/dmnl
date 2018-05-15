% How to inertpolate SP3 satellite trajectories using Matlab's SPLINE
% Written by Troy Bruggemann and Duncan Greer
% $Id: OrbitInerpolation.m 1884 2008-07-15 05:54:33Z n2523710 $
%
%

SP3DataFileName = 'data/sp3/NG05JA05.SP3';

if ~exist('NumberSVs_SP3')
    [NumberSVs_SP3, VehicleIDs_SP3, NumberEpochs_SP3, GPSTime_SP3, SV_X_Data_SP3, SV_Y_Data_SP3, SV_Z_Data_SP3, SV_T_Data_SP3] = readsp3(SP3DataFileName);
end



data = SV_Y_Data_SP3(1:90,1);
x = 1:90;
y = data;
xx = 1:0.1:90;
yy = spline(x,y,xx);

plot(xx,yy);
hold on;
stem(x,y,'r');