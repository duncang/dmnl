function PlotGPSBestXYZ(gps_data)
% function PlotGPSBestXYZ(gps_data)
% Written by Duncan Greer 18 June 2007
% $Id: PlotGPSBestXYZ.m 853 2007-06-18 03:26:32Z greerd $
%

n = length(gps_data.PosECEF);
lat = zeros(1,n);
long = zeros(1,n);
height = zeros(1,n);
for i=1:length(gps_data.PosECEF)
    [lat(i),long(i),height(i)] = ECEF2LLH(gps_data.PosECEF(i,:));
end
    
figure();plot(long*180/pi,lat*180/pi,'r');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
figure();plot(height,'r');
