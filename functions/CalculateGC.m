function [Bearing, Distance] = CalculateGC(lat1,lon1,lat2,lon2);
% function [Bearing, Distance] = CalculateGC(LatA,LongA,LatB,LongB);
% Calculates the great-circle bearing and distance between two points, A
% and B given by geodetic latitude and longitude coordinates (in Radians).
% Written by Duncan Greer January 2006
% Equations from Kayton page 40
%
% $Id: CalculateGC.m 2544 2009-05-08 03:40:04Z bruggema $
%

% load wgs84 constants
WGS84Constants;

% calculate Rg, gaussian radius of curvature
Rm = MeridianRadius(lat1);
Rp = PrimeRadius(lat1);
Rg = sqrt(Rm*Rp);

% cosDonRg = (sin(LatA)*sin(LatB) + cos(LatA)*cos(LatB)*cos(LongA - LongB));
% Distance = acos(cosDonRg) * Rg;
% 
% sinBt = (cos(LatB)/sin(Distance/Rg))*sin(LongA - LongB);
% Bearing = asin(sinBt);
% Bearing = 2*pi - Bearing; % for some reason i need to maek this correction 
% 
% % apply 2pi bound
% Bearing = mod(Bearing,2*pi);

% equations from http://williams.best.vwh.net/avform.htm

d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2));
Distance = d * Rg;

temp = (sin(lat2)-sin(lat1)*cos(d))/(sin(d)*cos(lat1));

if abs(temp) > 1   
    temp = sign(temp);    %i.e. make it either -1 or 1 but not greater than.   
end 


if( sin(lon2-lon1)<0)       
   tc1=acos(temp);  
else
   tc1=2*pi-acos(temp);    
end

%Fix by Troy B. with acos if  ABS(x) > 1.0 then it results in complex!! due to
%computational error there was 1.0000000000056 so it gave complex!! so make
%sure it is not larger than 1
%acos(1) 

Bearing = mod(tc1,2*pi);
Bearing = 2*pi - Bearing;
