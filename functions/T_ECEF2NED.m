function TMatrix = T_ECEF2NED(Latitude,Longitude);
% generate transfer matrix ECEF to LTP (North-East-Down) coordinates
% from page 202 Kayton
%
% TMatrix = T_ECEF2NED(Latitude,Longitude)
%
% Note. Latitude and Longitude in radians
%
% See also
% http://www.mathworks.com/access/helpdesk/help/toolbox/aeroblks/directioncosinematrixeceftoned.html
%
% $Id: T_ECEF2NED.m 1863 2008-07-14 07:02:29Z greerd $
%

if abs(Latitude) > pi/2
    argumentString = sprintf('%f rad = %f degrees', Latitude, Latitude * 180/pi);
    warning(strcat('[T_ECEF2NED] Latitude seems large - Check arguments: ', argumentString));
end

if abs(Longitude) > pi
    argumentString = sprintf('%f rad = %f degrees', Longitude, Longitude * 180/pi);
    warning(strcat('[T_ECEF2NED] Longitude seems large - Check arguments: ', argumentString));
end


clong = cos(Longitude);
clat = cos(Latitude);
slong = sin(Longitude);
slat = sin(Latitude);

TMatrix(1,1) = -clong*slat;
TMatrix(1,2) = -slong*slat;
TMatrix(1,3) =  clat;

TMatrix(2,1) = -slong;
TMatrix(2,2) = clong;
TMatrix(2,3) = 0;

TMatrix(3,1) = -clong*clat;
TMatrix(3,2) = -slong*clat;
TMatrix(3,3) = -slat;
