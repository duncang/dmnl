function TMatrix = T_ECEF2ENU(Longitude,Latitude);
% TMatrix = T_ECEF2ENU(Longitude,Latitude)
%
% generate transfer matrix ECEF to LTP coordinates (East-North-Up)
%
% $Id: T_ECEF2ENU.m 1863 2008-07-14 07:02:29Z greerd $
%


clong = cos(Longitude);
clat = cos(Latitude);
slong = sin(Longitude);
slat = sin(Latitude);



TMatrix(1,1) = -slong*clat;
TMatrix(1,2) =  clong*clat;
TMatrix(1,3) =  slat;

TMatrix(2,1) = -clong*slat;
TMatrix(2,2) = -slong*slat;
TMatrix(2,3) = clat;

TMatrix(3,1) = clong*clat;
TMatrix(3,2) = slong*clat;
TMatrix(3,3) = slat;

