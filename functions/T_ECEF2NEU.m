function TMatrix = T_ECEF2NEU(Longitude,Latitude);
% generate transfer matrix ECEF to LTP coordinates (North-East-Up)
% This is taken from Kayton  p 202 and has been validated. 28 Oct 05
% Note.  Longitude and Latitude in radians.
%
% $Id: T_ECEF2NEU.m 1863 2008-07-14 07:02:29Z greerd $
%

clong = cos(Longitude);
clat = cos(Latitude);
slong = sin(Longitude);
slat = sin(Latitude);

%this is the code that is correct
TMatrix(1,1) = -clong*slat;
TMatrix(1,2) = -slong*slat;
TMatrix(1,3) =  clat;

TMatrix(2,1) = -slong;
TMatrix(2,2) = clong;
TMatrix(2,3) = 0;

TMatrix(3,1) = clong*clat;
TMatrix(3,2) = slong*clat;
TMatrix(3,3) = slat;






%
% TMatrix(1,1) = -slong*clat;
% TMatrix(1,2) =  clong*clat;
% TMatrix(1,3) =  slat;
% 
% TMatrix(2,1) = -clong*slat;
% TMatrix(2,2) = -slong*slat;
% TMatrix(2,3) = clat;
% 
% TMatrix(3,1) = clong*clat;
% TMatrix(3,2) = slong*clat;
% TMatrix(3,3) = slat;

