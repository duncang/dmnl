function TMatrix = T_ECEF2LTP(Longitude,Latitude);
% generate transfer matrix ECEF to LTP coordinates (NEU)
% from page 202 Kayton

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

TMatrix(3,1) = clong*clat;
TMatrix(3,2) = slong*clat;
TMatrix(3,3) = slat;
