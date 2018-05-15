function [Latitude,Longitude,Height] = ECEF2LLH(Position)
% This has been validated by TB 28 October 2005 Note: function output in radians not degrees
% [Latitude,Longitude,Height] = ECEF2LLH(Position)
%
% $Id: ECEF2LLH.m 1874 2008-07-15 04:42:16Z n2523710 $
%
global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;

X = Position(1);
Y = Position(2);
Z = Position(3);

% WGS-84 ellipsoid parameters
a = 6378137.0; % semi-major axis
b = 6356752.3142; % semi-minor axis


P = sqrt(X^2 + Y^2);

Theta = atan((Z * a) / (P * b));

Esq = 1 - b^2 / a^2;
EPsq = a^2 / b^2 - 1;

Latitude = atan(((Z + EPsq * b * sin(Theta)^3)) / (P- Esq * a * cos(Theta)^3));
Longitude = atan2(Y,X);
n = a^2 / sqrt(a^2 * cos(Latitude)^2 + b^2 * sin(Latitude)^2);
Height = P / cos(Latitude) - n;
