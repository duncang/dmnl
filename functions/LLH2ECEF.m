function [Position] = LLH2ECEF(Latitude,Longitude,Height);
% [Position] = LLH2ECEF(Latitude,Longitude,Height)
% LLH to ECEF derived from functions in Chapter 4 of 'Geodesy' by Wolfgang
% Torge
% Version 1.00
% By Duncan Greer
% Converts from LLH to ECEF
% INPUT
% Latitude (rad)  %note - must be in radians!!
% Longitude (rad)
% Height (m)
%========================================================================
% OUTPUT
% Position = [X, Y, Z] (metres)
%========================================================================


% WGS-84 ellipsoid parameters
a = 6378137.0; % semi-major axis
b = 6356752.3142; % semi-minor axis

% calculate the first eccentricity
e = sqrt(a^2 - b^2) / a;



N = a^2 / sqrt(a^2 * cos(Latitude)^2 + b^2 * sin(Latitude)^2);

X = (N + Height) * cos(Latitude) * cos (Longitude);
Y = (N + Height) * cos(Latitude) * sin (Longitude);
Z = ((1 - e^2)*N + Height) * sin(Latitude);

Position = [X, Y, Z];

%This function has been validated and it is correct.