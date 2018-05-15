function RP = PrimeRadius(Latitude)
% Function to determine the prime radius of curvature of the earth
% Written by Duncan Greer 22 January 2006
%
% $Id: PrimeRadius.m 1883 2008-07-15 05:53:55Z n2523710 $
%
% Usage:
%     RP = PrimeRadius(Latitude)
% Where Latitude is provided in radians
%

% load WGS84 constants
WGS84Constants;


RP = a / sqrt(1 - e2 * sin(Latitude)^2);