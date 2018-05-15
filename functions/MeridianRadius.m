function RM = MeridianRadius(Latitude);
% Function to find the meridian radius of curvature
% Written by Duncan Greer 22 January 2006
%
% $Id: MeridianRadius.m 1883 2008-07-15 05:53:55Z n2523710 $
%
% Usage:
%    RM = MeridianRadius(Latitude)
% where Latitude is in radians.
%

% load WGS84 constants
WGS84Constants;

RM = a * (1 - e2) / (sqrt(1 - e2 * sin(Latitude)^2))^3;