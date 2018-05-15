function stec = vtec2stec(vtec,z,R,R0,H)
% function stec = vtec2stec(vtec,z,R,R0,H)
% 
% Written by Duncan Greer 15 Aug 2005
%
% This function implements the modified single layer mapping function
% (MSLM) to convert from a vertical TEC (given by IONEX) to a slant TEC
% which can be used to calculate ionospheric delay and code advance effects
% on GPS signals.
%
%   stec is the slant TEC that is returned
%   vtec is the vertical TEC for the Ionospheric Pierce Point
%   z is the geocentric zenith distance at the height of the station
%   R is the radius of the earth at the station (height of the receiver?)
%   R0 is the radius of the earth (6371.8 km)
%   H is the ionospheric single layer height
%   
% References:
%  1) S. Schaer "How to use CODE's Global Ionosphere Maps" Astronomical
%  Institute, University of Berne, May 1997.

Fz = 1/cos(asin((R/(R0+H))*sin(z)));
stec = vtec * Fz;

