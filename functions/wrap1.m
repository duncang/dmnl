%-------------------------------------------------------------------------
% WRAP1: Wraps a signal between min and max values
% -------------------------------------------------------------------------
%author Lennon Cork 
% [y] = wrap1(x,xmin,xmax)

%note : TB:

%If angle is 270 degrees, the resultant angle will be -90. 


function [y] = wrap1(x,xmin,xmax)
y = (mod((x-xmin),(xmax-xmin))+xmin);