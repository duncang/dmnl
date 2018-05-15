function degrees = degmin2deg(degrees,minutes)
% function degrees = degmin2deg(degrees,minutes)
% convert degrees and decimal minute to decimal degrees
% written by Duncan Greer 24 January 2006
% $Id: degmin2deg.m 1874 2008-07-15 04:42:16Z n2523710 $

if(sign(degrees)==1)
    degrees = degrees + minutes/60;
else
    degrees = degrees - minutes/60;
end
