function [Rn, Re] = WGS84_calcRnRe(Latitude);
% Function to find the Rn and Re 
% Written by Troy 
% where Latitude is in radians.
% this function has been verified

% load WGS84 constants
WGS84Constants;



    %Rn - North -south , Meridian
    Rn = a * (1 - e2) / (sqrt(1 - e2 * sin(Latitude)^2))^3;

    %find the normal radius of curvature, east-west, Re, (Prime)
    Re = a / sqrt(1 - e2 * sin(Latitude)^2);


