function [pseudo_roll,gps_flight_path_angle,gps_track] = GARD_PseudoAttitude(velocity,acceleration)
% function [pseudo_roll,gps_flight_path_angle,gps_track] = GARD_PseudoAttitude(velocity,acceleration)
% Written by Duncan Greer 12 June 2007
%
% Inputs
% ======
% velocity - NED velocity vector (1x3)
% acceleration - NED acceleration vector (1x3)
%
% Outputs
% =======
% pseudo_roll - 
% gps_flight_path_angle - 
% gps_track - 
%
% $Id: GARD_PseudoAttitude.m 858 2007-06-19 23:04:27Z greerd $
%

gravity_vec = [0;0;-9.79];

% calculate flight path angle
gps_flight_path_angle = atan2(-velocity(3),sqrt(velocity(1)^2+velocity(2)^2)); 
gps_track = atan2(velocity(2),velocity(1));

a_tilde = acceleration' - ((acceleration*velocity'/norm(velocity)^2) * velocity');
g_tilde =  gravity_vec - (gravity_vec'*velocity'/norm(velocity)^2)'* velocity';

pseudo_lift = a_tilde - g_tilde;
p_tilde = cross(g_tilde,velocity);

pseudo_roll = asin(pseudo_lift'*p_tilde' / (norm(pseudo_lift) * norm(p_tilde)));

