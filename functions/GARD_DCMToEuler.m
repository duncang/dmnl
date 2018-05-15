function [phi,theta,psi] = GARD_DCMToEuler(DCM);
% Euler anglers from DCM 

% Syntax: DCM = GARD_DCMToEuler(phi,theta,psi)
%
%
%Troy B. 5.2.09
%%needs verifying



phi = atan2(DCM(3,2), DCM(3,3)); % rotation about x-axis = roll
theta = asin(-DCM(3,1));         % rotation about y-axis = pitch
psi = atan2(DCM(2,1), DCM(1,1)); % rotation about z-axis = yaw
            
