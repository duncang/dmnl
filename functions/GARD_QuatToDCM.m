function DCM = GARD_QuatToDCM(quat)
% DCM from Quaternion
% Generates the DCM from quaternion vector
% Syntax: DCM = GARD_QuatToDCM(quat)
%
%
% Written by Duncan Greer May 2006
%
% $Id: GARDSim_DCMfromQuat.m 708 2006-12-20 01:17:28Z n2523710 $

a = quat(1);
b = quat(2);
c = quat(3);
d = quat(4);

c11 = (a^2 + b^2 - c^2 - d^2);
c12 = 2 * (b*c - a*d);
c13 = 2 * (b*d + a*c);
c21 = 2 * (b*c + a*d);
c22 = (a^2 - b^2 + c^2 - d^2);
c23 = 2 * (c*d - a*b);
c31 = 2 * (b*d - a*c);
c32 = 2 * (c*d + a*b);
c33 = (a^2 - b^2 - c^2 + d^2);

DCM = [c11, c12, c13; c21, c22, c23; c31, c32, c33];



    
   