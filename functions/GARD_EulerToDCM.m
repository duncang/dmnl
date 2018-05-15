function DCM = GARD_EulerToDCM(phi,theta,psi);
% DCM from Euler Angles
% generates the direction cosine matrix from euler angles
% Syntax: DCM = GARD_EulerToDCM(phi,theta,psi)
%
%
% Written by Duncan Greer May 2006
%
% $Id: GARD_EulerToDCM.m 1850 2008-07-14 04:52:47Z greerd $

c11 = cos(theta)*cos(psi);
c12 = -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi);
c13 = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
c21 = cos(theta)*sin(psi);
c22 = cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi);
c23 = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
c31 = -sin(theta);
c32 = sin(phi)*cos(theta);
c33 = cos(phi)*cos(theta);



DCM = [c11, c12, c13; c21, c22, c23; c31, c32, c33];