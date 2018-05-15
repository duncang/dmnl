function DCM = GARDSim_DCMfromEuler(phi,theta,psi);
% DCM from Euler Angles
% generates the direction cosine matrix from euler angles
% Syntax: DCM = GARDSim_DCMfromEuler(phi,theta,psi)
%
%
% Written by Duncan Greer May 2006
%
% $Id: GARDSim_DCMfromEuler.m 1863 2008-07-14 07:02:29Z greerd $

disp('Warning! GARDSim_DCMfromEuler is DEPRECIATED!  Use GARD_EulerToDCM Instead');
DCM = GARD_EulerToDCM(phi,theta,psi);
