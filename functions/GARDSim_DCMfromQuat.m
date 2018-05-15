function DCM = GARDSim_DCMfromQuat(quat)
% DCM from Quaternion
% Generates the DCM from quaternion vector
% Syntax: DCM = GARDSim_DCMfromQuat([quat])
%
%
% Written by Duncan Greer May 2006
%
% $Id: GARDSim_DCMfromQuat.m 1863 2008-07-14 07:02:29Z greerd $

warning('This function is depreciated. Please use GARD_QuatToDCM instead');
DCM = GARD_QuatToDCM(quat);




    
   