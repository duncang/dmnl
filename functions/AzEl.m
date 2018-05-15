function [Azimuth, Elevation] = AzEl(Xu, Su)
%
% [Azimuth, Elevation] = AzEl(Xu, Su)
%
% returns the Azimuth and Elevation of teh Satellite given by Su to the
% user position Xu.
%
% $Id: AzEl.m 1874 2008-07-15 04:42:16Z n2523710 $
%


% get user geodetic latitude and longitude
[phi_u,lamda_u,h_u] = ECEF2LLH(Xu);


% calculate satellite line of sight vector
SV_LOS_ecef = (Su - Xu)';

% generate rotation matrix
TMatrix_ECEF2NEU = T_ECEF2NEU(lamda_u,phi_u);

% convert to LTP coordinates
SV_LOS_neu = TMatrix_ECEF2NEU * SV_LOS_ecef;

% calculate azimuth and elevation
Azimuth = atan2(SV_LOS_neu(2),SV_LOS_neu(1)); % radians
Elevation = atan2(SV_LOS_neu(3),sqrt(SV_LOS_neu(1)^2 + SV_LOS_neu(2)^2)); % radians


%note - corrected the Azimuth calculation 07.06.06 TB. it was
%atan2(SV_LOS_neu(1),SV_LOS_neu(2)) which is incorrect.
