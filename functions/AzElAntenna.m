function [AzimuthBody, ElevationBody] = AzElAntenna(Xu, Su,PHI,THETA,PSI)
% returns the Azimuth and Elevation of teh Satellite given by Su to the
% user position Xu relative to aircraft body frame i.e. aircraft's antenna
%useful for determining whether a satellite signal will be blocked by
%antenna or not due to aircraft turn. 
%by Troy Bruggemann

% get user geodetic latitude and longitude
[phi_u,lamda_u,h_u] = ECEF2LLH(Xu);


% calculate satellite line of sight vector
SV_LOS_ecef = (Su - Xu)';

% generate rotation matrix
TMatrix_ECEF2NEU = T_ECEF2NEU(lamda_u,phi_u);

% convert to LTP coordinates
SV_LOS_neu = TMatrix_ECEF2NEU * SV_LOS_ecef;


SV_LOS_ned = [SV_LOS_neu(1),SV_LOS_neu(2), -SV_LOS_neu(3)]';

TMatrix_Body2NED = T_Body2NED(PHI,THETA, PSI);

SV_LOS_body = TMatrix_Body2NED'*SV_LOS_ned;


% calculate azimuth and elevation
AzimuthBody = atan2(SV_LOS_body(2),SV_LOS_body(1)); % radians
ElevationBody = atan2(-SV_LOS_body(3),sqrt(SV_LOS_body(1)^2 + SV_LOS_body(2)^2)); % radians %note needs to be -SV_LOS_body(3) to convert it to 'up' direction.


%note - corrected the Azimuth calculation 07.06.06 TB. it was
%atan2(SV_LOS_neu(1),SV_LOS_neu(2)) which is incorrect.


% AzimuthBody = atan2(SV_LOS_neu(2),SV_LOS_neu(1)); % radians
% ElevationBody = atan2(SV_LOS_neu(3),sqrt(SV_LOS_neu(1)^2 + SV_LOS_neu(2)^2)); % radians