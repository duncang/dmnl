function [PR_Pred,PRR_Pred] = GARD_GPSMeasurementPrediction(UserPos_LLHC,UserVel_NEDC,SVPos,SVVel)
% Carries out pseudorange and pseudorange rate prediction
% Written by Duncan Greer 4 June 2008
% $Id$



% load constants used in this function
GPSConstants;


Tecef2ned= T_ECEF2NED(UserPos_LLHC(1),UserPos_LLHC(2)); 
Tned2ecef = Tecef2ned';

% find apriori estimate of pseudorange for each sigma point

%Get User position in ECEF
UserPos = LLH2ECEF(UserPos_LLHC(1),UserPos_LLHC(2),UserPos_LLHC(3));
UserPos(4) = UserPos_LLHC(4);

UserVel =  Tned2ecef * UserVel_NEDC(1:3);
UserVel(4) = UserVel_NEDC(4);

geo_range_to_sat = sqrt((SVPos(1) - UserPos(1))^2 + (SVPos(2) - UserPos(2))^2 + (SVPos(3) - UserPos(3))^2);

geo_vel_to_sat = (SVVel(1) - UserVel(1))*(SVPos(1)-UserPos(1)) + ...
                 (SVVel(2) - UserVel(2))*(SVPos(2)-UserPos(2)) + ...
                 (SVVel(3) - UserVel(3))*(SVPos(3)-UserPos(3));

Relative_Velocity = geo_vel_to_sat/geo_range_to_sat;

delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(1) *UserPos(2) - SVPos(2) * UserPos(1));


%% calculate measurement predicition
% pseudorange prediction

PR_Pred = geo_range_to_sat  - delta_pr_omegaedot + UserPos(4) - SVPos(4);


%predicted relative velocity of sv and receiver
PRR_Pred = Relative_Velocity + UserVel(4) - SVVel(4);


