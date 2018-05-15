function [PR_Vec_minus, PRR_Vec_minus] = GARD_GetMeasurementPrediction(x,NumberGPSMeasurements,SVPos,SVVel)
%
%
% This fucntion calculates the predicted pseudorange (prior) based on the
% predicted state (x).
% Written by Duncan Greer
% $Id: GARD_GetMeasurementPrediction.m 1850 2008-07-14 04:52:47Z greerd $
%
% INPUTS
% ======
% x - System State Vector - must consist of at least user position,
% velocity and clock
% 
%
%
%
%
% OUTPUTS
% =======
%
% PR_Vec_minus
% PRR_Vec_minus
%
%
%

GPSConstants;

PR_Vec_minus = zeros(NumberGPSMeasurements,1);
PRR_Vec_minus = zeros(NumberGPSMeasurements,1);


for k = 1:NumberGPSMeasurements
    Tecef2ned= T_ECEF2NED(x(1),x(2)); 
    Tned2ecef = Tecef2ned';

    %Get User position in ECEF
    UserPos = LLH2ECEF(x(1),x(2),x(3));
    UserPos(4) = x(17);

    UserVel =  Tned2ecef * x(4:6);
    UserVel(4) = x(18);

    geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
    geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + ...
                     (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + ...
                     (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
    delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) *UserPos(2) - SVPos(k,2) * UserPos(1));

    %% calculate measurement predicition
    % pseudorange prediction
    PR_Vec_minus(k) = geo_range_to_sat  - delta_pr_omegaedot(k) + UserPos(4);
    %predicted relative velocity of sv and receiver
    Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
    PRR_Vec_minus(k) = Relative_Velocity(k) + UserVel(4);



end  % for k = 1:NumberGPSMeasurements

