function [ys_kminus_0, ys_kminus_i] = GARD_GenerateMeasurementPredictionSigmaPoints(Na,xs_0_k,xs_i_k,NumberGPSMeasurements,NumberStates,ProcessNoiseStates,SVPos,SVVel,UseMeasurements)

% load GPS constants
GPSConstants;


if ~exist('UseMeasurements','var')
   UseMeasurements = NumberGPSMeasurements; 
end

% get the a-proiri measuremnt prediction (y_minus)
for k = 1:UseMeasurements
    Tecef2ned= T_ECEF2NED(xs_0_k(1),xs_0_k(2)); 
    Tned2ecef = Tecef2ned';
    % find apriori estimate of pseudorange for each sigma point

    % zero-th sigma point
    %Get User position in ECEF
    UserPos = LLH2ECEF(xs_0_k(1),xs_0_k(2),xs_0_k(3));
    UserPos(4) = xs_0_k(17);

    UserVel =  Tned2ecef * xs_0_k(4:6);
    UserVel(4) = xs_0_k(18);

    geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
    geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + ...
                     (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + ...
                     (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
    delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) *UserPos(2) - SVPos(k,2) * UserPos(1));

    %% calculate measurement predicition
    % pseudorange prediction
    PR_Vec_minus_0(k) = geo_range_to_sat + UserPos(4) - delta_pr_omegaedot(k) - SVPos(k,4) + xs_0_k(NumberStates+ProcessNoiseStates+k);
    %predicted relative velocity of sv and receiver
    Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
    PRR_Vec_minus_0(k) = Relative_Velocity(k) + UserVel(4) + xs_0_k(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+k) - SVVel(k,4);

    for i=1:2*Na
        Tecef2ned= T_ECEF2NED(xs_i_k(1,i),xs_i_k(2,i));
        Tned2ecef = Tecef2ned';

       % Get User position in ECEF
        UserPos = LLH2ECEF(xs_i_k(1,i),xs_i_k(2,i),xs_i_k(3,i));
        UserPos(4) = xs_i_k(17,i);

        UserVel =  Tned2ecef * xs_i_k(4:6,i);
        UserVel(4) = xs_i_k(18,i);

        for m = 1:3
            ele(m) =  SVPos(k,m) - UserPos(m);
        end

        geo_range_to_sat =  norm(ele);

        %geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
        geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + ...
                         (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + ...
                         (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
        delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) *UserPos(2) - SVPos(k,2) * UserPos(1));

        PR_Vec_minus_i(k,i) = geo_range_to_sat + UserPos(4) - delta_pr_omegaedot(k)  - SVPos(k,4) + xs_i_k(NumberStates+ProcessNoiseStates+k,i);

        %predicted relative velocity of sv and receiver
        Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
        PRR_Vec_minus_i(k,i) = Relative_Velocity(k) + UserVel(4) - SVVel(k,4) + xs_i_k(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+k,i) ;
    end % for i=1:2*Na



end  % for k = 1:NumberGPSMeasurements

ys_kminus_0 = [PR_Vec_minus_0'; PRR_Vec_minus_0'];
ys_kminus_i = [PR_Vec_minus_i(:,:); PRR_Vec_minus_i(:,:)];

