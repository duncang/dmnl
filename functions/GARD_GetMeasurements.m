function [SVData, PRData] = GARD_GetMeasurements(UserPos,gps_obs,SV_Ephemeris,SVDontUse,iono,ElevationMaskAngle)


% need GPS Constants
GPSConstants;

NumberGPSMeasurements = gps_obs.lNumberObservations;
SVIndex = 0;

% Elevation mask angle defaults to 7.5 deg if not supplied
if ~exist('ElevationMaskAngle','var')
    ElevationMaskAngle = 7.5*pi/180;
end

[UserPos_LLH(1) UserPos_LLH(2) UserPos_LLH(3) ] = ECEF2LLH(UserPos);

for Obs = 1:NumberGPSMeasurements
        
        SV = gps_obs.Obs(Obs).usPRN;
        
        if(SVDontUse(SV) == 0) 
            % add to PR vector
            SVIndex = SVIndex + 1;

            SV_Vec(SVIndex) = SV;
            PR_Vec(SVIndex) = gps_obs.Obs(Obs).dPseudorange;
            PRR_Vec(SVIndex) = gps_obs.Obs(Obs).fDoppler * -L1_Wavelength;

            
            
            [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(SV) URA(SV)] = ...
               GPSOrbitPropagator(gps_obs.GPSWeek, gps_obs.GPSSec - PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris, 7500);
           %  
           
           if ValidPosData(SV) == 0
               SVIndex = SVIndex - 1;
               continue;
           end
           
             [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
              SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(SV)] = ...
              GPSOrbitPropagatorVelocities(gps_obs.GPSWeek, gps_obs.GPSSec-PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris,7500);
          
          if ValidVelData(SV) == 0
              SVIndex = SVIndex - 1;
              continue;
          end
          
            SVPos(SVIndex,4) = SVPos(SVIndex,4) * Speedoflight;
            SVVel(SVIndex,4) = SVVel(SVIndex,4) * Speedoflight;
            SVAcc(SVIndex,4) = SVAcc(SVIndex,4) * Speedoflight;

            [SV_Azimuth(SV), SV_Elevation(SV)] = AzEl(UserPos(1:3), SVPos(SVIndex,1:3));
           
            
            % if below the mask angle, exclude this SV;
            if SV_Elevation(SV) < (ElevationMaskAngle)
               SVIndex = SVIndex - 1;
               continue;
            end
            
            % calculate the iono delay correction - single frequency user
            % model from ICD 200
            
            IonoDelay(SV) = ionomodel(gps_obs.GPSSec, UserPos(1:3), SVPos(SVIndex,1:3), iono.ALPHA, iono.BETA);
            PR_Vec(SVIndex) = PR_Vec(SVIndex) - IonoDelay(SV);
            
            
            % Calculate tropo delay per Blue book method
            TropoDelay(SV) =  GARD_TropoDelay(SV_Elevation(SV),UserPos_LLH(3));
            
            PR_Vec(SVIndex) = PR_Vec(SVIndex) - TropoDelay(SV);
            
            % calculate hte earth rotation correction as per Kayton pg 228
            % eq 5.67

            
            PR_Vec_raw(SVIndex) = PR_Vec(SVIndex);  % save a raw (uncorrected copy) of the PR vector for use in the LSQ algorithm later.
            delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) *UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
            PR_Vec(SVIndex) = PR_Vec(SVIndex) + delta_pr_omegaedot + SVPos(SVIndex,4);
            PRR_Vec(SVIndex) = PRR_Vec(SVIndex) + SVVel(SVIndex,4);

            
            % calculate the estimated PR sigma value from the URA and estimated user segment errors
            PR_sigma(SVIndex) = sqrt(URA(SV)^2 + 4.3^2);

        end
end
    

% output data

SVData.SV_Vec = SV_Vec;
SVData.SV_Pos = SVPos;
SVData.SV_Vel = SVVel;
SVData.SV_Acc = SVAcc;
SVData.ElevationMaskAngle = ElevationMaskAngle;
SVData.SV_Elevation = SV_Elevation;
SVData.SV_Azimuth = SV_Azimuth;
PRData.PR_Vec = PR_Vec;
PRData.PRR_Vec = PRR_Vec;
PRData.PR_Vec_raw = PR_Vec_raw;
PRData.TropoDelay = TropoDelay;
PRData.IonoDelay = IonoDelay;
PRData.PR_sigma = PR_sigma;


