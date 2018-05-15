UserPos = [LLH2ECEF(InitialPosition(1),InitialPosition(2),InitialPosition(3)), 0 ];
UserVel = [0 0 0 0];
NumberGPSEpochs = length(gps.RangeData);
NumberSVs = 40;
SVDontUse = zeros(1,NumberSVs);

%SVDontUse(24) = 1;

clear SVs_Used;

GPS_PR_UERE = 5.2;

SV_Ephemeris = GARD_GPSEphemStruct_to_Table(gps.GPSEphem);



    
PFalseAlarm = 1.06e-7; %
PMissedDetection = 0.0475;  % 95%

P_MD_Vert_Proportion = 0.95;
P_MD_H = (1-P_MD_Vert_Proportion) * PMissedDetection;
P_MD_V = P_MD_Vert_Proportion * PMissedDetection;

% [RAIM_a_H, RAIM_lambda_H] = GARD_CalculateThresholdPbias(1.06e-7,P_MD_H,[1:10]);
RAIM_a_H = [   28.2612   33.5059   37.5420   41.0297   44.1894   47.1251   49.8956   52.5377   55.0764   57.5294 ];
RAIM_lambda_H = [   66.0000   72.7000   77.5000   81.3000   84.7000   87.7000   90.4000   92.9000   95.3000   97.5000 ];

% [RAIM_a_V, RAIM_lambda_V] = GARD_CalculateThresholdPbias(1.06e-7,P_MD_V,[1:10]);
RAIM_a_V = [  28.2612   33.5059   37.5420   41.0297   44.1894   47.1251   49.8956   52.5377   55.0764   57.5294 ];
RAIM_lambda_V = [   49.2000   54.9000   59.0000   62.3000   65.2000   67.7000   70.1000   72.2000   74.3000   76.2000 ];




for Epoch_lo = 1:NumberGPSEpochs

    NumberGPSMeasurements = gps.RangeData(Epoch_lo).lNumberObservations;
    
     % format the input measurement vectors
    SVIndex = 0;
    clear SV_Vec SVPos SVVel PR_Vec PRR_Vec PR_Vec_raw SVAcc
    for Obs = 1:NumberGPSMeasurements
        
        SV = gps.RangeData(Epoch_lo).Obs(Obs).usPRN;
        
        if(SVDontUse(SV) == 0) 
            % add to PR vector
            SVIndex = SVIndex + 1;

            SV_Vec(SVIndex) = SV;
            PR_Vec(SVIndex) = gps.RangeData(Epoch_lo).Obs(Obs).dPseudorange;
            PRR_Vec(SVIndex) = gps.RangeData(Epoch_lo).Obs(Obs).fDoppler * -L1_Wavelength;

            
            
            [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV)] = ...
               GPSOrbitPropagator(gps.RangeData(Epoch_lo).GPSWeek, gps.RangeData(Epoch_lo).GPSSec - PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris, 7500);
           %  
           
           if ValidPosData(Epoch_lo,SV) == 0
               SVIndex = SVIndex - 1;
               continue;
           end
           
             [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
              SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(Epoch_lo,SV)] = ...
              GPSOrbitPropagatorVelocities(gps.RangeData(Epoch_lo).GPSWeek, gps.RangeData(Epoch_lo).GPSSec-PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris,7500);
          
          if ValidVelData(Epoch_lo,SV) == 0
              SVIndex = SVIndex - 1;
              continue;
          end
          
            SVPos(SVIndex,4) = SVPos(SVIndex,4) * Speedoflight;
            SVVel(SVIndex,4) = SVVel(SVIndex,4) * Speedoflight;
            SVAcc(SVIndex,4) = SVAcc(SVIndex,4) * Speedoflight;

            [SV_Azimuth(Epoch_lo,SV), SV_Elevation(Epoch_lo,SV)] = AzEl(UserPos(1:3), SVPos(SVIndex,1:3));
           
            if SV_Elevation(Epoch_lo,SV) < (ElevationMaskAngle * pi/180)
               SVIndex = SVIndex - 1;
               continue;
            end
            
            % calculate the iono delay correction - single frequency user
            % model from ICD 200
            
            ALPHA = [gps.IONUTCData(1).a0 gps.IONUTCData(1).a1 gps.IONUTCData(1).a2 gps.IONUTCData(1).a3];
            BETA = [gps.IONUTCData(1).b0 gps.IONUTCData(1).b1 gps.IONUTCData(1).b2 gps.IONUTCData(1).b3];

            ionodelay = ionomodel(gps.RangeData(Epoch_lo).GPSSec, UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
            PR_Vec(SVIndex) = PR_Vec(SVIndex) - ionodelay;
            
            
            if ~exist('x_save_llh','var')
                TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),InitialPosition(3));
            else    
                TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),x_save_llh(Epoch_lo-1,3));
            end
            PR_Vec(SVIndex) = PR_Vec(SVIndex) - TropoDelay(Epoch_lo,SV);
            
            % calculate hte earth rotation correction as per Kayton pg 228
            % eq 5.67

            
            PR_Vec_raw(SVIndex) = PR_Vec(SVIndex);  % save a raw (uncorrected copy) of the PR vector for use in the LSQ algorithm later.
            delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) *UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
            PR_Vec(SVIndex) = PR_Vec(SVIndex) + delta_pr_omegaedot + SVPos(SVIndex,4);
            PRR_Vec(SVIndex) = PRR_Vec(SVIndex) + SVVel(SVIndex,4);


        end
    end

    if exist('PR_Vec','var')
        NumberGPSMeasurements = SVIndex;%length(PR_Vec);
    else
        disp(sprintf('Epoch %d: No Observations'));
        continue;
    end

    NumberGPSMeasurementsAvailable(Epoch_lo) = NumberGPSMeasurements;
    SVs_Used(Epoch_lo,1:NumberGPSMeasurements) = SV_Vec(1:NumberGPSMeasurements);
    
    if NumberGPSMeasurements < 6
        disp(sprintf('Epoch %d: Insufficient measurements (%d) for FDE', Epoch_lo,NumberGPSMeasurements));
        continue;
    end
    
    if NumberGPSMeasurements > 8
        NumberGPSMeasurements = 8;
     end

    
   
   
   NumberGPSMeasurementsUsed(Epoch_lo) = NumberGPSMeasurements;
   
   
    
    %% calculate solution using least squares for comparisson
    [LSQ_Solution(Epoch_lo,:), LSQ_Variance(Epoch_lo,:), LSQ_NumIterations(Epoch_lo), ...
        LSQ_ResidualVector(Epoch_lo,1:NumberGPSMeasurements), LSQ_M, LSQ_Fail(Epoch_lo), ...
        LSQ_limit(Epoch_lo), LSQ_DOP(Epoch_lo,:)] = GARD_LSQ(UserPos,NumberGPSMeasurements, PR_Vec_raw(1:NumberGPSMeasurements),SVPos(1:NumberGPSMeasurements,:));
   
     [LSQVel_SolutionVec(Epoch_lo,:), LSQVel_VarSolutionVec(Epoch_lo,:), LSQVel_NumIterations(Epoch_lo,:),...
            LSQVel_ResidualVector(Epoch_lo,1:NumberGPSMeasurements), LSQVel_M, LSQVel_Fail(Epoch_lo,:), ...
            LSQVel_limit(Epoch_lo,:)] = GARD_LSQVel(UserPos,UserVel,NumberGPSMeasurements,PRR_Vec(1:NumberGPSMeasurements),SVPos(1:NumberGPSMeasurements,:), SVVel(1:NumberGPSMeasurements,:));

        
        UserPos = LSQ_Solution(Epoch_lo,1:4);
        UserVel = LSQVel_SolutionVec(Epoch_lo,1:4);
    
   [LSQ_Solution_LLH_noDGPS(Epoch_lo,1) LSQ_Solution_LLH_noDGPS(Epoch_lo,2) LSQ_Solution_LLH_noDGPS(Epoch_lo,3)] = ECEF2LLH(LSQ_Solution(Epoch_lo,1:3));
   Tecef2ned = T_ECEF2NED(LSQ_Solution_LLH_noDGPS(Epoch_lo,1),LSQ_Solution_LLH_noDGPS(Epoch_lo,2));
   Tecef2ned2 = Tecef2ned;
   Tecef2ned2(4,1:3) = [0 0 0];
   Tecef2ned2(1:4,4) = [0 0 0 1];
   
   LSQVel_Solution_NED_noDGPS(Epoch_lo,:) = Tecef2ned * LSQVel_SolutionVec(Epoch_lo,1:3)';
   
   
   % find span data for comparisson
   spanEpoch = find (round(gps.RangeData(Epoch_lo).GPSSec) < (span(:,4)),1);
   PosLLH_error_LSQ_noDGPS(Epoch_lo,:) = (LSQ_Solution_LLH_noDGPS(Epoch_lo,:) - span(spanEpoch,5:7) .* [pi/180 pi/180 1]) .* [RM RP*cos(-0.47) 1]; 
   
   
   if ~exist('a', 'var')
      load 'data/ThresholdPBias.mat' 
   end
   
   %% do RAIM Parity on LSQ solution - Note, you must run
   %% GARDSim_CalculateThresholdPbias before using this.
   [LSQ_BadGeometry(Epoch_lo), LSQ_RAIM_ALERT(Epoch_lo), LSQ_SLOPE_Max(Epoch_lo), LSQ_r(Epoch_lo), LSQ_Td(Epoch_lo), ...
       LSQ_RAIM_HPL(Epoch_lo),LSQ_RAIM_VPL(Epoch_lo), LSQ_FaultySatFDI(Epoch_lo)] = ...
       GARDSim_RAIMParity(RAIM_a_H, RAIM_lambda_H, NumberGPSMeasurements,PFalseAlarm,GPS_PR_UERE,556,LSQ_ResidualVector(Epoch_lo,1:NumberGPSMeasurements)',LSQ_M*Tecef2ned2');
       
   
   if mod(Epoch_lo,100) == 0
      disp(sprintf('Completed epoch %d', Epoch_lo)); 
   end
end