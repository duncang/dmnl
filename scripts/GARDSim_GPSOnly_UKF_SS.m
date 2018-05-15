
% GARDSim
% GPS Only Unscented Kalman Filter
% with Solution Separation Integrity




    
iono.ALPHA = [gps.IONUTCData(1).a0 gps.IONUTCData(1).a1 gps.IONUTCData(1).a2 gps.IONUTCData(1).a3];
iono.BETA = [gps.IONUTCData(1).b0 gps.IONUTCData(1).b1 gps.IONUTCData(1).b2 gps.IONUTCData(1).b3];
SV_Ephemeris = GARD_GPSEphemStruct_to_Table(gps.GPSEphem);


SVDontUse = zeros(1,32);

StartTime_GPS = 271848;
StopTime_GPS = 274718;

NumberGPSEpochs = size(gps.RangeData,2);

InitialPosition =  [ -0.481198373257095 ;        2.67055752492006  ;        58.7557216892019];
UserPos = [LLH2ECEF(InitialPosition(1),InitialPosition(2),InitialPosition(3)), 0 ];

for Epoch_lo = 1:NumberGPSEpochs
    
    [SVData(Epoch_lo) PRData(Epoch_lo)] = GARD_GetMeasurements(UserPos, gps.RangeData(Epoch_lo),SV_Ephemeris,SVDontUse,iono);    
    
    NumberGPSMeasurementsAvailable(Epoch_lo) = size(SVData(Epoch_lo).SV_Vec,2);
   
    
    % clear data from last pass
    clear z_Vec PR_Vec_minus R_k H_k
    
    % arrange matrices
    for k=1:NumberGPSMeasurementsAvailable(Epoch_lo)
        R_k(k,k) = PRData(Epoch_lo).PR_sigma(k);
        
        %Calculated slant ranges

        for m = 1:3
             ele(m) =  SVData(Epoch_lo).SV_Pos(k,m) - UserPos(m);
        end    

        r_VecCalc(k) =  norm(ele);   

        H_k(k,1) =  -ele(1)/r_VecCalc(k);
        H_k(k,2) =  -ele(2)/r_VecCalc(k);
        H_k(k,3) =  -ele(3)/r_VecCalc(k);
        H_k(k,16) = 1.0;   

        % find apriori estimate of pseudorange
        PR_Vec_minus(k) = r_VecCalc(k) + UserPos(4);  % note PR measurements have already been corrected for earth rotation and satellite clock

        
        
        
    end
    
    
    
    %% predict
    
    x_hat_minus
    
    
    
    %% correct
    
    
    
    % form measurements vector
    z_Vec = PRData(Epoch_lo).PR_Vec - PR_Vec_minus;
        
    
    
    
    
    
end





