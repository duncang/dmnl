

% obsfile = 'data/Flight_Data/00282031/00282031.09O';
% navfile = 'data/Flight_Data/00282031/00282031.09N';
% 
% BaseStation.BasePos = [-5041801.675, 2567494.879, -2934548.000];  %% from Auspos report
% 
% [BaseStation.GPSTime_Week, BaseStation.GPSTime_Sec,BaseStation.NumberRinexObsTypes,BaseStation.ValidDataRinexObs,...
%     BaseStation.ApproxPos, BaseStation.Novatel_C1,BaseStation.Novatel_L1, BaseStation.Novatel_D1, BaseStation.Novatel_S1,...
%     BaseStation.Novatel_P2, BaseStation.Novatel_L2, BaseStation.Novatel_D2, BaseStation.Novatel_S2] = ReadRinexNovatel(obsfile);
% 
% % update approx pos due to invalid approx pos in obs file
% [BaseStation.ApproxPos(1) BaseStation.ApproxPos(2) BaseStation.ApproxPos(3)] = ECEF2LLH(BasePos);
% 
% [BaseStation.NavData, BaseStation.IONO_ALPHA, BaseStation.IONO_BETA] = freadnav(navfile);
% 
% GPSConstants;


% pre-allocate

BaseStation.NumberSVs = size(BaseStation.Novatel_C1,1);
BaseStation.NumberEpochs = size(BaseStation.Novatel_C1,2);

BaseStation.GeometricRange = zeros(BaseStation.NumberSVs,BaseStation.NumberEpochs);
BaseStation.PRC = zeros(BaseStation.NumberSVs,BaseStation.NumberEpochs);
BaseStation.ValidPosData = zeros(BaseStation.NumberEpochs,BaseStation.NumberSVs);

BaseStation.NumberGPSMeasurements = zeros(BaseStation.NumberEpochs,1);
BaseStation.LSQ_Solution = zeros(BaseStation.NumberEpochs,4);
BaseStation.LSQ_Solution_LLH = zeros(BaseStation.NumberEpochs,3);
BaseStation.LSQ_Variance = zeros(BaseStation.NumberEpochs,4);
BaseStation.LSQ_NumIterations = zeros(BaseStation.NumberEpochs,1);
BaseStation.LSQ_ResidualVector = zeros(BaseStation.NumberEpochs,12);
BaseStation.LSQ_Fail = zeros(BaseStation.NumberEpochs,1);
BaseStation.LSQ_limit = zeros(BaseStation.NumberEpochs,1);
BaseStation.LSQ_DOP = zeros(BaseStation.NumberEpochs,5);

% initialise user pos
UserPos = [BaseStation.BasePos 0];
SVIndex = 0;


% for each epoch
for i = 1:BaseStation.NumberEpochs

    SVIndex = 0;
    % for each satellite PRN
    for j = 1:BaseStation.NumberSVs
        
        % if we have a measurement, calculate a PRC
        if BaseStation.Novatel_C1(j,i) ~= 0
            
            % add to PR vector
            SVIndex = SVIndex + 1;

            
            
            
            
            % get the satellite position
            [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), BaseStation.ValidPosData(i,j)] = ...
               GPSOrbitPropagator(BaseStation.GPSTime_Week(i), BaseStation.GPSTime_Sec(i) - BaseStation.Novatel_C1(j,i)/Speedoflight, j, BaseStation.NavData, 7500);
           
            delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) *UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
           
            if BaseStation.ValidPosData(i,j)         
                % get the psueodrange measurement
                BaseStation.GeometricRange(j,i) = sqrt((SVPos(SVIndex,1) - BaseStation.BasePos(1))^2 + ...
                                                       (SVPos(SVIndex,2) - BaseStation.BasePos(2))^2 + ...
                                                       (SVPos(SVIndex,3) - BaseStation.BasePos(3))^2);

                % calculate the PRC
                BaseStation.PRC(j,i) = BaseStation.GeometricRange(j,i) - BaseStation.Novatel_C1(j,i) - delta_pr_omegaedot;
            else
                disp(sprintf('Invalid Ephem for SV%d at Epoch %d',j,i));
            end
            
            
            %% apply correction
            PR_Vec(SVIndex) = BaseStation.Novatel_C1(j,i) + BaseStation.PRC(j,i);
            
            
            
        end % for each measurement available
    end % for each satellite j
    
    BaseStation.NumberGPSMeasurements(i) = SVIndex;
    
    
    %% now do a Least Squares solution on the PRs to check the performance
    [BaseStation.LSQ_Solution(i,:), BaseStation.LSQ_Variance(i,1:4), BaseStation.LSQ_NumIterations(i), ...
        BaseStation.LSQ_ResidualVector(i,1:BaseStation.NumberGPSMeasurements(i)), BaseStation.LSQ_M, BaseStation.LSQ_Fail(i), ...
        BaseStation.LSQ_limit(i), BaseStation.LSQ_DOP(i,:)] = GARD_LSQ(UserPos,BaseStation.NumberGPSMeasurements(i), ...
            PR_Vec(1:BaseStation.NumberGPSMeasurements(i)),SVPos(1:BaseStation.NumberGPSMeasurements(i),:),1);
   
    % update user pos
    UserPos = LSQ_Solution(Epoch_lo,1:4);
    
    % calculate NED errors
    Tecef2ned = T_ECEF2NED(BaseStation.LSQ_Solution_LLH(i,1),BaseStation.LSQ_Solution_LLH(i,2));
    BaseStation.NEDError(i,1:3) = Tecef2ned * (BaseStation.LSQ_Solution(i,1:3) - BaseStation.BasePos)';
    
    % calculate LLH
    [BaseStation.LSQ_Solution_LLH(i,1) BaseStation.LSQ_Solution_LLH(i,2) BaseStation.LSQ_Solution_LLH(i,3)] = ECEF2LLH(BaseStation.LSQ_Solution(i,1:3));
    
    if mod(i,100) == 0
       disp(sprintf('Completed Epoch %d',i)); 
    end
    
    
end % for each epoch i




