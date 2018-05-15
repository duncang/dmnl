%
% Perform a multi-hypothesis solution separation least squares solution with integrity 
%
% Written by Duncan Greer 26 May 2010
%
% this is a sister script to GARDSim_WeightedRaim.m


%RinexObsFile = 'data/Flight_Data/log_22Jul09/22Jul09.09O';
%RinexNavFile = 'data/Flight_Data/log_22Jul09/22Jul09.09N';
InputLogFile = 'data/Flight_Data/log_22Jul09/log_20090722_gps.out';




% set number formatting for display
format long g;

% load GPS constants
GPSConstants;


% read rinex gps data
% [gps.GPSTime_Week, gps.GPSTime_Sec,gps.NumberRinexObsTypes, ...
%     gps.ValidDataRinexObs,gps.ApproxPos, ...
%     gps.Novatel_C1, gps.Novatel_L1, gps.Novatel_D1, gps.Novatel_S1, ...
%     gps.Novatel_P2, gps.Novatel_L2, gps.Novatel_D2, gps.Novatel_S2] = ReadRinexNovatel(RinexObsFile);
% 
% 
% 
% 
% [gps.SV_Ephemeris, gps.IONO_ALPHA, gps.IONO_BETA] = freadnav(RinexNavFile);


% gps = GARD_ReadNovatelLogData(InputLogFile);
% 
% gps.SV_Ephemeris = GARD_GPSEphemStruct_to_Table(gps.GPSEphem);
% 
% for i=1:length(gps.RangeData)/2
%     gps.RangeData(i) = gps.RangeData(i*2);
% end
% gps.RangeData = gps.RangeData(1:length(gps.RangeData)/2);
% %trim early bad data off gps
% gps.RangeData = gps.RangeData(47:1574);

load 'data/Flight_Data/log_22Jul09/gps-only.mat';
load 'data/Flight_Data/log_22Jul09/cptruth.mat';



% load truth data
%span = dlmread('data/Flight_Data/log_22Jul09/log_20090722_span.out');

load 'data/Flight_Data/log_22Jul09/log_20090722_span.mat';


SVWasTracked = zeros(32,1);

for i=1:length(gps.RangeData)
   for j=1:gps.RangeData(i).lNumberObservations
   prn = gps.RangeData(i).Obs(j).usPRN;
   GPSTime(i) = gps.RangeData(i).GPSSec;
   CNo(i,prn) = gps.RangeData(i).Obs(j).fCNo;
   LockTime(i,prn) = gps.RangeData(i).Obs(j).fLockTime;
   SVWasTracked(prn) = 1;
   end
end

GPSTime_Start = GPSTime(1);



% initialise a list of satellites to exclude by PRN
SVDontUse = zeros(32,1);
SVDontUse(22) = 1; % PRN 22 has a high URA (4)
%SVDontUse(15) = 1;
%SVDontUse(9) = 1;
%SVDontUse(3) = 1;
%SVDontUse(27) = 1;


[InitialPosition(1),InitialPosition(2),InitialPosition(3)] = ECEF2LLH([gps.BestXYZData(1).dPosX,gps.BestXYZData(1).dPosY,gps.BestXYZData(1).dPosZ]);
% initialise user pos and user vel
UserPos = [gps.BestXYZData(1).dPosX,gps.BestXYZData(1).dPosY,gps.BestXYZData(1).dPosZ, 0 ];
UserVel = [0 0 0 0];

ElevationMaskAngle = 7.5;
sigma_URA = 5;
sigma_multipath = 0.22;
sigma_UIRE = 4;
sigma_SNR = 0.22;
sigma_tropo = 0.7;
sigma_avg = sqrt(sigma_URA^2 + ...
            sigma_UIRE^2 + ...
            sigma_SNR^2 + ...
            (sigma_tropo^2)/(sin(45*pi/180)^2) + ...
            (sigma_multipath^2)/(tan(45*pi/180)^2));
    
PseudorangeBias = zeros(1528,1);

for i=500:length(PseudorangeBias)
   PseudorangeBias(i) = (i-500) * 0.5;  
end


% RAIM parameters
raim.Pfa = 1e-06;
raim.Pmd = 0.001;
raim.Td = [              0 , ...
                         0, ...
                         0, ...
                         0, ...
                         0, ...
          5.25653633994685, ...
          5.53760311399332, ...
          5.77726972136028, ...
          5.99066007536539, ...
          6.18534626967233, ...
          6.36565833518697, ...
          6.53458562207849, ...
           6.6940832148665, ...
          6.84567620826415];
raim.pbias = [           0, ...
                         0, ...
                         0, ...
                         0, ...
                         0, ...
          8.27127896115334, ...
          8.48245754381905, ...
          8.65601930661606, ...
          8.80620848618089, ...
          8.94006078563725, ...
          9.06153385066711, ...
          9.17335192577064, ...
          9.27725968357297, ...
          9.37462374678087];

kPmd = norminv(1-raim.Pmd/2);

% do normal LSQ solution
%for Epoch_lo = 1:length(gps.RangeData)
StartIndex = 30;
StopIndex = 1528;
for Epoch_lo = StartIndex:StopIndex
    truth_index = Epoch_lo - 29;
    
    
    NumberGPSMeasurementsAvailable = gps.RangeData(Epoch_lo).lNumberObservations;
    Time(Epoch_lo) = gps.RangeData(Epoch_lo).GPSSec;
    
    % format the input measurement vectors
    SVIndex = 0;
    clear SV_Vec SVPos SVVel PR_Vec PRR_Vec PR_Vec_raw SVAcc W
    for Obs = 1:NumberGPSMeasurementsAvailable
        
        SV = gps.RangeData(Epoch_lo).Obs(Obs).usPRN;
        
        if(SVDontUse(SV) == 0) 
            % add to PR vector
            SVIndex = SVIndex + 1;

            SV_Vec(SVIndex) = SV;
            PR_Vec(SVIndex) = gps.RangeData(Epoch_lo).Obs(Obs).dPseudorange;
            PRR_Vec(SVIndex) = gps.RangeData(Epoch_lo).Obs(Obs).fDoppler * -L1_Wavelength;

            %if SV == 14 
            %    PR_Vec(SVIndex) = PR_Vec(SVIndex) + PseudorangeBias(Epoch_lo);
            %end
            
            [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV), URA(Epoch_lo,SV)] = ...
               GPSOrbitPropagator(gps.RangeData(Epoch_lo).GPSWeek, gps.RangeData(Epoch_lo).GPSSec - PR_Vec(SVIndex)/Speedoflight, SV, gps.SV_Ephemeris, 7500);
           %  
           
           if ValidPosData(Epoch_lo,SV) == 0
               disp(sprintf('SV%02d Invalid Ephemeris at Epoch %d',SV,Epoch_lo));
               SVIndex = SVIndex - 1;
               continue;
           end
           
             [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
              SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(Epoch_lo,SV)] = ...
              GPSOrbitPropagatorVelocities(gps.RangeData(Epoch_lo).GPSWeek, gps.RangeData(Epoch_lo).GPSSec-PR_Vec(SVIndex)/Speedoflight, SV, gps.SV_Ephemeris,7500);
          
          if ValidVelData(Epoch_lo,SV) == 0
              SVIndex = SVIndex - 1;
              continue;
          end
          
            SVPos(SVIndex,4) = SVPos(SVIndex,4) * Speedoflight;
            SVVel(SVIndex,4) = SVVel(SVIndex,4) * Speedoflight;
            SVAcc(SVIndex,4) = SVAcc(SVIndex,4) * Speedoflight;

            % AzEl ordered by PRN
            [SV_Azimuth(Epoch_lo,SV), SV_Elevation(Epoch_lo,SV)] = AzEl(UserPos(1:3), SVPos(SVIndex,1:3));
            
            if SV_Elevation(Epoch_lo,SV) < (ElevationMaskAngle * pi/180)
               disp(sprintf('SV%02d below mask angle (%f < %f) - Removing from solution',SV,SV_Elevation(Epoch_lo,SV),ElevationMaskAngle));
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

             % form weighting matrix

            W(SVIndex,SVIndex) = 1/(URA(Epoch_lo,SV)^2 + sigma_UIRE^2 + sigma_SNR^2 + (sigma_tropo^2)/(sin(SV_Elevation(Epoch_lo,SV))^2) + (sigma_multipath^2)/(tan(SV_Elevation(Epoch_lo,SV))^2)); 


        end
    end

    if exist('PR_Vec','var')
        NumberGPSMeasurementsValid = SVIndex;
    else
        disp(sprintf('Epoch %d: No Observations'));
        continue;
    end

    lsqresult.NumberGPSMeasurementsAvailable(Epoch_lo) = NumberGPSMeasurementsAvailable;
    lsqresult.NumberGPSMeasurementsValid(Epoch_lo) = NumberGPSMeasurementsValid;
    
    lsqresult.SVs_Used(Epoch_lo,1:NumberGPSMeasurementsValid) = SV_Vec(1:NumberGPSMeasurementsValid);
    
    if NumberGPSMeasurementsValid < 6
        disp(sprintf('Epoch %d: Insufficient measurements (%d) for FDE', Epoch_lo,NumberGPSMeasurementsValid));
        continue;
    end
    

    % override weights matrix
    %W = eye(NumberGPSMeasurementsValid);
    
     %% calculate solution using normal least squares 
    [lsqresult.LSQ_Solution(Epoch_lo,:), lsqresult.LSQ_Variance(Epoch_lo,:), lsqresult.LSQ_NumIterations(Epoch_lo), ...
        lsqresult.LSQ_ResidualVector(Epoch_lo,1:NumberGPSMeasurementsValid), lsqresult.LSQ_M, lsqresult.LSQ_Fail(Epoch_lo), ...
        lsqresult.LSQ_limit(Epoch_lo), lsqresult.LSQ_DOP(Epoch_lo,:)] = GARD_LSQ(UserPos,NumberGPSMeasurementsValid, PR_Vec_raw(1:NumberGPSMeasurementsValid),SVPos(1:NumberGPSMeasurementsValid,:),W);
   
     [lsqresult.LSQVel_SolutionVec(Epoch_lo,:), lsqresult.LSQVel_VarSolutionVec(Epoch_lo,:), lsqresult.LSQVel_NumIterations(Epoch_lo,:),...
            lsqresult.LSQVel_ResidualVector(Epoch_lo,1:NumberGPSMeasurementsValid), lsqresult.LSQVel_M, lsqresult.LSQVel_Fail(Epoch_lo,:), ...
            lsqresult.LSQVel_limit(Epoch_lo,:)] = GARD_LSQVel(UserPos,UserVel,NumberGPSMeasurementsValid,PRR_Vec(1:NumberGPSMeasurementsValid),SVPos(1:NumberGPSMeasurementsValid,:), SVVel(1:NumberGPSMeasurementsValid,:));

        
        UserPos = lsqresult.LSQ_Solution(Epoch_lo,1:4);
        UserVel = lsqresult.LSQVel_SolutionVec(Epoch_lo,1:4);
    
   [lsqresult.LSQ_Solution_LLH(Epoch_lo,1) lsqresult.LSQ_Solution_LLH(Epoch_lo,2) lsqresult.LSQ_Solution_LLH(Epoch_lo,3)] = ECEF2LLH(lsqresult.LSQ_Solution(Epoch_lo,1:3));
   Tecef2ned = T_ECEF2NED(lsqresult.LSQ_Solution_LLH(Epoch_lo,1),lsqresult.LSQ_Solution_LLH(Epoch_lo,2));
   Tecef2ned2 = Tecef2ned;
   Tecef2ned2(4,1:3) = [0 0 0];
   Tecef2ned2(1:4,4) = [0 0 0 1];
   lsqresult.LSQVel_Solution_NED(Epoch_lo,:) = Tecef2ned * lsqresult.LSQVel_SolutionVec(Epoch_lo,1:3)';
   
   
   
   
   % calculate accuracy
   lsqerror.Pos_ECEF(Epoch_lo,:) = lsqresult.LSQ_Solution(Epoch_lo,1:3) -  [cptruth(Epoch_lo-29,3),cptruth(Epoch_lo-29,5),cptruth(Epoch_lo-29,7)];
   lsqerror.Pos_NED(Epoch_lo,:) = Tecef2ned * lsqerror.Pos_ECEF(Epoch_lo,:)';

   
   
%% now calculate sub solutions
   % exclude one measurement
   SubFilters_1 = nchoosek(SV_Vec,NumberGPSMeasurementsValid-1);
   
   NumberSubMeasurements = size(SubFilters_1,2);
   dW = diag(W);
   if ~exist('SubFilterIds','var')
    SubFilterIds(1).SVList = 0;
   end
   
   for SubFilter=1:size(SubFilters_1,1)
       % we need to generate a unique SubFilterId to reflect a particular combination of satellites
       FoundSubFilter = 0;
       for k=1:length(SubFilterIds)
           
           if length(SubFilterIds(k).SVList) ~=  length(SubFilters_1(SubFilter,:))
               % we know these won't match
               continue;
           end
           
           % we need a better detection here because this won't detect if the 
           % order has changed but the satelites are actually the same set
           %if SubFilterIds(k).SVList == SubFilters_1(SubFilter,:);
           if issameset(SubFilterIds(k).SVList,SubFilters_1(SubFilter,:))
               SubFilterId = k;
               FoundSubFilter = 1;
               break;
           end
           

       end
       
       if FoundSubFilter == 0
           % if we have gotten to here, then we didn't find an existing Id
           % need to allocate a new sub filter Id
           SubFilterId = length(SubFilterIds) + 1;
           SubFilterIds(SubFilterId).SVList = SubFilters_1(SubFilter,:);
           disp(sprintf('Allocating new SubFilterId (%d)',SubFilterId));
       end
       
       
       clear PR_Vec_Sub SV_Vec_Sub dW_Sub W_Sub
       % construct the PRVec and SVVecs
       for i=1:NumberSubMeasurements
          ii = find (SV_Vec == SubFilters_1(SubFilter,i));
          PR_Vec_Sub(i) = PR_Vec_raw(ii);
          SV_Pos_Sub(i,:) = SVPos(ii,:);
          dW_Sub(i) = dW(ii);
       end
       % de-diagonalise dW
       W_Sub = diag(dW_Sub);
   
       [lsqresult_1(SubFilterId).LSQ_Solution(Epoch_lo,:), ...
        lsqresult_1(SubFilterId).LSQ_Variance(Epoch_lo,:), ...
        lsqresult_1(SubFilterId).LSQ_NumIterations(Epoch_lo), ...
        lsqresult_1(SubFilterId).LSQ_ResidualVector(Epoch_lo,1:NumberSubMeasurements), ...
        lsqresult_1(SubFilterId).LSQ_M, lsqresult.LSQ_Fail(Epoch_lo), ...
        lsqresult_1(SubFilterId).LSQ_limit(Epoch_lo), ...
        lsqresult_1(SubFilterId).LSQ_DOP(Epoch_lo,:)] = GARD_LSQ(UserPos, ...
                                                               NumberSubMeasurements, ...
                                                               PR_Vec_Sub,...
                                                               SV_Pos_Sub,...
                                                               W_Sub);
                                                           
                                                           
        % calculate the separations - sub filter minus full filter
        SubFilterIds(SubFilterId).NumberSVs = length(SubFilterIds(SubFilterId).SVList);
        SubFilterIds(SubFilterId).ECEF(Epoch_lo,:) = lsqresult_1(SubFilterId).LSQ_Solution(Epoch_lo,1:3) - lsqresult.LSQ_Solution(Epoch_lo,1:3);
        SubFilterIds(SubFilterId).NED(Epoch_lo,:) =  (Tecef2ned * SubFilterIds(SubFilterId).ECEF(Epoch_lo,:)')';
        SubFilterIds(SubFilterId).Valid(Epoch_lo) = 1;
        
   end % for SubFilter_1... 
   
   % exclude 2 measurements
   %SubFilters_2 = nchoosek(SV_Vec,NumberGPSMeasurementsValid-2);
   
   
   % Fault Detection - Calculate Test Statistics
   
   
   
   % Calculate Protection Levels
   
    
   if mod(Epoch_lo,100) == 0
      disp(sprintf('Completed epoch %d', Epoch_lo)); 
   end
    
   
   
   
    
    
    
end % for Epoch_lo

horizontal_error = sqrt(lsqerror.Pos_NED(:,1).^2 + lsqerror.Pos_NED(:,2).^2);
vertical_error = sqrt(lsqerror.Pos_NED(:,3).^2);


figure(); grid on; hold on;
for i=2:length(SubFilterIds)
    plot(SubFilterIds(i).NED,'+')
end
