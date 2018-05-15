%
% Perform a weighted least squares solution with integrity according to 
% "Weighted RAIM for Precision Approach" by Todd Walter and Per Enge ION
% GPS 95
%
% Written by Duncan Greer 17 May 2010
%



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

% % plot CNo v.s. bank angle
% for SV=1:size(CNo,2)
%     if SVWasTracked(SV) == 1
%         %figure(); grid on; hold on;
%         %plot(GPSTime-GPSTime_Start,CNo(:,SV));
%         %xlabel('GPS Seconds Since Start');
%         %ylabel('Carrier-Noise Ratio (dbHz)');
%         
%         % plot bank angle
%         %plot(span(:,4)-GPSTime_Start,span(:,11),'r');
%         
%         
%         %axis([0 1550 0 55]);
%         
%         xlabel = 'GPS Seconds Since Start';
%         ylabel1 = 'Carrier-Noise Ratio (dbHz)';
%         ylabel2 = 'Bank Angle (degrees)';
%         
%         %plotyy(GPSTime-GPSTime_Start,CNo(:,SV),span(:,4)-GPSTime_Start,span(:,11));
%         layerplot(GPSTime-GPSTime_Start,CNo(:,SV),span(:,4)-GPSTime_Start,span(:,11), {xlabel,ylabel1,ylabel2});
%         
%         title(sprintf('Carrier to Noise Ratio for PRN-%02d',SV));
%         
%     end
% end

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
   
%% RAIM solution
   % calculate 'G' which is the locally level (navigation frame) version of
   % the Least-squares geometry matrix, M
   G = lsqresult.LSQ_M * Tecef2ned2';
   
   % calculate the unweighted geometry matrices
   GG = inv(G' * G);
   
   A = GG * G';
   S = G * A;
   I = eye(size(S));
   S = I - S;

   % calculate the weighted geometry matrices
   GGw = inv(G' * W * G);
   
   sigma_V(Epoch_lo) = sqrt(GGw(3,3));
   sigma_H(Epoch_lo) = sqrt(GGw(1,1) + GGw(2,2));
   
   Aw = GGw * G' * W;
   Sw = G * Aw;
   I = eye(size(Sw));
   Sw = I - Sw;
   
   % calculate the unweighted slope
   for i=1:NumberGPSMeasurementsValid
    slope_V(i) = sqrt(A(3,i)^2) / sqrt(S(i,i));
    slope_H(i) = sqrt(A(1,i)^2 + A(2,i)^2) / sqrt(S(i,i));
   end
   [slopeMax_V(Epoch_lo),slopeMax_V_index] = max(slope_V);
   [slopeMax_H(Epoch_lo),slopeMax_H_index] = max(slope_H);
   
   slopeMax_SV_V(Epoch_lo) = SV_Vec(slopeMax_V_index);
   slopeMax_SV_H(Epoch_lo) = SV_Vec(slopeMax_H_index);
   
   % calculate the weighted slope
   for i=1:NumberGPSMeasurementsValid
    slope_Vw(i) = sqrt(Aw(3,i)^2) * (1/sqrt(W(i,i))) / sqrt(Sw(i,i));
    slope_Hw(i) = sqrt(Aw(1,i)^2 + Aw(2,i)^2) * (1/sqrt(W(i,i))) / sqrt(Sw(i,i));
   end
   [slopeMax_Vw(Epoch_lo),slopeMax_Vw_index] = max(slope_Vw);
   [slopeMax_Hw(Epoch_lo),slopeMax_Hw_index] = max(slope_Hw);
   
   slopeMax_SV_Vw(Epoch_lo) = SV_Vec(slopeMax_Vw_index);
   slopeMax_SV_Hw(Epoch_lo) = SV_Vec(slopeMax_Hw_index);
   
   
   % calculate hte parity matrix
   [Q,R] = qr(G);
   P = Q(:,5:NumberGPSMeasurementsValid)'; % lower n-4 rows of Q';
   
   % calculate the test statistic
   y = lsqresult.LSQ_ResidualVector(Epoch_lo,1:NumberGPSMeasurementsValid)';
   p = P*y;
   % the following values are equivelant for equal weighing (i.e. W = eye(N,N)  
   WSSE(Epoch_lo) = y' * W * S * y;
   sWSSE(Epoch_lo) = sqrt(WSSE(Epoch_lo));  %% this is the test statistic for weighted RAIM (T Walter, P Enge)
   pp(Epoch_lo) = norm(p); % or equivelantly, sqrt(p'*p);  % this is the test statistic for parity raim (RG Brown)
   
   
%% calculate protection levels using RG Brown pbias/rbias Method
   PBias(Epoch_lo) = raim.pbias(NumberGPSMeasurementsValid)*sigma_avg;
   rbias_H = PBias(Epoch_lo) / norm(P(:,slopeMax_H_index));
   rbias_V = PBias(Epoch_lo) / norm(P(:,slopeMax_V_index));
   
   e_H = zeros(NumberGPSMeasurementsValid,1);
   e_H(slopeMax_H_index) = rbias_H;
   
   E_H = A * e_H;
   HPL_rgb(Epoch_lo) = norm(E_H(1:2));
   
   e_V = zeros(NumberGPSMeasurementsValid,1);
   e_V(slopeMax_V_index) = rbias_V;
   E_V = A * e_V;
   VPL_rgb(Epoch_lo) = norm(E_V(3));
   
%% calculate protection level according to Van Dyke method - this is equivelant to the above
   HPL_vd(Epoch_lo) = PBias(Epoch_lo) * slopeMax_H(Epoch_lo);
   VPL_vd(Epoch_lo) = PBias(Epoch_lo) * slopeMax_V(Epoch_lo);
   
   
   
%% calculate threshold and protection level according to Todd Walter method
   % The Weighted raim threshold is the noramlalised chi-squared threshold as used by the partiy scheme (i.e. raim.Td)
   Threshold_tw(Epoch_lo) = raim.Td(NumberGPSMeasurementsValid);
   if sWSSE(Epoch_lo) > Threshold_tw(Epoch_lo)
      FaultDetected_tw(Epoch_lo) = 1;
      disp(sprintf('Fault detected (tw) at Epoch %d',Epoch_lo));
   else
       FaultDetected_tw(Epoch_lo) = 0;
   end
   
   HPL_tw(Epoch_lo) = slopeMax_Hw(Epoch_lo) * Threshold_tw(Epoch_lo) + sigma_H(Epoch_lo) * kPmd;
   VPL_tw(Epoch_lo) = slopeMax_Vw(Epoch_lo) * Threshold_tw(Epoch_lo) + sigma_V(Epoch_lo) * kPmd;
   
   
%% calculate threshold for detection - assume we have 6 or more satellites here
   Threshold_rgb(Epoch_lo) = raim.Td(NumberGPSMeasurementsValid) * sigma_avg;
   if (pp(Epoch_lo) > Threshold_rgb(Epoch_lo))
      % fault is declared
      FaultDetected_rgb(Epoch_lo) = 1;
      disp(sprintf('Fault detected (rgb) at Epoch %d',Epoch_lo));
   else
       FaultDetected_rgb(Epoch_lo) = 0;
   end
   
   
   
%% now; find that minimum slope max for a combination of 6 satellites from the full set
%    clear C Gc Ac Sc
%    C = nchoosek(SV_Vec,5); % returns a matrix whose rows are the available combinations
%    for i=1:size(C,1) % for each combination
%        
%        % build a G matrix (Gc) which includes just the satellites in this set
%        for j = 1:length(SV_Vec)
%            if find(C(i,:) == SV_Vec(j))
%                Gc(j,:) = G(j,:);
%            else
%                Gc(j,:) = zeros(1,4);
%            end
%        end
%        
%        % calculate the slope
%        Ac = inv(Gc'*Gc)*Gc';
%        Sc = eye(length(Gc),length(Gc)) - Gc * Ac;
%        
% 
%        for k=1:NumberGPSMeasurementsValid
%         slope_Vc(i,k) = sqrt(Ac(3,k)^2)  / sqrt(Sc(k,k));
%         slope_Hc(i,k) = sqrt(Ac(1,k)^2 + Ac(2,k)^2)  / sqrt(Sc(k,k));
%        end
%        
%        [slopeMax_Vc(i),slopeMax_Vc_index(i)] = max(slope_Vc(i,:));
%        [slopeMax_Hc(i),slopeMax_Hc_index(i)] = max(slope_Hc(i,:));
%        
%        
%    end
%    
%    [minSlopeMax_Hc(Epoch_lo),minSlopeMax_Hc_index(Epoch_lo)] = min(slopeMax_Hc);
%    [minSlopeMax_Vc(Epoch_lo),minSlopeMax_Vc_index(Epoch_lo)] = min(slopeMax_Vc);
% 
%    
   


%% calculate a HPL and VPL for multiple bias/faults following J.E. Angus "RAIM with Multiple Faults"
    MaxBiases = 1;
    Hh = [1 0 0 0;
          0 1 0 0];
    Hv = [0 0 1 0];
   
    clear C DQSh DQSheigs DQSheigMax DQSv DQSveigs DQSveigMax Dh Dv Qc HPL_mf_ VPL_mf_
    C = nchoosek(SV_Vec,MaxBiases);  % how many ways can we put MaxBiases into the number of satellites we have in view
    % loop over all possible combinations of the Q (bias connections) matrix
    for i=1:size(C,1)
        
        % construct a candidate Q matrix - it is n x r i.e. number SVs vs Number biases
        Qc = zeros(length(SV_Vec),MaxBiases);
        % now for each column, put a 1 in the row corresponding to the 'biased' satellite
        for j=1:MaxBiases
            
            k = find(SV_Vec == C(i,j));
            Qc(k,j) = 1;
            
        end % for j=1:MaxBiases
        
        Dh = W * G * inv(G'*W*G) * Hh' * Hh * inv(G'*W*G) * G' * W;
        Dv = W * G * inv(G'*W*G) * Hv' * Hv * inv(G'*W*G) * G' * W;
        
        
        DQSh = Qc' * Dh * Qc * inv(Qc' * W * S * Qc);
        DQSheigs = eigs(DQSh);
        DQSheigMax = max(DQSheigs);
        
        DQSv = Qc' * Dv * Qc * inv(Qc' * W * S * Qc);
        DQSveigs = eigs(DQSv);
        DQSveigMax = max(DQSveigs);

        HPL_mf_(i) = sqrt(DQSheigMax) * Threshold_tw(Epoch_lo) + sigma_H(Epoch_lo) * kPmd;
        VPL_mf_(i) = sqrt(DQSveigMax) * Threshold_tw(Epoch_lo) + sigma_V(Epoch_lo) * kPmd;
        

    end % for C =  1:size(C,1)

    HPL_mf(Epoch_lo) = max(HPL_mf_);
    VPL_mf(Epoch_lo) = max(VPL_mf_);
    

    
    
   if mod(Epoch_lo,100) == 0
      disp(sprintf('Completed epoch %d', Epoch_lo)); 
   end
    
   
   
   
    
    
    
end % for Epoch_lo

horizontal_error = sqrt(lsqerror.Pos_NED(:,1).^2 + lsqerror.Pos_NED(:,2).^2);
vertical_error = sqrt(lsqerror.Pos_NED(:,3).^2);

figure(); grid on; hold on;
plot(Time-Time(StartIndex),horizontal_error,'b')
plot(Time-Time(StartIndex),HPL_rgb,'g')
plot(Time-Time(StartIndex),FaultDetected_rgb*40,'r');
xlabel('Time (sec)');
ylabel('Horizontal Domain (meters)');
legend('Position Error','HPL (rgb method)','Fault Flag');
axis([0 1500 0 50]);

figure(); grid on; hold on;
plot(Time-Time(StartIndex),pp,'b');
plot(Time-Time(StartIndex),Threshold_rgb,'g');
plot(Time-Time(StartIndex),FaultDetected_rgb*40,'r');
xlabel('Time (sec)');
ylabel('Parity Domain');
legend('Test Statistic |p|','Threshold','Fault Flag');
axis([0 1500 0 50]);

figure(); grid on; hold on;
plot(Time-Time(StartIndex),vertical_error,'b')
plot(Time-Time(StartIndex),VPL_rgb,'g')
plot(Time-Time(StartIndex),FaultDetected_rgb*40,'r');
xlabel('Time (sec)');
ylabel('Vertical Domain (meters)');
legend('Vertical Error','VPL (rgb method)','Fault Flag');
axis([0 1500 0 50]);

% % plot slope max 
%     avi = avifile('slope.avi');
%     fig = figure(); grid on; 
%     set(fig,'DoubleBuffer','on');
%     for Epoch_lo = 30:700
% 
%         % plot current test statistic vs position error
%         %plot(pp(:),sqrt(lsqerror.Pos_NED(:,1).^2 + lsqerror.Pos_NED(:,2).^2),'r*');
%         plot(pp(Epoch_lo),sqrt(lsqerror.Pos_NED(Epoch_lo,1).^2 + lsqerror.Pos_NED(Epoch_lo,2).^2),'b*');
% 
%         % plot the slope
%         line([0 60],[0 60*slopeMax_H(Epoch_lo)],'Color','blue','LineStyle','-','LineWidth',2);
% 
%         % plot the threshold
%         line([Threshold_rgb(Epoch_lo) Threshold_rgb(Epoch_lo)],[0 60],'Color','red','LineStyle','--');
% 
%         % plot pbias
%         line([PBias(Epoch_lo) PBias(Epoch_lo)],[0 60],'Color','red','LineStyle','-');
% 
%         % plot HPL
%         line([0 60],[HPL_rgb(Epoch_lo) HPL_rgb(Epoch_lo)],'Color','blue','LineStyle','-','LineWidth',2);
% 
%         % plot HPL
%         line([0 60],[HPL_vd(Epoch_lo) HPL_vd(Epoch_lo)],'Color','green','LineStyle','-','LineWidth',2);
% 
%         axis([0 80 0 80]);
% 
%         frame = getframe(gca);
%         avi = addframe(avi,frame);
%         disp(sprintf('added frame %d',Epoch_lo));
%     end
%     avi = close(avi);

% figure(); grid on; hold on;
% plot(minSlopeMax_Hc,'b')
% plot(minSlopeMax_Vc,'b--')
% plot(slopeMax_H,'g')
% plot(slopeMax_V,'g--')
% xlabel('Time');
% ylabel('SlopeMax');

figure(); grid on; hold on;
plot(HPL_rgb,'b');
plot(VPL_rgb,'b--');
plot(HPL_tw,'r');
plot(VPL_tw,'r--');
xlabel('Time');
ylabel('Protection Level (meters');
legend('HPL (Parity RAIM)','VPL (Parity RAIM)','HPL (Weighted RAIM)','VPL (Weighted RAIM');

figure(); grid on; hold on;
plot(HPL_mf,'b');
plot(VPL_mf,'b--');
plot(HPL_tw,'r');
plot(VPL_tw,'r--');
xlabel('Time');
ylabel('Protection Level (meters');
legend('HPL (Multiple Fault)','VPL (Multiple Fault)','HPL (Weighted RAIM)','VPL (Weighted RAIM');


GARD_PlotStanford(horizontal_error,HPL_rgb',50,1,100,'Horizontal',1);
GARD_PlotStanford(vertical_error,VPL_rgb',50,1,100,'Vertical',1);

 
