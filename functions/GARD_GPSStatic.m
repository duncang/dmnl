function [solution] = GARD_GPSStatic(datafile);
% function [solution] = GARD_GPSStatic(datafile);
% GPS Position Solution for static case using UKF, EKF and LSQ
% Written by Duncan Greer 2 Feb 2007
%
% $Id: GARD_GPSStatic.m 1850 2008-07-14 04:52:47Z greerd $
%
% INPUTS
% ======
% datafile - a matlab data file containing observation data, GPSTime, ION
% ALPHA and BETA, and SV_Ephemeris.
% 
% OUTPUTS
% =======
% solution - a data struct containing state solutions
%
%
% STATE VECTOR
% ============
% x1  = ECEF X
% x2  = ECEF Y
% x3  = ECEF Z
% x4  = ECEF X Velocity
% x5  = ECEF Y Velocity
% x6  = ECEF Z Velocity
% x7  = Receiver Clock Bias (Rb)
% x8  = Receiver Clock Drift Rate (Rf)
%
%




%% setup simulation

% set number formatting for display
format long g;

% load GPS constants
GPSConstants;

%% Load data - Note, this mat file is generated by the commented code
%% below.

disp(sprintf('Opening %s',datafile));
load(datafile);

% check for ApproxPos
if ~exist('ApproxPos','var')
    disp('Warning - Approx Pos Data Not Found - Using Defaults!');
    ApproxPos = [-5046719.69001578,2568403.35951166,-2925318.76002602];
end
    disp('Warning - meanECEF Data Not Found - Using Defaults: S-Block Roof');
if ~exist('meanECEF','var')
    meanECEF =  [-5046773.35774802  2568446.08440315 -2925289.01760974];
end

if ~exist('ALPHA','var')
    disp('Warning - ION ALPHA Data Not Found')
    ALPHA = [0 0 0 0];
end

if ~exist('BETA','var')
    disp('Warning - ION BETA Data Not Found')
    ALPHA = [0 0 0 0];
end

%% load observation and nav data
if ~exist('Novatel_C1','var')
    if exist('DATA_STRUCT','var')
        Novatel_C1 = DATA_STRUCT.C1(:,:);
    else
        disp('Error - No Observation Data Found!');
        return;
    end
end


if ~exist('Novatel_D1','var')
    Novatel_D1 = DATA_STRUCT.D1(:,:);
end
if ~exist('Novatel_L1','var')
    Novatel_L1 = DATA_STRUCT.L1(:,:);
end
if ~exist('Novatel_S1','var')
    Novatel_S1 = DATA_STRUCT.S1(:,:);
end



%% setup
NumberEpochsGPS = size(Novatel_C1,2);
gps_dt = 1;

NumberSVs = size(Novatel_C1,1);
SVDontUse = zeros(1,NumberSVs);




%% arrange pseudorange measurements

PR_Sim = Novatel_C1';
PRR_Sim = -Novatel_D1' * L1_Wavelength;
CP_Sim = Novatel_L1';

% position
[InitialPosition(1),InitialPosition(2),InitialPosition(3)] = ECEF2LLH(ApproxPos);
% velocity
InitialVelocity = [0 0 0]';  % NED velocities in m/s

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(2));
RP = PrimeRadius(InitialPosition(2));

% 3 position, 3 velocity and 2 clock
NumberStates = 8;

% process noise states - 8
ProcessNoiseStates = 8;  

% use a maximum of 6 measurements (PR+PRR)
MeasurementNoiseStates=12;


Na = NumberStates+ProcessNoiseStates+MeasurementNoiseStates;

% UKF scaling parameters
alpha = 0.1; % range: 1e-3 < alpha <= 1
beta = 2;  % 2 is optimal for gaussian priors
kapa = 0;  %

LAMBDA = alpha^2 * (Na+kapa) - Na;



TimeGPS = [0:gps_dt:NumberEpochsGPS-1];


%% initialise covariances
Px_kminus = eye(NumberStates,NumberStates);
Px_kminus(1:3,1:3) = eye(3,3)*100^2;
Px_kminus(4:6,4:6) = eye(3,3)*10^2;
Px_kminus(7,7) = 100^2;
Px_kminus(8,8) = 10^2;




GPS_PR_UERE = 7.5;
GPS_PRR_UERE = 0.2;

R = eye(MeasurementNoiseStates,MeasurementNoiseStates);
R(1:6,1:6) = eye(6,6)*GPS_PR_UERE^2;
R(7:12,7:12) = eye(6,6)*GPS_PRR_UERE^2;



x_hat_kminus = zeros(NumberStates,1);
%x_hat_kminus(1:3,1) = ApproxPos;
%x_hat_kminus(7,1) = 12;
GRASOn = 0;


x_hat_kminus(1:3,1) = [ -5046780.28152946          2568450.53282855         -2925294.35054563 ]';
x_hat_kminus(7,1) = 12.9024664188715;

SV_Azimuth = zeros(NumberEpochsGPS,NumberSVs);
SV_Elevation =  zeros(NumberEpochsGPS,NumberSVs);


SVDontUse(25) = 1;
SVDontUse(16) = 1;
SVDontUse(20) = 1;

URALimit = 1;

%% loop through and check for SV unhealthy (Ephemeris param 27)
for i=1:size(SV_Ephemeris,1)
   
   % check for unhealthy flag 
   if(SV_Ephemeris(i,27) ~= 0)
       disp(sprintf('SV %d Ephemeris Unhealthy',SV_Ephemeris(i,1)));
      SVDontUse(SV_Ephemeris(i,1)) = 1;
   end
   
   % check for large URA
   
   if(SV_Ephemeris(i,26) > URALimit)
       disp(sprintf('SV %d URA Exceeds Limit (%d)',SV_Ephemeris(i,1),URALimit));
      SVDontUse(SV_Ephemeris(i,1)) = 1;
   end
   
end



%% values for RAIM FDE

% PFalseAlarm = 2.22e-8;
% a =  [31.2919346425215; ...
%           35.2463471007126; ...
%           38.4953686364599; ...
%           41.4010185618922; ...
%            44.090770631714; ...
%           46.6276102566246; ...
%           49.0483514687134; ...
%           51.3768201044926; ...
%           53.6295337202138; ...
%           55.8185260702011; ...
%           57.9529000257331; ...
%           60.0397472340809; ...
%           62.0847253044164; ...
%            64.092436838505; ...
%            66.066688033301; ...
%            68.010670670063; ...
%           69.9270936151405; ...
%           71.8182800246577; ...
%           73.6862406187531; ...
%           75.5327299572836];
% 
% lambdatrue = [75.4140000000001; ...
%           80.2599999999999; ...
%           83.9499999999999; ...
%           87.0779999999999; ...
%           89.8469999999999; ...
%           92.3639999999999; ...
%           94.6879999999999; ...
%           96.8589999999999; ...
%                     98.906; ...
%                    100.848; ...
%                    102.698; ...
%                    104.468; ...
%                     106.17; ...
%                     107.81; ...
%                    109.396; ...
%                     110.93; ...
%                     112.42; ...
%                    113.865; ...
%                    115.275; ...
%                    116.647];
PFalseAlarm = 2e-5;
PMissedDetection = 0.047;

a =   [         18.189293484087;...
          21.6395565688235;...
          24.4623581112611;...
          26.9869539367368;...
          29.3272057089061;...
          31.5385558129051];

lambdatrue = [35.3000000000002;...
          38.9000000000003;...
          41.6000000000003;...
          43.8000000000004;...
          45.7000000000004;...
          47.5000000000004];
               
               
               
               
               
%% parameters for EKF

% power spectral densitites of the process noise - 
% Sp - position error
% Sf - clock bias error
% Sg - clock drift error

Sp = 10^2;  %% Static receiver
Sf = 2*2e-16*Speedoflight^2;
Sg = 8*pi^2*2e-17*Speedoflight^2;
% Sf = 10;
% Sg = 1;

% setup state transition - time invariant if dt is constant
EKF_phi = eye(NumberStates,NumberStates);
EKF_phi(1,5) = gps_dt;
EKF_phi(2,6) = gps_dt;
EKF_phi(3,7) = gps_dt;
EKF_phi(4,8) = gps_dt;


% setup Q and R matrices (process and noise covariance)
EKF_Q = zeros(NumberStates,NumberStates);

Sp_dt3 = Sp * (gps_dt ^ 3) / 3;
Sp_dt2 = Sp * (gps_dt ^ 2) / 2;
Sp_dt = Sp * gps_dt;
Sf_dt = Sf * gps_dt;
Sg_dt3 = Sg * (gps_dt ^ 3) / 3;
Sg_dt2 = Sg * (gps_dt ^ 2) / 2;
Sg_dt = Sg * gps_dt;

EKF_Q(1,1) = Sp_dt3;
EKF_Q(2,2) = Sp_dt3;
EKF_Q(3,3) = Sp_dt3;

EKF_Q(1,5) = Sp_dt2;
EKF_Q(2,6) = Sp_dt2;
EKF_Q(3,7) = Sp_dt2;
EKF_Q(4,8) = Sg_dt2;
EKF_Q(5,1) = Sp_dt2;
EKF_Q(6,2) = Sp_dt2;
EKF_Q(7,3) = Sp_dt2;
EKF_Q(8,4) = Sg_dt2;

EKF_Q(4,4) = Sf_dt + Sg_dt3;
EKF_Q(5,5) = Sp_dt;
EKF_Q(6,6) = Sp_dt;
EKF_Q(7,7) = Sp_dt;
EKF_Q(8,8) = Sg_dt;


EKF_R = R;

EKF_x_hat_kplus = x_hat_kminus;
EKF_Px_kplus = eye(NumberStates,NumberStates);
EKF_Px_kplus(1:3,1:3) = eye(3,3)*100^2;
EKF_Px_kplus(5:7,5:7) = eye(3,3)*10^2;
EKF_Px_kplus(4,4) = 100^2;
EKF_Px_kplus(8,8) = 10^2;
               



UKF_Q = zeros(ProcessNoiseStates,ProcessNoiseStates);
UKF_Q(1,1) = Sp_dt3;
UKF_Q(2,2) = Sp_dt3;
UKF_Q(3,3) = Sp_dt3;
UKF_Q(4,4) = Sp_dt;
UKF_Q(5,5) = Sp_dt;
UKF_Q(6,6) = Sp_dt;

UKF_Q(1,4) = Sp_dt2;
UKF_Q(2,5) = Sp_dt2;
UKF_Q(3,6) = Sp_dt2;

UKF_Q(4,1) = Sp_dt2;
UKF_Q(5,2) = Sp_dt2;
UKF_Q(6,3) = Sp_dt2;

UKF_Q(8,7) = Sg_dt2;
UKF_Q(7,8) = Sg_dt2;

UKF_Q(7,7) = Sf_dt + Sg_dt3;
UKF_Q(8,8) = Sg_dt;

%% start gps solutions

% initialise vectors
x_hat_save = zeros(NumberEpochsGPS,NumberStates);
P_save = zeros(NumberEpochsGPS,NumberStates);

EKF_x_hat_save = zeros(NumberEpochsGPS,NumberStates);
EKF_P_save = zeros(NumberEpochsGPS,NumberStates);

ECEFError_UKF = zeros(NumberEpochsGPS,3);
ECEFError_EKF = zeros(NumberEpochsGPS,3);
ECEFError_LSQ = zeros(NumberEpochsGPS,3);

NEDError_UKF = zeros(NumberEpochsGPS,3);
NEDError_EKF = zeros(NumberEpochsGPS,3);
NEDError_LSQ = zeros(NumberEpochsGPS,3);

x_save_vel = zeros(NumberEpochsGPS,3);
P_save_vel = zeros(NumberEpochsGPS,3);


for Epoch_lo = 1:NumberEpochsGPS

    % get the gps time - used for the iono correction
    %GPSTime  = 259199 + Epoch_lo;
    GPSTime = GPSTime_Sec(Epoch_lo);
    
    

    % perform state prediction
    % generate augmented state vector
    xa_hat_kminus = [x_hat_kminus;zeros(ProcessNoiseStates,1);zeros(MeasurementNoiseStates,1)];

    Pxa_kminus = [Px_kminus          zeros(NumberStates,NumberStates) zeros(NumberStates,MeasurementNoiseStates);
              zeros(ProcessNoiseStates,NumberStates) UKF_Q           zeros(ProcessNoiseStates,MeasurementNoiseStates);
              zeros(MeasurementNoiseStates,NumberStates) zeros(MeasurementNoiseStates,ProcessNoiseStates)  R];
          
          
          
          
          
    blah = sqrt(Na+LAMBDA) * chol(Pxa_kminus);
    for i=0:Na
        if i==0
            xs_0 = xa_hat_kminus;
            W_0_m = LAMBDA / (Na + LAMBDA);
            W_0_c = W_0_m + (1 - alpha^2 + beta);
        else
            xs_i(:,i) = xa_hat_kminus + blah(:,i);
            xs_i(:,i+Na) = xa_hat_kminus - blah(:,i);
            W_i_m = 1 / (2 * (Na + LAMBDA));
            W_i_c = W_i_m;
        end
    end

    
    
    
    
    % propagate system dynamics
    for i=0:2*Na
        if i==0
            %% zeroth point
            
            % position
            xs_0_k(1) = xs_0(1) + xs_0(4)*gps_dt;
            xs_0_k(2) = xs_0(2) + xs_0(5)*gps_dt;
            xs_0_k(3) = xs_0(3) + xs_0(6)*gps_dt;
            xs_0_k(7) = xs_0(7) + xs_0(8)*gps_dt;

            
            % velocity
            xs_0_k(4) = xs_0(4);
            xs_0_k(5) = xs_0(5);
            xs_0_k(6) = xs_0(6);
            xs_0_k(8) = xs_0(8);
%             xs_0_k(4) = xs_0(4) + xs_0(9)*gps_dt;
%             xs_0_k(5) = xs_0(5) + xs_0(10)*gps_dt;
%             xs_0_k(6) = xs_0(6) + xs_0(11)*gps_dt;
%             xs_0_k(8) = xs_0(8);
%             
%             % acceleration
%             xs_0_k(9) = xs_0(9);
%             xs_0_k(10) = xs_0(10);
%             xs_0_k(11) = xs_0(11);
%             
            
            
            %% add process noise
            V = xs_0(NumberStates+1:NumberStates+ProcessNoiseStates)';
            xs_0_k(1:NumberStates) = xs_0_k(1:NumberStates) + V;
            
            
        else
            
            % i-points
            
            % position
            xs_i_k(1,i) = xs_i(1,i) + xs_i(4,i)*gps_dt;
            xs_i_k(2,i) = xs_i(2,i) + xs_i(5,i)*gps_dt;
            xs_i_k(3,i) = xs_i(3,i) + xs_i(6,i)*gps_dt;
            xs_i_k(7,i) = xs_i(7,i) + xs_i(8,i)*gps_dt;

            
            % velocity
%             xs_i_k(4,i) = xs_i(4,i) + xs_i(9,i)*gps_dt;
%             xs_i_k(5,i) = xs_i(5,i) + xs_i(10,i)*gps_dt;
%             xs_i_k(6,i) = xs_i(6,i) + xs_i(11,i)*gps_dt;
%             xs_i_k(8,i) = xs_i(8,i);
            xs_i_k(4,i) = xs_i(4,i);
            xs_i_k(5,i) = xs_i(5,i);
            xs_i_k(6,i) = xs_i(6,i);
            xs_i_k(8,i) = xs_i(8,i);
            % acceleration
%             xs_i_k(9,i) = xs_i(9,i);
%             xs_i_k(10,i) = xs_i(10,i);
%             xs_i_k(11,i) = xs_i(11,i);
            
            
            
            %% add process noise
            V = xs_i(NumberStates+1:NumberStates+ProcessNoiseStates,i);
            xs_i_k(1:NumberStates,i) = xs_i_k(1:NumberStates,i) + V;
        end
    end
    
    
    
    
    
    
    
    
    W = eye(1,2*Na+1);
    W(1,1) = W_0_m;
    W(1,2:2*Na+1) = W_i_m;
    X = zeros(NumberStates,2*Na+1);
    X(:,1) = xs_0_k(1:NumberStates);
    X(:,2:2*Na+1) = xs_i_k(1:NumberStates,:);
    x_hat_kminus = (W*X')';

    nx = length(x_hat_kminus); nw = length(W);
    Px_kminus=((ones(nx,1)*W).*(X-x_hat_kminus*ones(1,nw)))*(X-x_hat_kminus*ones(1,nw))';
        
    % initial guess at user pos required to estimate earth rotation
    % correction
    UserPos(1:3) = [x_hat_kminus(1);x_hat_kminus(2);x_hat_kminus(3)];
 %   UserPos(4) = xs_0_k(7);
        
    % format the input measurement vectors
    SVIndex = 0;
    for SV=1:NumberSVs
        if((PR_Sim(Epoch_lo,SV) > 100) && SVDontUse(SV) == 0)  %%SV ~= ExcludeSVPRN
            % add to PR vector
            SVIndex = SVIndex + 1;

            SV_Vec(SVIndex) = SV;
            PR_Vec(SVIndex) = PR_Sim(Epoch_lo,SV);
            PRR_Vec(SVIndex) = PRR_Sim(Epoch_lo,SV);
            if(Epoch_lo == 1)
                PRR_Vec(SVIndex) = 0;
                %PR_csc(SVIndex) = PR_Vec(SVIndex);
            else
                if(CP_Sim(Epoch_lo-1,SV) == 0)
                    %PRR_Vec(SVIndex) = (CP_Sim(Epoch_lo+1,SV) - CP_Sim(Epoch_lo,SV));
                    %PR_csc(SVIndex) = PR_Vec(SVIndex);
                else
                    %PRR_Vec(SVIndex) = (CP_Sim(Epoch_lo,SV) - CP_Sim(Epoch_lo-1,SV));

                 %%% PERFORM CARRIER PHASE SMOOTHING of PR
                 %%% HERE


                   % alpha = gps_dt/100;
                    %P_proj = PR_prev + (L1_Wavelength/(2*pi))*(CP_Sim(Epoch_lo) - CP_Sim(Epoch_lo-1));
                    %PR_csc(SVIndex) = alpha * PR_Vec(SVIndex) + (1-alpha)*P_proj;                            

                end
            end
            % save previous value
            %PR_prev = PR_csc(SVIndex);

            % dont use the carrier smoothed code
            %PR_Vec = PR_csc;                  



%             SVPos(SVIndex,:) = [SV_X_Data(Epoch_lo,SV) SV_Y_Data(Epoch_lo,SV) SV_Z_Data(Epoch_lo,SV) SV_T_Data(Epoch_lo,SV)*c];
%             SVVel(SVIndex,:) = [SV_Xvel_Data(Epoch_lo,SV) SV_Yvel_Data(Epoch_lo,SV) SV_Zvel_Data(Epoch_lo,SV) SV_Tvel_Data(Epoch_lo,SV)*c];
%             SVAcc(SVIndex,:) = [SV_Xacc_Data(Epoch_lo,SV) SV_Yacc_Data(Epoch_lo,SV) SV_Zacc_Data(Epoch_lo,SV) SV_Tacc_Data(Epoch_lo,SV)*c];

            [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV)] = ...
               GPSOrbitPropagator(GPSTime_Week(Epoch_lo), GPSTime_Sec(Epoch_lo) - PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris, 7500);
           
             [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
              SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(Epoch_lo,SV)] = ...
              GPSOrbitPropagatorVelocities(GPSTime_Week(Epoch_lo),GPSTime_Sec(Epoch_lo)-PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris);
          
            SVPos(SVIndex,4) = SVPos(SVIndex,4) * Speedoflight;
            SVVel(SVIndex,4) = SVVel(SVIndex,4) * Speedoflight;
            SVAcc(SVIndex,4) = SVAcc(SVIndex,4) * Speedoflight;

            [SV_Azimuth(Epoch_lo,SV), SV_Elevation(Epoch_lo,SV)] = AzEl(UserPos(1:3), SVPos(SVIndex,1:3));
           
            % calculate the iono delay correction - single frequency user
            % model from ICD 200
            if(GRASOn == 0)
                ionodelay = ionomodel(GPSTime, UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
                PR_Vec(SVIndex) = PR_Vec(SVIndex) - ionodelay;
            end

            % calculate hte earth rotation correction as per Kayton pg 228
            % eq 5.67

            
            PR_Vec_raw(SVIndex) = PR_Vec(SVIndex);  % save a raw (uncorrected copy) of the PR vector for use in the LSQ algorithm later.
            delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) *UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
            PR_Vec(SVIndex) = PR_Vec(SVIndex) + delta_pr_omegaedot + SVPos(SVIndex,4);
            PRR_Vec(SVIndex) = PRR_Vec(SVIndex) + SVVel(SVIndex,4);


        end
    end

    NumberGPSMeasurements = length(PR_Vec);



    if NumberGPSMeasurements > 6
        NumberGPSMeasurements = 6;
    end
    
    if NumberGPSMeasurements < 6;
        disp(sprintf('Warning - Only %d measurements available on epoch %d of %d: skipping',NumberGPSMeasurements,Epoch_lo,NumberEpochsGPS));
        continue;
    end

    % formulate the measurment vector
    y_k = [PR_Vec(1:NumberGPSMeasurements)';PRR_Vec(1:NumberGPSMeasurements)'];
    

    for k = 1:NumberGPSMeasurements
        % find apriori estimate of pseudorange for each sigma point

        % zero-th sigma point

        
        % Get User position in ECEF
        UserPos(1:3) = [xs_0_k(1);xs_0_k(2);xs_0_k(3)];
        UserPos(4) = xs_0_k(7);
        UserVel(1:3) =  [xs_0_k(4);xs_0_k(5);xs_0_k(6)];
        UserVel(4) = xs_0_k(8);
        
        geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
        geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));

    %% calculate measurement predicition
        % pseudorange prediction
        PR_Vec_minus_0(k) = geo_range_to_sat + UserPos(4) + xs_0(NumberStates+ProcessNoiseStates+k);% + UserVel(4) * gps_dt;  % geometric range + c * delta_T
        %predicted relative velocity of sv and receiver
        Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
        PRR_Vec_minus_0(k) = Relative_Velocity(k) + UserVel(4) + xs_0(NumberStates+ProcessNoiseStates+MeasurementNoiseStates/2+k);

        for i=1:2*Na
            % Get User position in ECEF
            UserPos(1:3) = [xs_i_k(1,i);xs_i_k(2,i);xs_i_k(3,i)];
            UserPos(4) = xs_i_k(7,i);
            UserVel(1:3) =  [xs_i_k(4,i);xs_i_k(5,i);xs_i_k(6,i)];
            UserVel(4) = xs_i_k(8,i);
            
            
            geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
            geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));

            PR_Vec_minus_i(k,i) = geo_range_to_sat + UserPos(4) + xs_i(NumberStates+ProcessNoiseStates+k,i);% + UserVel(4) * gps_dt + xs_i(NumberStates+ProcessNoiseStates+k-1,i);  % geometric range + c * delta_T

            %predicted relative velocity of sv and receiver
            Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
            PRR_Vec_minus_i(k,i) = Relative_Velocity(k) + UserVel(4) + xs_i(NumberStates+ProcessNoiseStates+MeasurementNoiseStates/2+k,i);
        end % for i=1:2*Na

        % calculate EKF H-matrix
        UserPos(1:4) = EKF_x_hat_kplus(1:4);
        UserVel(1:4) = EKF_x_hat_kplus(5:8);
        for m = 1:3
             ele(m) =  SVPos(k,m) - UserPos(m);
        end    

        r_VecCalc(k) =  norm(ele);   

        EKF_H(k,1) =  -ele(1)/r_VecCalc(k);
        EKF_H(k,2) =  -ele(2)/r_VecCalc(k);
        EKF_H(k,3) =  -ele(3)/r_VecCalc(k);
        EKF_H(k,4) = 1;   
        EKF_H(k,5) = 0;
        EKF_H(k,6) = 0;
        EKF_H(k,7) = 0;
        EKF_H(k,8) = 0;

        EKF_H(k+NumberGPSMeasurements,1) = 0;
        EKF_H(k+NumberGPSMeasurements,2) = 0;
        EKF_H(k+NumberGPSMeasurements,3) = 0;
        EKF_H(k+NumberGPSMeasurements,4) = 0;   
        EKF_H(k+NumberGPSMeasurements,5) = -ele(1)/r_VecCalc(k);
        EKF_H(k+NumberGPSMeasurements,6) = -ele(2)/r_VecCalc(k);
        EKF_H(k+NumberGPSMeasurements,7) = -ele(3)/r_VecCalc(k);
        EKF_H(k+NumberGPSMeasurements,8) = 1;

        % find apriori estimate of pseudorange
        PR_Vec_minus(k) = r_VecCalc(k) + UserPos(4) + UserVel(4) * gps_dt;  % geometric range + c * delta_T
        %predicted relative velocity of sv and receiver
        r_VecCalcVel(k) = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
        Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);
        PRR_Vec_minus(k) = Relative_Velocity(k) + UserVel(4);
        

    end  % for k = 1:NumberGPSMeasurements

    ys_kminus_0 = [PR_Vec_minus_0'; PRR_Vec_minus_0'];
    ys_kminus_i = [PR_Vec_minus_i(:,:); PRR_Vec_minus_i(:,:)];
    
    % find the sum of sigma points for the measurement prediction
   

    W = ones(1,2*Na+1);
    W(1,1) = W_0_m;
    W(1,2:2*Na+1) = W_i_m;

    Y = zeros(MeasurementNoiseStates,2*Na+1);
    Y(:,1) = ys_kminus_0(1:MeasurementNoiseStates);
    Y(:,2:2*Na+1) = ys_kminus_i(1:MeasurementNoiseStates,:);
    y_hat_kminus = (W*Y')';

    nx = length(x_hat_kminus); ny = length(y_hat_kminus); nw = length(W);
    Py_kminus=((ones(ny,1)*W).*(Y-y_hat_kminus*ones(1,nw)))*(Y-y_hat_kminus*ones(1,nw))';
    Pxy_kminus =((ones(nx,1)*W).*(X-x_hat_kminus*ones(1,nw)))*(Y-y_hat_kminus*ones(1,nw))';

    %% Py_kminus is equivalent to HPxH'+R 
    %% Pxy_kminus PxH' 
    %% therefor Hequiv = (inv(Px_kminus)*Pxy_kminus)'
    % calculate 'H' by normalizing Pxy
    Hequiv = (inv(Px_kminus) * Pxy_kminus)';
    Hequiv2(1:NumberGPSMeasurements,1:3) = Hequiv(1:NumberGPSMeasurements,1:3);
    Hequiv2(1:NumberGPSMeasurements,4) = Hequiv(1:NumberGPSMeasurements,7);
    
    
    
    [lat,lon] = ECEF2LLH(ApproxPos);
    Tecef2ned2 = T_ECEF2NED(lon,lat);
    Tecef2ned2(4,1:3) = [0 0 0];
    Tecef2ned2(1:4,4) = [0 0 0 1];
    H_ltp = Hequiv2 * Tecef2ned2';
        
    % calculate covariance matrix (AA) and transform to local coords
    AA = inv(H_ltp' * H_ltp);
    
    % calculate DOPs
    var_x = AA(1,1);
    var_y = AA(2,2);
    var_z = AA(3,3);
    var_dt = AA(4,4);
    DOP_UKF(Epoch_lo,1) = sqrt(var_x + var_y + var_z + var_dt);
    DOP_UKF(Epoch_lo,2) = sqrt(var_x + var_y + var_z);
    DOP_UKF(Epoch_lo,3) = sqrt(var_x + var_y);
    DOP_UKF(Epoch_lo,4) = sqrt(var_z);
    DOP_UKF(Epoch_lo,5) = sqrt(var_dt);
    
    
    % calculate kalman gain
    K_k = Pxy_kminus * inv(Py_kminus);
    %K_k = (Pxy_kminus / Sy_kminus') / Sy_kminus;

    % apply correction
    z_k = y_k - y_hat_kminus;
    z_save(Epoch_lo,:) = z_k;
    

    
    x_hat_kplus = x_hat_kminus + K_k * (z_k);

    Px_kplus = Px_kminus - K_k * Py_kminus * K_k';


    
        
    %% calculate solution using least squares for comparisson
    [LSQ_Solution(Epoch_lo,:), LSQ_Variance(Epoch_lo,:), LSQ_NumIterations(Epoch_lo), ...
        LSQ_ResidualVector(Epoch_lo,:), LSQ_M, LSQ_Fail(Epoch_lo), ...
        LSQ_limit(Epoch_lo), LSQ_DOP(Epoch_lo,:)] = GARD_LSQ([ApproxPos 0],NumberGPSMeasurements, PR_Vec_raw(1:NumberGPSMeasurements),SVPos(1:NumberGPSMeasurements,:));
    
   %% do RAIM Parity on LSQ solution - Note, you must run
   %% GARDSim_CalculateThresholdPbias before using this.
   [BadGeometry(Epoch_lo), RAIM_ALERT(Epoch_lo), SLOPE_Max(Epoch_lo), r(Epoch_lo), Td(Epoch_lo), HPL(Epoch_lo),VPL(Epoch_lo), FaultySatFDI(Epoch_lo)] = ...
       GARDSim_RAIMParity(a, lambdatrue, NumberGPSMeasurements,PFalseAlarm,GPS_PR_UERE,556,LSQ_ResidualVector(Epoch_lo,:)',LSQ_M);
       
    
    SSE(Epoch_lo) = LSQ_ResidualVector(Epoch_lo,:)*LSQ_ResidualVector(Epoch_lo,:)';
    
    if(RAIM_ALERT(Epoch_lo))
       disp(sprintf('RAIM Error Detected at epoch %d',Epoch_lo)); 
    end
    
    
    %% calculate EKF solution for comparisson
    
    EKF_x_hat_kminus = EKF_x_hat_kplus;
    EKF_Px_kminus = EKF_Px_kplus;
    
    EKF_z = [PR_Vec(1:NumberGPSMeasurements)';PRR_Vec(1:NumberGPSMeasurements)'];
    EKF_z = ([PR_Vec(1:NumberGPSMeasurements) PRR_Vec(1:NumberGPSMeasurements)]' - [PR_Vec_minus(1:NumberGPSMeasurements) PRR_Vec_minus(1:NumberGPSMeasurements)]');
    
    EKF_z_save(Epoch_lo,:) = EKF_z;
    
    [EKF_x_hat_kplus, EKF_Px_kplus, EKF_v_out, EKF_s2_out] = ...
        GARD_EvaluateKF(gps_dt, EKF_x_hat_kminus, EKF_Px_kminus, ...
        EKF_phi, EKF_H, EKF_z, EKF_Q, EKF_R);
    
    
    EKF_x_hat_save(Epoch_lo,:) = EKF_x_hat_kplus;
    EKF_P_save(Epoch_lo,:) = diag(EKF_Px_kplus);
    

    
    %save results
    t_save(Epoch_lo) = Epoch_lo;
    x_hat_save(Epoch_lo,:) = x_hat_kplus;
    P_save(Epoch_lo,:) = diag(Px_kplus);


    % copy updated state to high-speed loop state for next round
    x_hat_kminus = x_hat_kplus;
    Px_kminus = Px_kplus;
    
    
        
    %% convert results to LLH
    [x_save_llh(Epoch_lo,1),x_save_llh(Epoch_lo,2),x_save_llh(Epoch_lo,3)] = ECEF2LLH(x_hat_save(Epoch_lo,1:3));
    Tecef2ned = T_ECEF2NED(x_save_llh(Epoch_lo,2),x_save_llh(Epoch_lo,1));
    

    P_save_NED(Epoch_lo,1:3) = abs(Tecef2ned*P_save(Epoch_lo,1:3)');
    S_save_NED(Epoch_lo,1:3) = sqrt(P_save_NED(Epoch_lo,1:3));

    EKF_P_save_NED(Epoch_lo,1:3) = abs(Tecef2ned*EKF_P_save(Epoch_lo,1:3)');
    EKF_S_save_NED(Epoch_lo,1:3) = sqrt(EKF_P_save_NED(Epoch_lo,1:3));
    
    EKF_P_save_NED(Epoch_lo,1:3) = abs(Tecef2ned*EKF_P_save(Epoch_lo,1:3)');
    EKF_S_save_NED(Epoch_lo,1:3) = sqrt(EKF_P_save_NED(Epoch_lo,1:3));
    
    %% calculate velocity in North-East-Down

    x_save_vel(Epoch_lo,:) = Tecef2ned*x_hat_save(Epoch_lo,4:6)';
    P_save_vel(Epoch_lo,:) = abs(Tecef2ned*P_save(Epoch_lo,4:6)');
    S_save_vel(Epoch_lo,:) = sqrt(P_save_vel(Epoch_lo,1:3));
    
    %% calculate position error
    ECEFError_UKF(Epoch_lo,:) = x_hat_save(Epoch_lo,1:3) - meanECEF;
    ECEFError_LSQ(Epoch_lo,:) = LSQ_Solution(Epoch_lo,1:3) - meanECEF;
    ECEFError_EKF(Epoch_lo,:) = EKF_x_hat_save(Epoch_lo,1:3) - meanECEF;
    
    %% calculate velocity in North-East-Down
    [lat,lon] = ECEF2LLH(meanECEF);
    Tecef2ned = T_ECEF2NED(lon,lat);
    
    NEDError_UKF(Epoch_lo,:) = Tecef2ned * ECEFError_UKF(Epoch_lo,:)';
    NEDError_LSQ(Epoch_lo,:) = Tecef2ned * ECEFError_LSQ(Epoch_lo,:)';
    NEDError_EKF(Epoch_lo,:) = Tecef2ned * ECEFError_EKF(Epoch_lo,:)';
    
   
    
    disp(sprintf('Completed Epoch %d of %d',Epoch_lo,NumberEpochsGPS));    
    
end



solution.x_hat_save = x_hat_save;
solution.P_save = P_save;
solution.EKF_x_hat_save = EKF_x_hat_save;
solution.EKF_P_save = EKF_P_save;
solution.ECEFError_UKF = ECEFError_UKF;
solution.ECEFError_EKF = ECEFError_EKF;
solution.ECEFError_LSQ = ECEFError_LSQ;
solution.NEDError_UKF = NEDError_UKF;
solution.NEDError_EKF = NEDError_EKF;
solution.NEDError_LSQ = NEDError_LSQ;
solution.x_save_vel = x_save_vel;
solution.P_save_vel = P_save_vel;


%% calculate the RMS error
%Height_RMS = sqrt(sum(Height_error.^2)/NumberEpochsINS)'
North_RMS_LSQ = sqrt(sum(NEDError_LSQ(:,1).^2)/NumberEpochsGPS);
North_RMS_EKF = sqrt(sum(NEDError_EKF(:,1).^2)/NumberEpochsGPS);
North_RMS_UKF = sqrt(sum(NEDError_UKF(:,1).^2)/NumberEpochsGPS);

East_RMS_LSQ = sqrt(sum(NEDError_LSQ(:,2).^2)/NumberEpochsGPS);
East_RMS_EKF = sqrt(sum(NEDError_EKF(:,2).^2)/NumberEpochsGPS);
East_RMS_UKF = sqrt(sum(NEDError_UKF(:,2).^2)/NumberEpochsGPS);

Down_RMS_LSQ = sqrt(sum(NEDError_LSQ(:,3).^2)/NumberEpochsGPS);
Down_RMS_EKF = sqrt(sum(NEDError_EKF(:,3).^2)/NumberEpochsGPS);
Down_RMS_UKF = sqrt(sum(NEDError_UKF(:,3).^2)/NumberEpochsGPS);

solution.RMSError_UKF = [North_RMS_UKF;East_RMS_UKF;Down_RMS_UKF];
solution.RMSError_EKF = [North_RMS_EKF;East_RMS_EKF;Down_RMS_EKF];
solution.RMSError_LSQ = [North_RMS_LSQ;East_RMS_LSQ;Down_RMS_LSQ];


North_RMS_LSQ_sub = sqrt(sum(NEDError_LSQ(30:length(NEDError_LSQ),1).^2)/NumberEpochsGPS);
North_RMS_EKF_sub = sqrt(sum(NEDError_EKF(30:length(NEDError_LSQ),1).^2)/NumberEpochsGPS);
North_RMS_UKF_sub = sqrt(sum(NEDError_UKF(30:length(NEDError_LSQ),1).^2)/NumberEpochsGPS);

East_RMS_LSQ_sub = sqrt(sum(NEDError_LSQ(30:length(NEDError_LSQ),2).^2)/NumberEpochsGPS);
East_RMS_EKF_sub = sqrt(sum(NEDError_EKF(30:length(NEDError_LSQ),2).^2)/NumberEpochsGPS);
East_RMS_UKF_sub = sqrt(sum(NEDError_UKF(30:length(NEDError_LSQ),2).^2)/NumberEpochsGPS);

Down_RMS_LSQ_sub = sqrt(sum(NEDError_LSQ(30:length(NEDError_LSQ),3).^2)/NumberEpochsGPS);
Down_RMS_EKF_sub = sqrt(sum(NEDError_EKF(30:length(NEDError_LSQ),3).^2)/NumberEpochsGPS);
Down_RMS_UKF_sub = sqrt(sum(NEDError_UKF(30:length(NEDError_LSQ),3).^2)/NumberEpochsGPS);