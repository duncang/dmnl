% Unscented Kalman Filter Implementation for integrated GPS-INS
%
% Implements an integrated GPS-Inertial position solution
%
% Written by Duncan Greer
%
% $Id: GARDSim_GPSINS_UKF_withCompass.m 2724 2009-07-26 12:48:03Z greerd $
%
% The GPS data is at 1 second epochs.  The INS data is calculated at 100Hz,
% or 0.01 second epochs.  The UKF will run at 1Hz, calculating the INS
% system errors
%
%
% x1  = Latitude position (rad)
% x2  = Longitude position (rad)
% x3  = Height position (m)
% x4  = North velocity (m/s)
% x5  = East velocity (m/s)
% x6  = Down velicty (m/s)
% x7  = q0
% x8  = q1
% x9  = q2
% x10 = q3
% x11 = x acc bias (m/s/s)
% x12 = y acc bias (m/s/s)
% x13 = z acc bias (m/s/s)
% x14 = x gyro bias (m/s/s)
% x15 = y gyro bias (m/s/s)
% x16 = z gyro bias (m/s/s)
% x17 = Clock Bias (m or sec?)
% x18 = Clock Drift (m/s or sec/sec?)
%
% There are two processes going on - a high speed loop evaluating the INS
% output, and a low-speed loop which integrates the GPS measurements.
%
%
%


% check if we are in the gardsim directory
if isempty(strfind(pwd,'gardsim'))
    warning('Are in in the GARDSim Directory?');
end

% set number formatting for display
format long g;

% load GPS constants
GPSConstants;

% decide if we want to make a movie
MakeMovie = 0;

ModuleName = '[GARDSim_GPSINS_UKF_withCompass] ';

% load data path
%DataPath = 'data/rnav_approach/';
%DataPath = 'data/DebugAeromodelpqrchanged11.1.08/';
%DataPath = 'data/DebugAeroModel7.1.08/';
%DataPath = 'data/Feb0108FlightNoWind/';

DataPath = 'data/Simulator_Data/Apr_08/data/rnav_approach_long/';


CalculateSubSolutions = 0;


%
gps_dt = 1;
ins_dt = 0.01;

GPSRATE = 1/gps_dt;
INSRATE = 1/ins_dt;




disp(strcat(ModuleName,'Loading data from: "',DataPath,'"'));



if CalculateSubSolutions
    disp('Calculating sub-solutions');
else
    disp('NOT calculating sub-solutions');
end


% save data path
dtg = clock;
ukf_results_filename = strcat(DataPath,'ukf_compass_results_');
ukf_results_filename = sprintf('%s%d%02d%02d_%02d%02d.mat',...
    ukf_results_filename,...
    dtg(1),...
    dtg(2),...
    dtg(3),...
    dtg(4),...
    dtg(5));
    

progress_filename = strcat(DataPath,'progress.txt');
% check if we are loading Spirent Sim data or GARDSim data
if isempty(strfind(DataPath,'Simulator_Data'))
    Sim_Data = 0;
else
    Sim_Data = 1;
end

if Sim_Data == 0
    % load the GPS and INS measurements
    if ~exist('PR_Sim','var')
        GPSMeasurementData = strcat(DataPath,'PR_Simulation.mat');
        load(GPSMeasurementData);
    end

    % load IMU measurements data
    UseNoisy = 0;

    GRASOn = 1;  % if GRAS is enabled then there is no iono bias

    if ~exist('sensors','var')

        if UseNoisy == 1
            load(strcat(DataPath,'sensors_noisy.mat'));
            sensors = sensors_noisy;
            clear 'sensors_noisy';
        else
            load(strcat(DataPath,'sensors_clean.mat'));
            sensors = sensors_clean;
            clear 'sensors_clean';
        end
    end

    %sensors(:,5:7) = sensors(:,5:7)*d2r - 0.0044;


    % load gravity data
    % if ~exist('gravity','var') 
    %     load(strcat(DataPath,'gravity.mat'));
    % end
    % iono model parameters for the above Nav file
    ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA           
    BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA  


    TruthOffset = 0;
else
    load(strcat(DataPath,'flexpak.mat'));
        GPSTime_Start = 259200;
    TruthOffset = round(GPSTime_Sec(1)-GPSTime_Start);


    PR_Sim = Novatel_C1';
    PRR_Sim = -Novatel_D1' * L1_Wavelength;
    CP_Sim = Novatel_L1';
    
    GRASOn = 0;
    
    % load IMU measurements data
    UseNoisy = 0;
    
    if ~exist('sensors','var')
        if UseNoisy == 1
            load(strcat(DataPath,'sensors_noisy.mat'));
            sensors = sensors_noisy;
            clear 'sensors_noisy';
        else
            load(strcat(DataPath,'sensors_clean.mat'));
            sensors = sensors_clean;
            clear 'sensors_clean';
        end
    end
end

% load truth data
if ~exist('pos_truth_llh','var')
    load(strcat(DataPath,'pos_truth_llh.mat'));
    load(strcat(DataPath,'pos_truth_ecef.mat'));
    
    load(strcat(DataPath,'vel_truth.mat'));
    load(strcat(DataPath,'att_truth.mat'));
end

disp(strcat(ModuleName,'Finished loading data'));


StartTime = TruthOffset+300;%1900;%  % start at 60 seconds into the data set
StopTime = TruthOffset+400;%size(PR_Sim,1);%


% initialise epochs
Epoch_lo = StartTime / gps_dt - 1;
Epoch_hi = StartTime / ins_dt - 1;

NumberStates = 18;

% process noise states - 6 for IMU sensors, 2 for gps clock
ProcessNoiseStates = 8;  

% measurement noise states - pseudorange and pseuorange rate for each
% satellite; plus compass
MeasurementNoiseStates=17;


% setup constants

WGS84Constants;
d2r = pi/180;
r2d = 180/pi;

ApproxPos = LLH2ECEF(-0.27,2.67,1000);

% position
InitialPosition = [pos_truth_llh(2,Epoch_hi) pos_truth_llh(3,Epoch_hi) pos_truth_llh(4,Epoch_hi)]';  % somewhere near brisbane in LLH, rads and meters
% velocity
InitialVelocity = [vel_truth(2,Epoch_hi) vel_truth(3,Epoch_hi) vel_truth(4,Epoch_hi)]';  % NED velocities in m/s
% attitude
InitialAttitude = [att_truth(2,Epoch_hi) att_truth(3,Epoch_hi) att_truth(4,Epoch_hi)]';  % roll, pitch, yaw in radians

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));

% gravity model - Kayton eqn 2.6
g = -9.79;  % m/s/s
%g = -0.01 * (978.049 * (1 + .00529 * sin(InitialPosition(1))^2));


% for the purpose of this simulation, the number of measurements is
% fixed at 16 = 2*8 - the best 8 measurements based on GDOP(??) should be chosen, plus there are PRR measurments for each SV.
NumberMeasurements = 17;
P0 = eye(NumberStates,NumberStates);

% initial position uncertainty
RMh = RM + InitialPosition(3);
RPh = RP + InitialPosition(3);
P0(1,1) = (5.0/RMh)^2;
P0(2,2) = (5.0/RPh*cos(InitialAttitude(1)))^2;

P0(3,3) = 5^2;

% initial velocity uncertainty
P0(4:6,4:6) = (0.5*eye(3,3))^2;

% initial attitude uncertainty
P0(7:10,7:10) = eye(4,4) * 0.0000001;

% initial sensor bias uncertainty
P0(11:13,11:13) = 0.000001*eye(3,3); % accelerometer
P0(14:16,14:16) = 0.0001*eye(3,3)*pi/180; % gyroscope

% initial clock uncertainty
P0(17,17) = 100;
P0(18,18) = 10;



% setup state transition - time invariant if dt is constant
% power spectral densitites of the process noise - 
x_gyro_beta = 1/300;
y_gyro_beta = 1/300;
z_gyro_beta = 1/300;
x_accel_beta = 1/300;
y_accel_beta = 1/300;
z_accel_beta = 1/300;
x_gyro_Q = 0.0552*pi/180; % PSD of gyro noise in radians
y_gyro_Q = 0.0552*pi/180; % PSD of gyro noise in radians
z_gyro_Q = 0.0552*pi/180; % PSD of gyro noise in radians
x_accel_Q = 0.0124; % PSD of accelerometer noise
y_accel_Q = 0.0124; % PSD of accelerometer noise
z_accel_Q = 0.0124; % PSD of accelerometer noise

% gps clock bias and frequency error PSD
Sc = 1*0.14;
Sf = 1*0.0359;

% initialise INS process noise
UKF_Q = zeros(ProcessNoiseStates,ProcessNoiseStates);

UKF_Q(1,1) = 2*x_accel_beta*x_accel_Q^2;
UKF_Q(2,2) = 2*y_accel_beta*y_accel_Q^2;
UKF_Q(3,3) = 2*z_accel_beta*z_accel_Q^2;
UKF_Q(4,4) = 2*x_gyro_beta*x_gyro_Q^2;
UKF_Q(5,5) = 2*y_gyro_beta*y_gyro_Q^2;
UKF_Q(6,6) = 2*z_gyro_beta*z_gyro_Q^2;

UKF_Q = UKF_Q * ins_dt;

UKF_Q(7,7) = (Sc)*ins_dt + Sf*(ins_dt^3)/3;
UKF_Q(8,8) = Sf*ins_dt;
UKF_Q(7,8) = Sf*(ins_dt^2)/2;
UKF_Q(8,7) = Sf*(ins_dt^2)/2;




% initialise measurement niose - pseudorange, then doppler 
sigma_pr_gnd = 2.0;
sigma_prr = 0.5;

if GRASOn == 0
    sigma_iono = 1.5;
else
    sigma_iono = 0.1;
end

sigma_pr = sqrt(sigma_pr_gnd^2 + sigma_iono^2);

UKF_R = zeros(NumberMeasurements,NumberMeasurements);
UKF_R(1:8,1:8) = sigma_pr^2*eye(8,8);
UKF_R(9:16,9:16) = sigma_prr^2*eye(8,8);

NumberSubMeasurements = 7;
UKF_R_Sub(1:NumberSubMeasurements,1:NumberSubMeasurements) = sigma_pr^2*eye(NumberSubMeasurements,NumberSubMeasurements);
UKF_R_Sub(NumberSubMeasurements+1:2*NumberSubMeasurements,NumberSubMeasurements+1:2*NumberSubMeasurements)...
    = sigma_prr^2*eye(NumberSubMeasurements,NumberSubMeasurements);

% simulate compass data from attitude truth
CompassVariance = (1.5*pi/180)^2;
CompassSensor = att_truth(4,:) + sqrt(CompassVariance)*randn(size(att_truth(4,:)));
% normalise to +- 180 deg
for i=1:length(CompassSensor)
    if CompassSensor(i) > pi
        CompassSensor(i) = CompassSensor(i) - 2*pi;
    elseif CompassSensor(i) < -pi
        CompassSensor(i) = CompassSensor(i) + 2*pi;
    end
end

UKF_R(17,17) = CompassVariance;
UKF_R_Sub(2*NumberSubMeasurements+1,2*NumberSubMeasurements+1) = CompassVariance;

Na = NumberStates+ProcessNoiseStates+MeasurementNoiseStates;

% UKF scaling parameters
alpha = 1; % range: 1e-3 < alpha <= 1
beta = 2;  % 2 is optimal for gaussian priors
kapa = 0;  %

LAMBDA = alpha^2 * (Na+kapa) - Na;

% nana is a variable which describes the convergence speed of the
% quaternion to unity
nana=0.5;


NumberEpochsGPS = min(StopTime,length(PR_Sim));

TimeINS = 0:ins_dt:NumberEpochsGPS;
TimeGPS = 0:gps_dt:NumberEpochsGPS-1;

NumberSVs = 32;
SVDontUse = zeros(1,NumberSVs);

%dt = gps_dt;

NumberEpochsINS = length(TimeINS);



disp(sprintf('[GARDSim_GPSINS_UKF] Beginning navigation loop for %d Epochs',NumberEpochsINS));



%% initialise storage vars to save time
Lat_error = zeros(1,NumberEpochsINS);
Long_error = zeros(1,NumberEpochsINS);
Height_error = zeros(1,NumberEpochsINS);

Clock_error = zeros(NumberEpochsGPS,1);
Vel_error = zeros(NumberEpochsINS,3);

UKF_HPL= zeros(NumberEpochsGPS,1);
UKF_VPL= zeros(NumberEpochsGPS,1);

% use these vectors to save the state vector at the end of each cycle
P_xn = zeros(1,NumberEpochsINS);
P_yn = zeros(1,NumberEpochsINS);
P_zn = zeros(1,NumberEpochsINS);
V_xn = zeros(1,NumberEpochsINS);
V_yn = zeros(1,NumberEpochsINS);
V_zn = zeros(1,NumberEpochsINS);

% save the covariance matrix at the end of each cycle

phi_q = zeros(NumberEpochsINS,1);
theta_q = zeros(NumberEpochsINS,1);
psi_q = zeros(NumberEpochsINS,1);

phi_err = zeros(NumberEpochsINS,1);
theta_err = zeros(NumberEpochsINS,1);
psi_err = zeros(NumberEpochsINS,1);


t_save = zeros(NumberEpochsINS,1);
UKF_x_hat_save = zeros(NumberEpochsINS,NumberStates);
UKF_P_save = zeros(NumberEpochsINS,NumberStates);



UKF_x_hat_save_Sub = zeros(NumberEpochsGPS,3,8);
UKF_P_save_Sub = zeros(NumberEpochsGPS,3,8);



UKF_z_save = zeros(NumberEpochsGPS,NumberMeasurements);
UKF_z_post_save = zeros(NumberEpochsGPS,NumberMeasurements);

% the size of the augmented state vector is given by A and includes the
% original states, plus the process and covariance noise states
A = NumberStates + NumberStates + NumberMeasurements;

LSQ_Solution = zeros(NumberEpochsGPS,4);
LSQ_Variance = zeros(NumberEpochsGPS,4);
LSQ_NumIterations = zeros(NumberEpochsGPS,1);
LSQ_ResidualVector = zeros(NumberEpochsGPS,8);
LSQ_Fail = zeros(NumberEpochsGPS,1);
LSQ_limit = zeros(NumberEpochsGPS,1);
LSQ_DOP = zeros(NumberEpochsGPS,5);
LSQ_RAIM_HPL = zeros(NumberEpochsGPS,1);
LSQ_RAIM_VPL = zeros(NumberEpochsGPS,1);

LSQ_FaultySatFDI = zeros(NumberEpochsGPS,1);
LSQ_SSE = zeros(NumberEpochsGPS,1);
LSQ_Td = zeros(NumberEpochsGPS,1);
LSQ_r = zeros(NumberEpochsGPS,1);
LSQ_SLOPE_Max = zeros(NumberEpochsGPS,1);
LSQ_RAIM_ALERT = zeros(NumberEpochsGPS,1);
LSQ_BadGeometry = zeros(NumberEpochsGPS,1);

LSQ_LLH_error = zeros(NumberEpochsGPS,3);
LSQ_ECEF_error = zeros(NumberEpochsGPS,3);
LSQVel_error = zeros(NumberEpochsGPS,3);
LSQ_Clock_error = zeros(NumberEpochsGPS,2);

UKF_HPL_sub = zeros(NumberEpochsGPS,8);
UKF_VPL_sub = zeros(NumberEpochsGPS,8);

UKF_lambda_ss_H = zeros(NumberEpochsGPS,8);
UKF_lambda_ss_V = zeros(NumberEpochsGPS,8);

UKF_TD_H = zeros(NumberEpochsGPS,8);
UKF_TD_V = zeros(NumberEpochsGPS,8);



ValidPosData = zeros(NumberEpochsGPS,NumberSVs);
ValidVelData = zeros(NumberEpochsGPS,NumberSVs);

SV_Azimuth = zeros(NumberEpochsGPS,1);
SV_Elevation = zeros(NumberEpochsGPS,1);

TropoDelay = zeros(NumberEpochsGPS,NumberSVs);
delta_pr_omegaedot = zeros(1,8);
ele = zeros(1,3);
    
PFalseAlarm = 1e-5; %
PMissedDetection = 0.048;  % 95%

RAIM_a =   [         18.189293484087;...
          21.6395565688235;...
          24.4623581112611;...
          26.9869539367368;...
          29.3272057089061;...
          31.5385558129051];

RAIM_lambdatrue = [35.3000000000002;...
          38.9000000000003;...
          41.6000000000003;...
          43.8000000000004;...
          45.7000000000004;...
          47.5000000000004];
               
      
P_MD_Vert_Proportion = 0.50;
P_MD_H = (1-P_MD_Vert_Proportion) * PMissedDetection;
P_MD_V = P_MD_Vert_Proportion * PMissedDetection;


    [lat,lon] = ECEF2LLH(ApproxPos);
    Tecef2ned2 = T_ECEF2NED(lat,lon);
    Tecef2ned2(4,1:3) = [0 0 0];
    Tecef2ned2(1:4,4) = [0 0 0 1];

    
%GyroBiasTruth=[0.01,-0.02,0.005];
GyroBiasTruth = zeros(3,NumberEpochsINS);
GyroNoiseTruth = randn(3,NumberEpochsINS)*0.00;

for i=1:NumberEpochsINS
%     GyroBiasTruth(1,i) = sin((1/10000)*i)*0.001;
%     GyroBiasTruth(2,i) = sin((1/15000)*i)*-0.002;
%     GyroBiasTruth(3,i) = sin((1/9000)*i)*0.00001;
end

%plot(GyroNoiseTruth(1,:)+GyroBiasTruth(1,:));



%AccelBiasTruth=[0.0,0.0,0.0];
AccelBiasTruth = zeros(3,NumberEpochsINS);
AccelNoiseTruth = randn(3,NumberEpochsINS)*0.0;

for i=1:NumberEpochsINS
    AccelBiasTruth(1,i) = sin((1/6000)*i)*0.0;
    AccelBiasTruth(2,i) = sin((1/12000)*i)*-0.00;
    AccelBiasTruth(3,i) = sin((1/13000)*i)*0.00000;
end

    
if MakeMovie == 1
    aviobj = avifile('solution_sep.avi','fps',5,'quality',100,'compression','None','videoname','Solution Separation Vectors');
end

% write to progress log file
dtg = clock;
fid = fopen(progress_filename,'a');
fprintf(fid,'\n==========\nStarting %s at %d%02d%02d - %02d:%02d\n==========\n',...
    ModuleName, ...
    dtg(1),...
    dtg(2),...
    dtg(3),...
    dtg(4),...
    dtg(5));
fclose(fid);

% initial value for number of sub filters
NumberSubFilters = 8;

            
            


            
%% high-speed loop - INS evaluation
for Epoch_hi = StartTime/ins_dt+1:StopTime/ins_dt+1
    
    
   if Epoch_hi == StartTime/ins_dt+1
       % first time step
  
       % initialise geodetic position
       %Pos_LLH(Epoch_hi,:) = InitialPosition;       
        Tecef2ned= T_ECEF2NED(InitialPosition(1),InitialPosition(2));
        Tned2ecef = Tecef2ned';

             
        % initialise the augmented state vector
        x_hat = zeros(NumberStates,1);
        x_hat(1)  = InitialPosition(1);
        x_hat(2)  = InitialPosition(2);
        x_hat(3)  = InitialPosition(3);
        x_hat(4)  = InitialVelocity(1);
        x_hat(5)  = InitialVelocity(2);
        x_hat(6)  = InitialVelocity(3);
%         x_hat(7)  = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
%         x_hat(8)  = sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) - cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
%         x_hat(9)  = cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
%         x_hat(10) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2) - sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2);
        x_hat(7:10) = EulerToQuat(InitialAttitude);

% initial bias are estimated as zero
        
   
        %x_hat_kplus = 
        UKF_x_hat_kminus =  x_hat;
        UKF_Px_kminus = P0;

        
        
        % initialise sub filters
        for SubSolution=1:NumberSubFilters
            UKF_x_hat_kminus_Sub(:,SubSolution) = UKF_x_hat_kminus;
            UKF_Px_kminus_Sub(:,:,SubSolution) = UKF_Px_kminus;

            [xs_Sub(:,:,SubSolution), Ws_Sub(:,SubSolution)] = ...
                GARD_GenerateSigmaPoints(UKF_x_hat_kminus,UKF_Px_kminus,UKF_Q,UKF_R_Sub,alpha,beta,kapa);
        end

        
        % initialise sigma points
        xs_0_k = zeros(Na,1);
        xs_i_k = zeros(Na,2*Na);

        % generate augmented state vector
        [xs,Ws] = GARD_GenerateSigmaPoints(UKF_x_hat_kminus,UKF_Px_kminus,UKF_Q,UKF_R,alpha,beta,kapa);
        xs_0 = xs(:,1);
        xs_i = xs(:,2:size(xs,2));
        W_0_m = Ws(1);
        W_0_c = Ws(2);
        W_i_m = Ws(3);
        W_i_c = Ws(4);
        
        
   else
        ins_dt = 0.01;%TimeINS(Epoch_hi) - TimeINS(Epoch_hi-1);
        %g = -0.01 * (978.049 * (1 + .00529 * sin(UKF_x_hat_kminus(1))^2));
        
        % get sensor measurements and correct for estimated biases
        A_xb_0 = sensors(2,Epoch_hi) + AccelBiasTruth(1,Epoch_hi) + AccelNoiseTruth(1,Epoch_hi);% + UKF_x_hat_kminus(11);
        A_yb_0 = sensors(3,Epoch_hi) + AccelBiasTruth(2,Epoch_hi) + AccelNoiseTruth(2,Epoch_hi);% + UKF_x_hat_kminus(12);
        A_zb_0 = sensors(4,Epoch_hi) + AccelBiasTruth(3,Epoch_hi) + AccelNoiseTruth(3,Epoch_hi);% + UKF_x_hat_kminus(13);
        omega_x_0 = sensors(5,Epoch_hi) + GyroBiasTruth(1,Epoch_hi) + GyroNoiseTruth(1,Epoch_hi);% + UKF_x_hat_kminus(14);
        omega_y_0 = sensors(6,Epoch_hi) + GyroBiasTruth(2,Epoch_hi) + GyroNoiseTruth(2,Epoch_hi);% + UKF_x_hat_kminus(15);
        omega_z_0 = sensors(7,Epoch_hi) + GyroBiasTruth(3,Epoch_hi) + GyroNoiseTruth(3,Epoch_hi);% + UKF_x_hat_kminus(16);
        
        Accel_in(:,Epoch_hi) = [A_xb_0;A_yb_0;A_zb_0];
        Gyro_in(:,Epoch_hi) = [omega_x_0;omega_y_0;omega_z_0];
        
        %% check if we've just come out of an update loop - if so we
        %% need to regenerate sigma points, if not, we propagate the
        %% existing sigma points.
        if mod(Epoch_hi-2,100) == 0
            Regen = 1;
%             [xs,Ws] = GARD_GenerateSigmaPoints(UKF_x_hat_kminus,UKF_Px_kminus,UKF_Q,UKF_R,alpha,beta,kapa);
%             xs_0 = xs(:,1);
%             xs_i = xs(:,2:size(xs,2));
%             W_0_m = Ws(1);
%             W_0_c = Ws(2);
%             W_i_m = Ws(3);
%             W_i_c = Ws(4);
        else
            Regen = 0;
        end
       

        
%         
%         % propogate sigma-points through system dynamics
%         for i=0:2*Na
%             if i==0  % the mean sigma point
%                 %% add process noise sigma points to sensor inputs
%                 Acc_in(1) = A_xb_0 - xs_0(11) + xs_0(NumberStates+1);
%                 Acc_in(2) = A_yb_0 - xs_0(12) + xs_0(NumberStates+2);
%                 Acc_in(3) = A_zb_0 - xs_0(13) + xs_0(NumberStates+3);
%                 Omega_in(1) = omega_x_0 - xs_0(14) + xs_0(NumberStates+4);
%                 Omega_in(2) = omega_y_0 - xs_0(15) + xs_0(NumberStates+5);
%                 Omega_in(3) = omega_z_0 - xs_0(16) + xs_0(NumberStates+6);
%                 
%                 g = GravityModel(xs_0(1:3));
%                 [xs_0_k,output] = GARD_INSMechanisation(xs_0,ins_dt,Acc_in,Omega_in,g);
%                 
%                 % update clock
%                 xs_0_k(17) = xs_0(17) + xs_0(18)*ins_dt + xs_0(NumberStates+7);
%                 xs_0_k(18) = xs_0(18) + xs_0(NumberStates+8);
%                                                      
% %                 % update random walk process for sensor biases
%                     
%                  % 1st order GM model
% %                 xs_0_k(11) = (1.0-x_accel_beta)*xs_0(11) + xs_0(19);
% %                 xs_0_k(12) = (1.0-y_accel_beta)*xs_0(12) + xs_0(20);
% %                 xs_0_k(13) = (1.0-z_accel_beta)*xs_0(13) + xs_0(21);
% %                 xs_0_k(14) = (1.0-x_gyro_beta)*xs_0(14) + xs_0(22);
% %                 xs_0_k(15) = (1.0-x_gyro_beta)*xs_0(15) + xs_0(23);
% %                 xs_0_k(16) = (1.0-x_gyro_beta)*xs_0(16) + xs_0(24);
% %                 
% %                 xs_0_k(11) = GaussMarkov_Process2(xs_0(11), x_accel_beta,xs_0(19)/ins_dt,ins_dt);
% %                 xs_0_k(12) = GaussMarkov_Process2(xs_0(12), y_accel_beta,xs_0(20)/ins_dt,ins_dt);
% %                 xs_0_k(13) = GaussMarkov_Process2(xs_0(13), z_accel_beta,xs_0(21)/ins_dt,ins_dt);
% %                 xs_0_k(14) = GaussMarkov_Process2(xs_0(14), x_gyro_beta,xs_0(22)/ins_dt,ins_dt);
% %                 xs_0_k(15) = GaussMarkov_Process2(xs_0(15), y_gyro_beta,xs_0(23)/ins_dt,ins_dt);
% %                 xs_0_k(16) = GaussMarkov_Process2(xs_0(16), z_gyro_beta,xs_0(24)/ins_dt,ins_dt);
%                 
%                  % random walk model
%                 xs_0_k(11) = xs_0(11) + xs_0(19)*ins_dt;
%                 xs_0_k(12) = xs_0(12) + xs_0(20)*ins_dt;
%                 xs_0_k(13) = xs_0(13) + xs_0(21)*ins_dt;
%                 xs_0_k(14) = xs_0(14) + xs_0(22)*ins_dt;
%                 xs_0_k(15) = xs_0(15) + xs_0(23)*ins_dt;
%                 xs_0_k(16) = xs_0(16) + xs_0(24)*ins_dt;
%                 % random constant model
% %                 xs_0_k(11) = xs_0(11);
% %                 xs_0_k(12) = xs_0(12);
% %                 xs_0_k(13) = xs_0(13);
% %                 xs_0_k(14) = xs_0(14);
% %                 xs_0_k(15) = xs_0(15);
% %                 xs_0_k(16) = xs_0(16);
% 
%                 % propogate process noise for sensor bias
%                  for m=19:24
%                      xs_0_k(m) = xs_0(m);
%                  end
%                 
%                 % propogate process noise for clock states
%                 xs_0_k(25) = xs_0(25);
%                 xs_0_k(26) = xs_0(26);
%                 
%                 % propogate measurement noise
%                 for m=27:Na
%                     xs_0_k(m) = xs_0(m);
%                 end
%                 
%                                                      
%             else
% 
%                 %% add process noise sigma points to sensor inputs
%                 Acc_in(1) = A_xb_0 - xs_i(11,i) + xs_i(NumberStates+1,i);
%                 Acc_in(2) = A_yb_0 - xs_i(12,i) + xs_i(NumberStates+2,i);
%                 Acc_in(3) = A_zb_0 - xs_i(13,i) + xs_i(NumberStates+3,i);
%                 Omega_in(1) = omega_x_0 - xs_i(14,i) + xs_i(NumberStates+4,i);
%                 Omega_in(2) = omega_y_0 - xs_i(15,i) + xs_i(NumberStates+5,i);
%                 Omega_in(3) = omega_z_0 - xs_i(16,i) + xs_i(NumberStates+6,i);
%                 
%                 g = GravityModel(xs_i(1:3,i));
%                 [xs_i_k(1:10,i),output] = GARD_INSMechanisation(xs_i(:,i),ins_dt,Acc_in,Omega_in,g);
%                 
%                 
% % 
%                 % update clock
%                 xs_i_k(17,i) = xs_i(17,i) + xs_i(18,i)*ins_dt + xs_i(NumberStates+7,i);
%                 xs_i_k(18,i) = xs_i(18,i) + xs_i(NumberStates+8,i);
%                                                      
%                 % update random walk process for sensor biases
% %                  xs_i_k(11,i) = (1.0-x_accel_beta)*xs_i(11,i) + xs_i(19,i);
% %                  xs_i_k(12,i) = (1.0-y_accel_beta)*xs_i(12,i) + xs_i(20,i);
% %                  xs_i_k(13,i) = (1.0-z_accel_beta)*xs_i(13,i) + xs_i(21,i);
% %                  xs_i_k(14,i) = (1.0-x_gyro_beta)*xs_i(14,i) + xs_i(22,i);
% %                  xs_i_k(15,i) = (1.0-y_gyro_beta)*xs_i(15,i) + xs_i(23,i);
% %                  xs_i_k(16,i) = (1.0-z_gyro_beta)*xs_i(16,i) + xs_i(24,i);
% % 
% %                 xs_i_k(11,i) = GaussMarkov_Process2(xs_i(11,i), x_accel_beta,xs_i(19,i)/ins_dt,ins_dt);
% %                 xs_i_k(12,i) = GaussMarkov_Process2(xs_i(12,i), y_accel_beta,xs_i(20,i)/ins_dt,ins_dt);
% %                 xs_i_k(13,i) = GaussMarkov_Process2(xs_i(13,i), z_accel_beta,xs_i(21,i)/ins_dt,ins_dt);
% %                 xs_i_k(14,i) = GaussMarkov_Process2(xs_i(14,i), x_gyro_beta,xs_i(22,i)/ins_dt,ins_dt);
% %                 xs_i_k(15,i) = GaussMarkov_Process2(xs_i(15,i), y_gyro_beta,xs_i(23,i)/ins_dt,ins_dt);
% %                 xs_i_k(16,i) = GaussMarkov_Process2(xs_i(16,i), z_gyro_beta,xs_i(24,i)/ins_dt,ins_dt);                 
%          
%                  xs_i_k(11,i) = xs_i(11,i) + xs_i(19,i)*ins_dt;
%                  xs_i_k(12,i) = xs_i(12,i) + xs_i(20,i)*ins_dt;
%                  xs_i_k(13,i) = xs_i(13,i) + xs_i(21,i)*ins_dt;
%                  xs_i_k(14,i) = xs_i(14,i) + xs_i(22,i)*ins_dt;
%                  xs_i_k(15,i) = xs_i(15,i) + xs_i(23,i)*ins_dt;
%                  xs_i_k(16,i) = xs_i(16,i) + xs_i(24,i)*ins_dt;
% 
% %                  xs_i_k(11,i) = xs_i(11,i);
% %                  xs_i_k(12,i) = xs_i(12,i);
% %                  xs_i_k(13,i) = xs_i(13,i);
% %                  xs_i_k(14,i) = xs_i(14,i);
% %                  xs_i_k(15,i) = xs_i(15,i);
% %                  xs_i_k(16,i) = xs_i(16,i);
% 
%                  % propogate process noise for sensor bias
%                  for m=19:24
%                      xs_i_k(m,i) = xs_i(m,i);
%                  end
%                  
%                  % propogate process noise for clock states
%                  xs_i_k(25,i) = xs_i(25,i);
%                  xs_i_k(26,i) = xs_i(26,i);
%                
%                 % propogate measurement noise
%                  for m=27:Na 
%                      xs_i_k(m,i) = xs_i(m,i);
%                  end
%             end
%         end
% 
%         
%         % update xs_i
%         xs_0 = xs_0_k;
%         xs_i = xs_i_k;
%         
%         % calculate the predicted 
%         % mean and state covariance
%          xs = [xs_0_k,xs_i_k];
%          [UKF_x_hat_kminus,UKF_Px_kminus] = GARD_GenerateMeanFromSigmaPoints(xs,Ws,NumberStates);

         UKF_Params.Q = UKF_Q;
         UKF_Params.R = UKF_R;
         UKF_Params.alpha = alpha;
         UKF_Params.beta = beta;
         UKF_Params.kapa = kapa;
         
         Sensor_Params.dt = ins_dt;
         Sensor_Params.Axb = A_xb_0;
         Sensor_Params.Ayb = A_yb_0;
         Sensor_Params.Azb = A_zb_0;
         Sensor_Params.p = omega_x_0;
         Sensor_Params.q = omega_y_0;
         Sensor_Params.r = omega_z_0;
         
         [UKF_x_hat_kminus,UKF_Px_kminus,xs,Ws]  = GARD_SigmaPointPropagation(...
                        UKF_x_hat_kminus,...
                        UKF_Px_kminus,...
                        UKF_Params,...
                        Sensor_Params,xs,Ws,Regen);
                    
          xs_0_k = xs(:,1);
          xs_i_k = xs(:,2:size(xs,2));          
                    
         


         
if CalculateSubSolutions
        UKF_Params.R = UKF_R_Sub;
        %% for each sub filter calculate the INS mechanisation
        for SubSolution = 1:NumberSubFilters
            
               [UKF_x_hat_kminus_Sub(:,SubSolution),UKF_Px_kminus_Sub(:,:,SubSolution),...
                   xs_Sub(:,:,SubSolution),Ws_Sub(:,SubSolution)]  = ...
                   GARD_SigmaPointPropagation(...
                        UKF_x_hat_kminus_Sub(:,SubSolution),...
                        UKF_Px_kminus_Sub(:,:,SubSolution),...
                        UKF_Params,...
                        Sensor_Params,xs_Sub(:,:,SubSolution),Ws_Sub(:,SubSolution),Regen);
            
        end
end % if calculatesubsolutions
        

         
        % save sigma points
        %xs_save_x_gyro_bias(Epoch_hi,:) = xs(14,:);
        % xs_save_y_gyro_bias(Epoch_hi,:) = xs(15,:);
        %xs_save_z_gyro_bias(Epoch_hi,:) = xs(16,:);
        
         % calculate euler angles
         euler = QuatToEuler(UKF_x_hat_kminus(7:10));
         phi_q(Epoch_hi) = euler(1);
         theta_q(Epoch_hi) = euler(2);
         psi_q(Epoch_hi) = euler(3);
         




   end
    
%% low-speed loop - kalman filter
    if(mod(Epoch_hi-1,100) == 0)

        Epoch_lo = Epoch_lo + 1;

        if Epoch_lo < TruthOffset
            warning('no gps data.. continuing');
            continue;
        end
        %disp(sprintf('Running Low Speed Loop %d',Epoch_lo));
        
          
        % Get User position in ECEF
        UserPos = LLH2ECEF(UKF_x_hat_kminus(1),UKF_x_hat_kminus(2),UKF_x_hat_kminus(3));
        UserPos(4) = UKF_x_hat_kminus(17);
        
        Tecef2ned= T_ECEF2NED(UKF_x_hat_kminus(1),UKF_x_hat_kminus(2));
        Tned2ecef = Tecef2ned';
        
        UserVel =  Tned2ecef * UKF_x_hat_kminus(4:6);
        UserVel(4) = UKF_x_hat_kminus(18);
        
        
        % format the input vectors
        SVIndex = 0;
        for SV=1:NumberSVs
            if((PR_Sim(Epoch_lo,SV) ~= 0) && SVDontUse(SV) == 0)  %%SV ~= ExcludeSVPRN
                % add to PR vector
                SVIndex = SVIndex + 1;

                SV_Vec(SVIndex) = SV;
                
if Sim_Data == 1
                PR_Vec(SVIndex) = PR_Sim(Epoch_lo-TruthOffset+1,SV);
else
                PR_Vec(SVIndex) = PR_Sim(Epoch_lo,SV);    
end

if Sim_Data == 1
                % random hack
                if(Epoch_lo == 1)
                    PRR_Vec(SVIndex) = PRR_Sim(2,SV);
                else
                    PRR_Vec(SVIndex) = PRR_Sim(Epoch_lo-TruthOffset+1,SV);
                end
                    
                if(PRR_Vec(SVIndex) == 0)
                    PRR_Vec(SVIndex) = PRR_Sim(Epoch_lo+1-TruthOffset+1,SV);
                end
else
                    PRR_Vec(SVIndex) = PRR_Sim(Epoch_lo,SV);    
end

if Sim_Data == 1
                [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV)] = ...
                    GPSOrbitPropagator(GPSTime_Week(Epoch_lo-TruthOffset+1), GPSTime_Sec(Epoch_lo-TruthOffset+1)-PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris, 7500);
else
                [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV)] = ...
                    GPSOrbitPropagator(GPSTime_Week(Epoch_lo), GPSTime_Sec(Epoch_lo), SV, SV_Ephemeris, 7500);    
end
                if ValidPosData(Epoch_lo,SV) ~= 1
                    SVDontUse(SV) = 1;
                    SVIndex = SVIndex - 1;
                    disp(sprintf('Warning: SV%d pos not valid',SV));
                    continue;
                end
                
if Sim_Data == 1                
                [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
                    SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(Epoch_lo,SV)] = ...
                    GPSOrbitPropagatorVelocities(GPSTime_Week(Epoch_lo-TruthOffset+1),GPSTime_Sec(Epoch_lo-TruthOffset+1)-PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris);
else
                [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
                    SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(Epoch_lo,SV)] = ...
                    GPSOrbitPropagatorVelocities(GPSTime_Week(Epoch_lo),GPSTime_Sec(Epoch_lo), SV, SV_Ephemeris);    
end
                
                if ValidVelData(Epoch_lo,SV) ~= 1
                    SVDontUse(SV) = 1;
                    SVIndex = SVIndex - 1;
                    disp('Warning: SV%d vel not valid',SV);
                    continue;
                end
                
if Sim_Data == 1                
                 SVPos(SVIndex,4) = SVPos(SVIndex,4) * Speedoflight;
                 SVVel(SVIndex,4) = SVVel(SVIndex,4) * Speedoflight;
                 SVAcc(SVIndex,4) = SVAcc(SVIndex,4) * Speedoflight;
else
    
end
                [SV_Azimuth(Epoch_lo,SV), SV_Elevation(Epoch_lo,SV)] = AzEl(UserPos(1:3), SVPos(SVIndex,1:3));

                % calculate the iono delay correction - single frequency user
                % model from ICD 200
                if GRASOn == 0
                    ionodelay = ionomodel(GPSTime_Sec(Epoch_lo-TruthOffset+1), UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
                    PR_Vec(SVIndex) = PR_Vec(SVIndex) - ionodelay;
                    sigma_iono = 1.5;
                else
                    PR_Vec(SVIndex) = PR_Vec(SVIndex);
                    sigma_iono = 0.1;
                end
                
                if Epoch_lo == StartTime
                    TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),InitialPosition(3));
                else    
                    TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),UKF_x_hat_save(Epoch_hi-1,3));
                end
                
                PR_Vec(SVIndex) = PR_Vec(SVIndex) - TropoDelay(Epoch_lo,SV);
                PR_Vec_raw(SVIndex) = PR_Vec(SVIndex);  % save a raw (uncorrected copy) of the PR vector for use in the LSQ algorithm later.
                

            end
        end

        NumberGPSMeasurements = length(PR_Vec);

          %if(Epoch_lo > 120)
          %   PR_Vec(1) =  PR_Vec(1) + (Epoch_lo-StartTime); % 1 m/s ramp fault
          %end
        
        if NumberGPSMeasurements < 8
            continue;
        end
        if NumberGPSMeasurements > 8
            NumberGPSMeasurements = 8;
        end

        % formulate the measurment vector
        UKF_y_k = [PR_Vec(1:8)';PRR_Vec(1:8)';CompassSensor(Epoch_hi)];
        
        PR_Vec_minus_0 = zeros(NumberGPSMeasurements,1);
        PRR_Vec_minus_0 = zeros(NumberGPSMeasurements,1);
        
        PR_Vec_minus_i = zeros(NumberGPSMeasurements,2*Na);
        PRR_Vec_minus_i = zeros(NumberGPSMeasurements,2*Na);
        
        Relative_Velocity = zeros(NumberGPSMeasurements,1);
        
        
        % get the a-proiri measuremnt prediction (y_minus)
        for k = 1:NumberGPSMeasurements
   
            [PR_Vec_minus_0(k),PRR_Vec_minus_0(k)] = GARD_GPSMeasurementPrediction(...
                   [xs_0_k(1:3); xs_0_k(17)],...
                   [xs_0_k(4:6); xs_0_k(18)],...
                   SVPos(k,:),...
                   SVVel(k,:));
                
            PR_Vec_minus_0(k) = PR_Vec_minus_0(k) + xs_0(NumberStates+ProcessNoiseStates+k);
            PRR_Vec_minus_0(k) = PRR_Vec_minus_0(k) + xs_0(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+k);
            
            
            
            for i=1:2*Na

            
            [PR_Vec_minus_i(k,i),PRR_Vec_minus_i(k,i)] = GARD_GPSMeasurementPrediction(...
                   [xs_i_k(1:3,i); xs_i_k(17,i)],...
                   [xs_i_k(4:6,i); xs_i_k(18,i)],...
                   SVPos(k,:),...
                   SVVel(k,:));
                
            PR_Vec_minus_i(k,i) = PR_Vec_minus_i(k,i) + xs_i(NumberStates+ProcessNoiseStates+k,i);
            PRR_Vec_minus_i(k,i) = PRR_Vec_minus_i(k,i) + xs_i(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+k,i);
            
            
            end % for i=1:2*Na



        end  % for k = 1:NumberGPSMeasurements
        
        %% compass measurement
        CompassHeading_minus_0 = GARD_HeadingFromQuaternion(xs_0_k(7:10)) + xs_0(NumberStates+ProcessNoiseStates+2*NumberGPSMeasurements+1);
        for i=1:2*Na
            CompassHeading_minus_i(i) = GARD_HeadingFromQuaternion(xs_i_k(7:10,i)) + xs_i(NumberStates+ProcessNoiseStates+2*NumberGPSMeasurements+1,i);
        end
        
        ys_kminus_0 = [PR_Vec_minus_0; PRR_Vec_minus_0;CompassHeading_minus_0];
        ys_kminus_i = [PR_Vec_minus_i(:,:); PRR_Vec_minus_i(:,:);CompassHeading_minus_i];


        % calculate measurement covariance
         ys = [ys_kminus_0,ys_kminus_i];
         [UKF_y_hat_kminus,UKF_Py_kminus,UKF_Pxy_kminus] = GARD_GenerateMeasurementMeanFromSigmaPoints(xs,ys,Ws,NumberStates,NumberMeasurements);


%         Hequiv = (inv(UKF_Px_kminus) * UKF_Pxy_kminus)';
%         
% %        [ObsRank(Epoch_lo),ObsMatrix] = GARD_CalculateObservability(F,Hequiv) 
%         
%         H_ltp(:,1:3) = Hequiv(9:16,4:6);
%         H_ltp(:,4) = Hequiv(9:16,18);
%         
%         H_ltp2(:,1:3) = Hequiv(1:8,1:3);
%         H_ltp2(:,1) = H_ltp2(:,1)/RMh;
%         H_ltp2(:,2) = H_ltp2(:,2)/(RPh*cos(InitialPosition(1)));
%         H_ltp2(:,3) = H_ltp2(:,3)*-1.0;
%         
%         H_ltp2(:,4) = Hequiv(1:8,17);
%         
%        [UKF_GDOP(Epoch_lo) UKF_PDOP(Epoch_lo) UKF_HDOP(Epoch_lo) UKF_VDOP(Epoch_lo) UKF_TDOP(Epoch_lo)] = GARD_CalculateDOPS(H_ltp);
       %[UKF_GDOP2(Epoch_lo) UKF_PDOP2(Epoch_lo) UKF_HDOP2(Epoch_lo) UKF_VDOP2(Epoch_lo) UKF_TDOP2(Epoch_lo)] = GARD_CalculateDOPS(H_ltp2);
        

        
        % calculate kalman gain
        UKF_K_k = UKF_Pxy_kminus * inv(UKF_Py_kminus);

        
        % apply correction
        UKF_z = (UKF_y_k - UKF_y_hat_kminus);
        
        % wrap the 17th element
        UKF_z(17) = pibound(UKF_z(17));

        UKF_z_save(Epoch_lo,:) = UKF_z;

        % Uncomment this for GPS outage
%         if Epoch_lo > 100 && Epoch_lo < 110
%             UKF_x_hat_kplus = UKF_x_hat_kminus;
%             UKF_Px_kplus = UKF_Px_kminus;
%         else
            UKF_x_hat_kplus = UKF_x_hat_kminus + UKF_K_k * UKF_z;
            UKF_Px_kplus = UKF_Px_kminus - UKF_K_k * UKF_Py_kminus * UKF_K_k';
%         end
        
        % normalise the attitude quaternion
          UKF_x_hat_kplus(7:10) = UKF_x_hat_kplus(7:10)/norm(UKF_x_hat_kplus(7:10));

         %%%% REMOVE THIS  %%%%%
         % quat truth
%          quat_truth = EulerToQuat(att_truth(2:4,Epoch_hi))';
%          UKF_x_hat_kplus(7:10) = quat_truth;
         
         
        %save results
        t_save(Epoch_hi) = Epoch_hi/100;
        UKF_x_hat_save(Epoch_hi,:) = UKF_x_hat_kplus;
        UKF_P_save(Epoch_hi,:) = diag(UKF_Px_kplus);
    
        % save the post-update residual (shoudl be whiteish)
        [PR_Vec_plus, PRR_Vec_plus] = GARD_GetMeasurementPrediction(UKF_x_hat_kplus,NumberGPSMeasurements,SVPos,SVVel);
        CompassHeading_plus = GARD_HeadingFromQuaternion(UKF_x_hat_kplus(7:10));
        UKF_y_hat_kplus = [PR_Vec_plus; PRR_Vec_plus; CompassHeading_plus];
        UKF_z_post_save(Epoch_lo,:) = UKF_y_k - UKF_y_hat_kplus;
        
        % copy updated state to high-speed loop state for next round
        UKF_x_hat_kminus = UKF_x_hat_kplus;
        UKF_Px_kminus = UKF_Px_kplus;

if Sim_Data == 0
        Clock_error(Epoch_lo) = UKF_x_hat_kplus(17) - UserPos_T(Epoch_lo);
else
        Clock_error(Epoch_lo) = UKF_x_hat_kplus(17);
end
        
        
%% calculate sub solutions
if CalculateSubSolutions

        NumberSubMeasurements = NumberGPSMeasurements-1;
        NumberSubFilters = nchoosek(NumberGPSMeasurements,NumberSubMeasurements);
        
        

                
                
        for SubSolution = 1:NumberSubFilters
            PR_SubVec(1:SubSolution-1) = PR_Vec(1:SubSolution-1);
            PR_SubVec(SubSolution:NumberSubMeasurements) = PR_Vec(SubSolution+1:NumberGPSMeasurements);

            PRR_SubVec(1:SubSolution-1) = PRR_Vec(1:SubSolution-1);
            PRR_SubVec(SubSolution:NumberSubMeasurements) = PRR_Vec(SubSolution+1:NumberGPSMeasurements);
            
            % formulate the measurment vector
            UKF_y_k_Sub = [PR_SubVec(1:NumberSubMeasurements)';PRR_SubVec(1:NumberSubMeasurements)';CompassSensor(Epoch_hi)];

%             
%             PR_Vec_minus_0_Sub(1:SubSolution-1) = PR_Vec_minus_0(1:SubSolution-1);
%             PR_Vec_minus_0_Sub(SubSolution:NumberSubMeasurements) = PR_Vec_minus_0(SubSolution+1:NumberGPSMeasurements);
%             
%             PRR_Vec_minus_0_Sub(1:SubSolution-1) = PRR_Vec_minus_0(1:SubSolution-1);
%             PRR_Vec_minus_0_Sub(SubSolution:NumberSubMeasurements) = PRR_Vec_minus_0(SubSolution+1:NumberGPSMeasurements);
%             
%             PR_Vec_minus_i_Sub(1:SubSolution-1,:) = PR_Vec_minus_i(1:SubSolution-1,:);
%             PR_Vec_minus_i_Sub(SubSolution:NumberSubMeasurements,:) = PR_Vec_minus_i(SubSolution+1:NumberGPSMeasurements,:);
%             
%             PRR_Vec_minus_i_Sub(1:SubSolution-1,:) = PRR_Vec_minus_i(1:SubSolution-1,:);
%             PRR_Vec_minus_i_Sub(SubSolution:NumberSubMeasurements,:) = PRR_Vec_minus_i(SubSolution+1:NumberGPSMeasurements,:);
%              

            %% calculate PR measurement sigma points by excluding the
            %% measurement from the sub filter number.  i.e. if sub
            %% filter=1 then exlude PR_Vec(1);
            % k is the index of the SV and PR in the SV_Vec and PR_Vec
            % j is the index of where the PR is stored
            
            for j=1:NumberSubMeasurements
                if j<SubSolution
                    k = j;
                else
                    k = j+1;
                end
                
                % SVPos_Sub needed for dops
                SVPos_Sub(j,:) = SVPos(k,:);
                
                
                [PR_Vec_minus_0_Sub(j,SubSolution),PRR_Vec_minus_0_Sub(j,SubSolution)] = ...
                    GARD_GPSMeasurementPrediction(...
                       [xs_Sub(1:3,1,SubSolution); xs_Sub(17,1,SubSolution)],...
                       [xs_Sub(4:6,1,SubSolution); xs_Sub(18,1,SubSolution)],...
                       SVPos(k,:),...
                       SVVel(k,:));

                PR_Vec_minus_0_Sub(j,SubSolution) = PR_Vec_minus_0_Sub(j,SubSolution) + ...
                                                    xs_0(NumberStates+ProcessNoiseStates+j);
                PRR_Vec_minus_0_Sub(j,SubSolution) = PRR_Vec_minus_0_Sub(j,SubSolution) + ...
                                                    xs_0(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+j);

               for i=1:2*(Na-2) %% two measurements are missing
                   [PR_Vec_minus_i_Sub(j,i,SubSolution),PRR_Vec_minus_i_Sub(j,i,SubSolution)] = ...
                        GARD_GPSMeasurementPrediction(...
                           [xs_Sub(1:3,i+1,SubSolution); xs_Sub(17,i+1,SubSolution)],...
                           [xs_Sub(4:6,i+1,SubSolution); xs_Sub(18,i+1,SubSolution)],...
                           SVPos(k,:),...
                           SVVel(k,:));

                    PR_Vec_minus_i_Sub(j,i,SubSolution) = PR_Vec_minus_i_Sub(j,i,SubSolution) + ...
                                                        xs_Sub(NumberStates+ProcessNoiseStates+j,i+1,SubSolution);
                    PRR_Vec_minus_i_Sub(j,i,SubSolution) = PRR_Vec_minus_i_Sub(j,i,SubSolution) + ...
                                                        xs_Sub(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+j,i+1,SubSolution);     
               end
                                                
            end
            
            CompassHeading_minus_0_Sub(SubSolution) = GARD_HeadingFromQuaternion(xs_Sub(7:10,1,SubSolution)) + ...
                    xs_Sub(NumberStates+ProcessNoiseStates+2*NumberSubMeasurements+1,1,SubSolution);
            for i=1:2*(Na-2) %% two measurements are missing
                CompassHeading_minus_i_Sub(i,SubSolution) = GARD_HeadingFromQuaternion(xs_Sub(7:10,i+1,SubSolution)) + ...
                    xs_Sub(NumberStates+ProcessNoiseStates+2*NumberSubMeasurements+1,i+1,SubSolution);
            end
            
            ys_kminus_0_Sub = [PR_Vec_minus_0_Sub(:,SubSolution); PRR_Vec_minus_0_Sub(:,SubSolution);CompassHeading_minus_0_Sub(:,SubSolution)];
            ys_kminus_i_Sub = [PR_Vec_minus_i_Sub(:,:,SubSolution); PRR_Vec_minus_i_Sub(:,:,SubSolution);CompassHeading_minus_i_Sub(:,SubSolution)'];

            
            
            ys_Sub = [ys_kminus_0_Sub,ys_kminus_i_Sub];
            [UKF_y_hat_kminus_Sub,UKF_Py_kminus_Sub,UKF_Pxy_kminus_Sub] = ...
                GARD_GenerateMeasurementMeanFromSigmaPoints(...
                    xs_Sub(:,:,SubSolution),...
                    ys_Sub,...
                    Ws_Sub(:,SubSolution),...
                    NumberStates,...
                    NumberSubMeasurements);
            
            % calculate kalman gain
            UKF_K_k_Sub = UKF_Pxy_kminus_Sub * inv(UKF_Py_kminus_Sub);

        
            % apply correction
            UKF_z_Sub = (UKF_y_k_Sub - UKF_y_hat_kminus_Sub);
            
            UKF_x_hat_kplus_Sub(:,SubSolution) = UKF_x_hat_kminus_Sub(:,SubSolution) + UKF_K_k_Sub * UKF_z_Sub;

            UKF_Px_kplus_Sub(:,:,SubSolution) = UKF_Px_kminus_Sub(:,:,SubSolution) - UKF_K_k_Sub * UKF_Py_kminus_Sub * UKF_K_k_Sub';
            
            
            %% calculate DOPS of this sub soluiton
            DOP_Sub(Epoch_lo,:,SubSolution) = GARD_CalculateDOPS2(UserPos,SVPos_Sub);
            
            %% save a copy of the results for analysis
            UKF_x_hat_save_Sub(Epoch_lo,:,SubSolution) = UKF_x_hat_kplus_Sub(1:3,SubSolution);
            UKF_P_save_Sub(Epoch_lo,:,SubSolution) = diag(UKF_Px_kplus_Sub(1:3,1:3,SubSolution));
            
            
            %% update kminus for next loop
            UKF_x_hat_kminus_Sub(:,SubSolution) = UKF_x_hat_kplus_Sub(:,SubSolution);
            UKF_Px_kminus_Sub(:,:,SubSolution) = UKF_Px_kplus_Sub(:,:,SubSolution);
        end % end sub solutions


        
        for j = 1:NumberSubFilters
            % note that only the position estimates are used - not
            % clock or velocity
            P_ned = UKF_Px_kplus(1:3,1:3);
            P_ned_Sub = UKF_Px_kplus_Sub(1:3,1:3,j);

            NED_ss = (UKF_x_hat_kplus(1:3) - UKF_x_hat_kplus_Sub(1:3,j));
            UKF_Beta_ss_H(:,j) = NED_ss(1:2);
            UKF_B_ss_H(:,:,j) = abs(P_ned_Sub(1:2,1:2) - P_ned(1:2,1:2));

            UKF_lambda_ss_H(Epoch_lo,j) = UKF_Beta_ss_H(:,j)' * UKF_Beta_ss_H(:,j);

            UKF_B_lambda_H = eigs(UKF_B_ss_H(:,:,j));
            UKF_TD_H(Epoch_lo,j) = sqrt(max(UKF_B_lambda_H)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));

            if(UKF_lambda_ss_H(Epoch_lo,j) > UKF_TD_H(Epoch_lo,j))
                disp(sprintf('%s H-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
            end

            UKF_Beta_ss_V(:,j) = NED_ss(3);
            UKF_B_ss_V(:,:,j) = abs(P_ned_Sub(3,3) - P_ned(3,3));

            UKF_lambda_ss_V(Epoch_lo,j) = UKF_Beta_ss_V(:,j)' * UKF_Beta_ss_V(:,j);

            UKF_B_lambda_V = eigs(UKF_B_ss_V(:,:,j));
            UKF_TD_V(Epoch_lo,j) = sqrt(max(UKF_B_lambda_V)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));

            if(UKF_lambda_ss_V(Epoch_lo,j) > UKF_TD_V(Epoch_lo,j))
                disp(sprintf('%s V-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
            end

            %% this HPL represents the fault-free (H0) hypothesis. 
            Sigma_MD_H = abs(norminv(P_MD_H,0,1));
            Sigma_MD_V = abs(norminv(P_MD_V,0,1));
            UKF_HPL_sub(Epoch_lo,j) = Sigma_MD_H * sqrt(P_ned(1,1) + P_ned(2,2)) + UKF_TD_H(Epoch_lo,j);
            UKF_VPL_sub(Epoch_lo,j) = Sigma_MD_V * sqrt(P_ned(3,3)) + UKF_TD_V(Epoch_lo,j);

            %% todo - calculate the Fault-in progress (h1) hypothesis based
            %% HPL_H1.


            

        end

        UKF_HPL(Epoch_lo) = max(UKF_HPL_sub(Epoch_lo,:));
        UKF_VPL(Epoch_lo) = max(UKF_VPL_sub(Epoch_lo,:));
end %% if calculatesubsolutions        
        
        %% calculate solution using least squares for comparisson
    [LSQ_SolutionVec(Epoch_lo,:), LSQ_Variance(Epoch_lo,:), LSQ_NumIterations(Epoch_lo), ...
        LSQ_ResidualVector(Epoch_lo,:), LSQ_M, LSQ_Fail(Epoch_lo), ...
        LSQ_limit(Epoch_lo), LSQ_DOP(Epoch_lo,:)] = GARD_LSQ([ApproxPos 0],NumberGPSMeasurements, PR_Vec_raw(1:NumberGPSMeasurements),SVPos(1:NumberGPSMeasurements,:));
   
    [LSQVel_SolutionVec(Epoch_lo,:), LSQVel_VarSolutionVec(Epoch_lo,:), LSQVel_NumIterations(Epoch_lo,:),...
            LSQVel_ResidualVector(Epoch_lo,1:NumberGPSMeasurements), LSQVel_M, LSQVel_Fail(Epoch_lo,:), ...
            LSQVel_limit(Epoch_lo,:)] = GARD_LSQVel(UserPos,UserVel,NumberGPSMeasurements,PRR_Vec(1:NumberGPSMeasurements),SVPos(1:NumberGPSMeasurements,:), SVVel(1:NumberGPSMeasurements,:));

        [LSQ_BadGeometry(Epoch_lo), LSQ_RAIM_ALERT(Epoch_lo), LSQ_SLOPE_Max(Epoch_lo), LSQ_r(Epoch_lo), LSQ_Td(Epoch_lo), ...
       LSQ_RAIM_HPL(Epoch_lo),LSQ_RAIM_VPL(Epoch_lo), LSQ_FaultySatFDI(Epoch_lo)] = ...
       GARDSim_RAIMParity(RAIM_a, RAIM_lambdatrue, NumberGPSMeasurements,PFalseAlarm,sigma_pr,556,LSQ_ResidualVector(Epoch_lo,:)',LSQ_M*Tecef2ned2');
       
    
    LSQ_SSE(Epoch_lo) = LSQ_ResidualVector(Epoch_lo,:)*LSQ_ResidualVector(Epoch_lo,:)';
    
        LSQVel_NED(Epoch_lo,:) = (Tecef2ned*LSQVel_SolutionVec(Epoch_lo,1:3)');
        
        [LSQ_LLH(Epoch_lo,1) LSQ_LLH(Epoch_lo,2) LSQ_LLH(Epoch_lo,3)] = ECEF2LLH(LSQ_SolutionVec(Epoch_lo,1:3));

        LSQ_LLH_error(Epoch_lo,:) = LSQ_LLH(Epoch_lo,:) - pos_truth_llh(2:4,Epoch_hi)';
        LSQ_ECEF_error(Epoch_lo,:) = LSQ_SolutionVec(Epoch_lo,1:3) - pos_truth_ecef(2:4,Epoch_hi)';
        LSQVel_error(Epoch_lo,:) = LSQVel_NED(Epoch_lo,:) - vel_truth(2:4,Epoch_hi)';
        
        if Sim_Data == 0
            LSQ_Clock_error(Epoch_lo,1) = LSQ_SolutionVec(Epoch_lo,4) - UserPos_T(Epoch_lo);
            LSQ_Clock_error(Epoch_lo,2) = LSQVel_SolutionVec(Epoch_lo,4);
        else
            LSQ_Clock_error(Epoch_lo,1) = LSQ_SolutionVec(Epoch_lo,4);
            LSQ_Clock_error(Epoch_lo,2) = LSQVel_SolutionVec(Epoch_lo,4);
        end
        if MakeMovie == 1
            plot_separation;            
            current_frame = getframe(gcf);
            aviobj = addframe(aviobj,current_frame);
            close(gcf);
        end


                    
    else
        
        % not a low-speed loop
        t_save(Epoch_hi) = Epoch_hi/100;
        UKF_x_hat_save(Epoch_hi,:) = UKF_x_hat_kminus;
        UKF_P_save(Epoch_hi,:) = diag(UKF_Px_kminus);
        
        
        
    end% end low speed loop

    hgt = UKF_x_hat_kminus(3);
    Lat_error(Epoch_hi) = (UKF_x_hat_kminus(1) - pos_truth_llh(2,Epoch_hi)) * (RM+hgt);
    Long_error(Epoch_hi) = (UKF_x_hat_kminus(2) - pos_truth_llh(3,Epoch_hi)) * ((RP+hgt) * cos(UKF_x_hat_kminus(1)));
    Height_error(Epoch_hi) = (UKF_x_hat_kminus(3) - pos_truth_llh(4,Epoch_hi));

    Vel_error(Epoch_hi,1:3) = UKF_x_hat_kminus(4:6) - vel_truth(2:4,Epoch_hi);
    
    
    % calculate attitde error
    phi_err(Epoch_hi) = (att_truth(2,Epoch_hi) - phi_q(Epoch_hi));
    theta_err(Epoch_hi) = (att_truth(3,Epoch_hi) - theta_q(Epoch_hi));
    psi_err(Epoch_hi) = (att_truth(4,Epoch_hi) - psi_q(Epoch_hi));

    if(psi_err(Epoch_hi) > pi)
        psi_err(Epoch_hi) = psi_err(Epoch_hi) - 2*pi;
    end
    if(psi_err(Epoch_hi) < -pi)
        psi_err(Epoch_hi) = psi_err(Epoch_hi),+ 2*pi;
    end

    if(mod(Epoch_hi,100) == 0)
        logline = sprintf('%s Completed Epoch %d of %d (%3.1f percent)\n',ModuleName,Epoch_lo, StopTime, 100*(Epoch_lo-StartTime)/(StopTime-StartTime));
        disp(logline);
        
        % write to progress log file
        fid = fopen(progress_filename,'a');
        fprintf(fid,logline);
        fclose(fid);
        
    end
    
    if(mod(Epoch_hi,1000) == 0)
        disp(sprintf('%s[%d] Saving Data',ModuleName,Epoch_hi));
        
        % save results
        save(ukf_results_filename);
    end
end % end high-speed loop

if MakeMovie == 1
    aviobj = close(aviobj);
end

% save results
save(ukf_results_filename);
logline = sprintf('%s Saving Results to %s\n',ModuleName,ukf_results_filename);
fid = fopen(progress_filename,'a');
fprintf(fid,logline);
fclose(fid);

disp('[GARDSim_GPSINS_UKF] Done!');

% plot results
GARDSim_PlotResults;


