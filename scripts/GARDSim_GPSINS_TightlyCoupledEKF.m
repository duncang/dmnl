% GPS-INS position solution using the EKF
% Written by Duncan Greer 29 May 2007
% $Id: GARDSim_GPSINS_TightlyCoupledEKF.m 1782 2008-06-19 04:32:16Z greerd $
%
%
% ==== STATES ====
%
% x1 = Latitude position (rad)
% x2 = Longitude position (rad)
% x3 = Height position (m)
% x4 = North velocity (m/s)
% x5 = East velocity (m/s)
% x6 = Down velicty (m/s)
% x7 = north tilt error (rad)
% x8 = east tilt error (rad)
% x9 = vertical tilt error (rad)
% x10 = x acc bias error (m/s/s)
% x11 = y acc bias error (m/s/s)
% x12 = z acc bias error (m/s/s)
% x13 = x gyro bias error (m/s/s)
% x14 = y gyro bias error (m/s/s)
% x15 = z gyro bias error (m/s/s)
% x16 = Clock Bias (m or sec?)
% x17 = Clock Drift (m/s or sec/sec?)
%
% Pos_LLH gives the position solution in LLH coordinates
% Vel_NED gives the velocity solution in NED coordinates
% The kalman filter states are tracked in ECEF coordinates (position and
% velocity) which are then converted to P and V
% UserPos and UserVel are the user position and veolocity in ECEF
% coordinates


ModuleName = '[GARDSim_GPSINS_TightlyCoupledEKF] ';


% setup constants
WGS84Constants;
d2r = pi/180;
r2d = 180/pi;
g = -9.79;  % m/s/s
GPSConstants;

NumberStates = 17;


% load data path
%DataPath = 'data/rnav_approach/';
%DataPath = 'data/DebugAeromodelpqrchanged11.1.08/';
%DataPath = 'data/DebugAeroModel7.1.08/';
%DataPath = 'data/Feb0108FlightNoWind/';
%DataPath = 'data/rnav_approach_long/';

DataPath = 'data/Simulator_Data/Apr_08/data/rnav_approach_long/';

% save data path
dtg = clock;
ekf_results_filename = strcat(DataPath,'ekf_results_');
ekf_results_filename = sprintf('%s%d%02d%02d_%02d%02d.mat',...
    ekf_results_filename,...
    dtg(1),...
    dtg(2),...
    dtg(3),...
    dtg(4),...
    dtg(5));
    

progress_filename = strcat(DataPath,'progress.txt');


disp(strcat(ModuleName,'Loading data from: "',DataPath,'"'));



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
if(~exist('pos_truth_llh','var'))
    load(strcat(DataPath,'pos_truth_llh.mat'));
    load(strcat(DataPath,'vel_truth.mat'));
    load(strcat(DataPath,'att_truth.mat'));
end

disp(strcat(ModuleName,'Finished loading data'));


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

gps_dt = 1;
ins_dt = 0.01;



StartTime = TruthOffset+300; % start at 60 seconds into the data set
StopTime = size(PR_Sim,1);%TruthOffset+1000;%


% initialise epochs
% initialise epochs
Epoch_lo = StartTime / gps_dt - 1;
Epoch_hi = StartTime / ins_dt - 1;

% position
%InitialPosition = [-27.4*pi/180 153.2*pi/180 4500/3.28]';  % somewhere near brisbane in LLH, rads and meters
InitialPosition = [pos_truth_llh(2,StartTime/ins_dt+1),pos_truth_llh(3,StartTime/ins_dt+1),pos_truth_llh(4,StartTime/ins_dt+1)]';
% velocity
%InitialVelocity = [0 -60 0]';  % NED velocities in m/s
InitialVelocity = [vel_truth(2,StartTime/ins_dt+1),vel_truth(3,StartTime/ins_dt+1),vel_truth(4,StartTime/ins_dt+1)]';
% attitude
%InitialAttitude = [0 0 270*pi/180]';  % roll, pitch, yaw in radians
InitialAttitude = [att_truth(2,StartTime/ins_dt+1),att_truth(3,StartTime/ins_dt+1),att_truth(4,StartTime/ins_dt+1)]';

C_BN = GARD_EulerToDCM(InitialAttitude(1),InitialAttitude(2),InitialAttitude(3));
% 

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));
RMh = RM + InitialPosition(3);
RPh = RP + InitialPosition(3);


NumberSVs = 32;

SVDontUse = zeros(1,NumberSVs);

%SVDontUse(11) = 1;

%% artificially inflate the PR noise at a certain time
NoiseOnTime = StartTime + 500;
NoiseOffTime = StopTime;
NoiseAmplification = 1;

%PR_Sim(NoiseOnTime:NoiseOffTime,:) = PR_Sim(NoiseOnTime:NoiseOffTime,:) + randn(NoiseOffTime-NoiseOnTime+1,NumberSVs) * NoiseAmplification;


LOCAL_GRAVITY = 9.80;

NumberEpochsGPS = min(StopTime,length(PR_Sim));

TimeINS = 0:ins_dt:NumberEpochsGPS;
TimeGPS = 0:gps_dt:NumberEpochsGPS-1;


NumberEpochsINS = min(length(TimeINS),StopTime/ins_dt + 1);

% initialise storage vars to svae time
Lat_error = zeros(1,NumberEpochsINS);
Long_error = zeros(1,NumberEpochsINS);
Height_error = zeros(1,NumberEpochsINS);
Clock_error = zeros(NumberEpochsGPS,1);

Vel_error = zeros(NumberEpochsINS,3);

Pos_LLH = zeros(NumberEpochsINS,3);
Vel_NED = zeros(NumberEpochsINS,3);
Acc_NED = zeros(NumberEpochsINS,3);

UserClock = zeros(NumberEpochsINS,2);

q_INS = zeros(NumberEpochsINS,4);

phi_q = zeros(NumberEpochsINS,1);
theta_q = zeros(NumberEpochsINS,1);
psi_q = zeros(NumberEpochsINS,1);

phi_uc = zeros(NumberEpochsINS,1);
theta_uc = zeros(NumberEpochsINS,1);
psi_uc = zeros(NumberEpochsINS,1);

gvec_phi = zeros(NumberEpochsINS,1);
gvec_theta = zeros(NumberEpochsINS,1);

lat_dot = zeros(NumberEpochsINS,1);
long_dot = zeros(NumberEpochsINS,1);
h_dot = zeros(NumberEpochsINS,1);

omega_b = zeros(NumberEpochsINS,3);
A_b = zeros(NumberEpochsINS,3);



Coriolis = zeros(NumberEpochsINS,3);

P_out_diag = zeros(NumberStates,NumberEpochsINS);
Pos_var = zeros(NumberEpochsINS,3);

SatellitesUsed = zeros(NumberEpochsGPS,1);

% setup error values
if(GRASOn == 1)
    PositionError = 2.1;
    VelocityError = 0.5;
    ClockBiasError = 1E2;
    ClockDriftError = 1.0;
    RangeNoiseVariance = 2.1^2; % m
    RangeRateNoiseVariance = 0.5^2; % m/s
else
    PositionError = 5.0;
    VelocityError = 2.0;
    ClockBiasError = 1E2;
    ClockDriftError = 1.0;
    RangeNoiseVariance = 5.0^2; % m
    RangeRateNoiseVariance = 0.8^2; % m/s
end



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


x_hat_out_full = zeros(NumberStates,1);
x_hat_minus = x_hat_out_full;

PositionVariance = PositionError ^ 2;
VelocityVariance = VelocityError ^ 2;
ClockBiasVariance = ClockBiasError ^ 2;
ClockDriftVariance = ClockDriftError ^2;
%

sigma_pr_gnd = 2.0;
sigma_prr = 0.5;

if GRASOn == 0
    sigma_iono = 1.5;
else
    sigma_iono = 0.1;
end

sigma_pr = sqrt(sigma_pr_gnd^2 + sigma_iono^2);


P_in = zeros(NumberStates,NumberStates);


P_in(1,1) = PositionVariance/RM;
P_in(2,2) = PositionVariance/RP*cos(-0.47);
P_in(3,3) = 2*PositionVariance;

P_in(4:6,4:6) = eye(3,3)*VelocityVariance;

P_in(7:9,7:9) = (1*eye(3,3)*pi/180)^2;
P_in(10:12,10:12) = 0.0001*eye(3,3);
P_in(13:15,13:15) = 0.0001*eye(3,3);
P_in(16,16) = ClockBiasVariance;
P_in(17,17) = ClockDriftVariance;

P_in_full = P_in;

%% initialise sub filter
for j=1:10
    P_in_sub(:,:,j) = P_in;
end





phi_true  = downsample(att_truth(2,:),1);
phi_true = phi_true(1:length(phi_q));

theta_true = downsample(att_truth(3,:),1);
theta_true = theta_true(1:length(theta_q));

psi_true = downsample(att_truth(4,:),1);
psi_true = psi_true(1:length(psi_q));


%% Pseudorange Error
PR_Error = 0;
PR_Error_Rate = 0.0;

% set the threshold for fault detection
%FaultThreshold = 100;

%% FD is not turned on until the solution has reached stability - about 1
%% minute
FDEnabled = 0;
ExcludeSV = 0;
ExcludeSVPRN = 0;

% GyroBiasTruth=[0.1,0.15,0.2];
% AccelBiasTruth=[0.0,0.0,0.0];
% GyroBias = zeros(NumberEpochsGPS,3);
% AccelBias = zeros(NumberEpochsGPS,3);
% 
% for i=1:NumberEpochsGPS
%     GyroBias(i,:) = GyroBiasTruth;%+[0.01,-0.01,0.02];
%     AccelBias(i,:) = AccelBiasTruth;%+[0.1,0.1,-0.1];
% end

GyroBias = zeros(NumberEpochsGPS,3);
AccelBias = zeros(NumberEpochsGPS,3);

%GyroBiasTruth=[0.01,-0.02,0.005];
GyroBiasTruth = zeros(3,NumberEpochsINS);
GyroNoiseTruth = randn(3,NumberEpochsINS)*0.00000;

for i=1:NumberEpochsINS
%     GyroBiasTruth(1,i) = sin((1/10000)*i)*0.00001;
%     GyroBiasTruth(2,i) = sin((1/15000)*i)*-0.00002;
%     GyroBiasTruth(3,i) = sin((1/9000)*i)*0.0000001;
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


PFalseAlarm = 1e-5; %
PMissedDetection = 0.048;  % 95%

P_MD_Vert_Proportion = 0.50;
P_MD_H = (1-P_MD_Vert_Proportion) * PMissedDetection;
P_MD_V = P_MD_Vert_Proportion * PMissedDetection;


% two loops will run - a high speed loop calculating teh INS solution at
% 50 Hz - The data is collected at 100Hz, and a Runge-Kutta integration
% scheme is used 

Pfa = 2.22e-8;
SV_TrackTime = zeros(NumberSVs,1);

disp(strcat(ModuleName,sprintf('Beginning navigation loop for %d Epochs',NumberEpochsINS)));

for Epoch_hi = StartTime/ins_dt+1:StopTime/ins_dt+1
    
%    LOCAL_GRAVITY = gravity(2,Epoch_hi);
    %g = -9.79; 
    
    % implement the inertial navigation solution here
    if Epoch_hi == StartTime/ins_dt+1
%        % first time step

        q_INS(Epoch_hi,:) = EulerToQuat(InitialAttitude);
        %C_BN = GARD_EulerToDCM(InitialAttitude(1),InitialAttitude(2),InitialAttitude(3));
        C_BN = GARD_QuatToDCM(q_INS(Epoch_hi,:));
        c32 = C_BN(3,2);
        c33 = C_BN(3,3);
        c31 = C_BN(3,1);
        c21 = C_BN(2,1);
        c11 = C_BN(1,1);
    
        % calculate euler angles for reference
        phi_q(Epoch_hi) = atan2(C_BN(3,2),C_BN(3,3));
        theta_q(Epoch_hi) = asin(-C_BN(3,1));
        psi_q(Epoch_hi) = atan2(C_BN(2,1),C_BN(1,1));
        
        
        Pos_LLH(Epoch_hi,:) = InitialPosition;       
        Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch_hi,1),Pos_LLH(Epoch_hi,2));
        Tned2ecef = Tecef2ned';
        
        Vel_NED(Epoch_hi,:) = InitialVelocity';% + Acc_NED(Epoch_hi,:)*ins_dt;

        lat_dot(Epoch_hi) = Vel_NED(Epoch_hi,1) / (RM + Pos_LLH(Epoch_hi,3));
        long_dot(Epoch_hi) = Vel_NED(Epoch_hi,2) / (cos(Pos_LLH(Epoch_hi,1))*(RP + Pos_LLH(Epoch_hi,3)));
        h_dot(Epoch_hi) =  -Vel_NED(Epoch_hi,3);

        Pos_LLH(Epoch_hi,:) = InitialPosition;% + [lat_dot(Epoch_hi);long_dot(Epoch_hi);h_dot(Epoch_hi)]*ins_dt;
        
        x_out = [InitialPosition;InitialVelocity;EulerToQuat(InitialAttitude)'];
        
        F_INS = zeros(NumberStates,NumberStates);
        F_INS(1:9,1:9) = GARD_GenerateINSFMatrix(Pos_LLH(Epoch_hi,:),Vel_NED(Epoch_hi,:),Acc_NED(Epoch_hi,:)-[0 0 LOCAL_GRAVITY]);

        F_INS(16,17) = 1.0;
        
        phi2 = eye(NumberStates,NumberStates);
%        phi2(1:9,1:9) = phi2(1:9,1:9) + F_INS * ins_dt;
        phi2 = expm(F_INS*ins_dt);
        
        
        
        % % augmented state matrix
        phi2(7:9,13:15) = -C_BN*ins_dt; %% attitude to gyro
        phi2(4:6,10:12) = C_BN*ins_dt; %%  velocity to acceleration

        phi2(10,10) = exp( - x_accel_beta * ins_dt);
        phi2(11,11) = exp( - y_accel_beta * ins_dt);
        phi2(12,12) = exp( - z_accel_beta * ins_dt);
        phi2(13,13) = exp( - x_gyro_beta * ins_dt);
        phi2(14,14) = exp( - y_gyro_beta * ins_dt);
        phi2(15,15) = exp( - z_gyro_beta * ins_dt);
        
        
        phi2(16,17) = 1.0*ins_dt;
        Sc = 1*0.14;
        Sf = 1*0.0359;
        
        Q = zeros(NumberStates,NumberStates);

        Q(10,10) = 2*x_accel_beta*x_accel_Q^2;
        Q(11,11) = 2*y_accel_beta*y_accel_Q^2;
        Q(12,12) = 2*z_accel_beta*z_accel_Q^2;
        Q(13,13) = 2*x_gyro_beta*x_gyro_Q^2;
        Q(14,14) = 2*y_gyro_beta*y_gyro_Q^2;
        Q(15,15) = 2*z_gyro_beta*z_gyro_Q^2;

        Q(16,16) = (Sc)*ins_dt + Sf*(ins_dt^3)/3;
        Q(17,17) = Sf*ins_dt;
        Q(16,17) = Sf*(ins_dt^2)/2;
        Q(17,16) = Sf*(ins_dt^2)/2;
        



        G = zeros(NumberStates,NumberStates);
        G(7:9,13:15) = -C_BN;  %% attitude to gyro
        G(4:6,10:12) = C_BN; %%  velocity to acceleration

        G(10:12,10:12) = eye(3,3);
        G(13:15,13:15) = eye(3,3);

        G(16,16) = 1;
        G(16,17) = 1;
        G(17,16) = 0;
        G(17,17) = 1;

        Q2d = phi2 * (G * Q * G') * phi2' * ins_dt;

        Q2d(16,16) = 1*Sc*ins_dt+Sf*(ins_dt^3)/3;
        Q2d(16,17) = 1*Sf*(ins_dt^2)/2;
        Q2d(17,16) = 1*Sf*(ins_dt^2)/2;
        Q2d(17,17) = 1*Sf*ins_dt;
        
        x_hat_minus = phi2 * x_hat_minus;
        % propagate the covariance matrix
        P_in_full = phi2 * P_in_full * phi2' + Q2d;
        
        P_out_diag(:,Epoch_hi) = diag(P_in_full);
        %% calculate error bounds
        Pos_var(Epoch_hi,:) = P_out_diag(1:3,Epoch_hi).*[RM^2;RP^2*cos(Pos_LLH(Epoch_hi,1))^2;1];
   else
        ins_dt = 0.01;% TimeINS(Epoch_hi) - TimeINS(Epoch_hi-1);
        
        % get sensor measurements
          if(Epoch_lo < 2)
             omega_x = sensors(5,Epoch_hi);
             omega_y = sensors(6,Epoch_hi);
             omega_z = sensors(7,Epoch_hi);
 
             A_xb = sensors(2,Epoch_hi);
             A_yb = sensors(3,Epoch_hi);
             A_zb = sensors(4,Epoch_hi);
         else
            omega_x = sensors(5,Epoch_hi) - GyroBias(Epoch_lo,1) + GyroBiasTruth(1,Epoch_hi) + GyroNoiseTruth(1,Epoch_hi);
            omega_y = sensors(6,Epoch_hi) - GyroBias(Epoch_lo,2) + GyroBiasTruth(2,Epoch_hi) + GyroNoiseTruth(2,Epoch_hi);
            omega_z = sensors(7,Epoch_hi) - GyroBias(Epoch_lo,3) + GyroBiasTruth(3,Epoch_hi) + GyroNoiseTruth(3,Epoch_hi);

            A_xb = sensors(2,Epoch_hi) - AccelBias(Epoch_lo,1) + AccelBiasTruth(1,Epoch_hi) + AccelNoiseTruth(1,Epoch_hi);
            A_yb = sensors(3,Epoch_hi) - AccelBias(Epoch_lo,2) + AccelBiasTruth(2,Epoch_hi) + AccelNoiseTruth(2,Epoch_hi);
            A_zb = sensors(4,Epoch_hi) - AccelBias(Epoch_lo,3) + AccelBiasTruth(3,Epoch_hi) + AccelNoiseTruth(3,Epoch_hi);
          end
         
          % save sensor inputs
          omega_b(Epoch_hi,:) = [omega_x,omega_y,omega_z];
          A_b(Epoch_hi,:) = [A_xb,A_yb,A_zb];
%         
        % generate current rotation matrix from earth to navigation frame
        % coordinates
        Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch_hi-1,1),Pos_LLH(Epoch_hi-1,2));
        Tned2ecef = Tecef2ned';
        
        % find gravity vector for attitude estimate
        gvec_phi(Epoch_hi) =  atan2(-A_yb,sqrt(A_xb^2 + A_zb^2));% roll
        gvec_theta(Epoch_hi) = atan2(A_xb,-A_zb); % pitch
        
        Acc_in(1) = A_xb;
        Acc_in(2) = A_yb;
        Acc_in(3) = A_zb;
        Omega_in(1) = omega_x;
        Omega_in(2) = omega_y;
        Omega_in(3) = omega_z;
                
%         %q = EulerToQuat([phi_q(Epoch_hi-1),theta_q(Epoch_hi-1),psi_q(Epoch_hi-1)]);
%          g = GravityModel(x_out(1:3));
%         [x_out,output] = GARD_INSMechanisation([Pos_LLH(Epoch_hi-1,:) Vel_NED(Epoch_hi-1,:) q_INS(Epoch_hi-1,:) zeros(1,6)]',ins_dt,Acc_in,Omega_in,g);
% 
%         lat_dot(Epoch_hi) = output.llhdot(1);
%         long_dot(Epoch_hi) = output.llhdot(2);
%         h_dot(Epoch_hi) = output.llhdot(3);
%         
%         q_INS(Epoch_hi,:) = x_out(7:10);
%         
%         C_BN = output.C_BN;
%         
%         phi_q(Epoch_hi) = atan2(C_BN(3,2),C_BN(3,3));
%         theta_q(Epoch_hi) = asin(-C_BN(3,1));
%         psi_q(Epoch_hi) = atan2(C_BN(2,1),C_BN(1,1));
%         
%         Vel_NED(Epoch_hi,:) = x_out(4:6);
%         Pos_LLH(Epoch_hi,:) = x_out(1:3);
%         
% 
%         
%         % clock update
%         UserClock(Epoch_hi,2) = UserClock(Epoch_hi-1,2);
%         UserClock(Epoch_hi,1) = UserClock(Epoch_hi-1,1) + UserClock(Epoch_hi-1,2)*ins_dt;
%      
%             % propagate the state error
%         % generate PHI matrix
%         F_INS = zeros(NumberStates,NumberStates);
%         F_INS(1:9,1:9) = GARD_GenerateINSFMatrix(Pos_LLH(Epoch_hi-1,:),Vel_NED(Epoch_hi-1,:),C_BN*Acc_in');
% 
%         F_INS(16,17) = 1.0;
%         
%         phi2 = eye(NumberStates,NumberStates);
% %        phi2(1:9,1:9) = phi2(1:9,1:9) + F_INS * ins_dt;
%         phi2 = expm(F_INS*ins_dt);
%         
%         
%         
%         % % augmented state matrix
%         phi2(7:9,13:15) = -C_BN*ins_dt; %% attitude to gyro
%         phi2(4:6,10:12) = C_BN*ins_dt; %%  velocity to acceleration
% 
%         phi2(10,10) = exp( - x_accel_beta * ins_dt);
%         phi2(11,11) = exp( - y_accel_beta * ins_dt);
%         phi2(12,12) = exp( - z_accel_beta * ins_dt);
%         phi2(13,13) = exp( - x_gyro_beta * ins_dt);
%         phi2(14,14) = exp( - y_gyro_beta * ins_dt);
%         phi2(15,15) = exp( - z_gyro_beta * ins_dt);
%         
%         
%         phi2(16,17) = 1.0*ins_dt;
% 
% 
% 
%         %Sf = 8*pi^2*2e-20*Speedoflight^2;  % frequency noise psd
%         %Sc = 2*2e-19*Speedoflight^2;  % bias noise psd
%         Sc = 1*0.14;
%         Sf = 1*0.0359;
%         
%         Q = zeros(NumberStates,NumberStates);
% 
%         Q(10,10) = 2*x_accel_beta*x_accel_Q^2;
%         Q(11,11) = 2*y_accel_beta*y_accel_Q^2;
%         Q(12,12) = 2*z_accel_beta*z_accel_Q^2;
%         Q(13,13) = 2*x_gyro_beta*x_gyro_Q^2;
%         Q(14,14) = 2*y_gyro_beta*y_gyro_Q^2;
%         Q(15,15) = 2*z_gyro_beta*z_gyro_Q^2;
% 
%         Q(16,16) = (Sc)*ins_dt + Sf*(ins_dt^3)/3;
%         Q(17,17) = Sf*ins_dt;
%         Q(16,17) = Sf*(ins_dt^2)/2;
%         Q(17,16) = Sf*(ins_dt^2)/2;
%         
% 
% 
% 
%         G = zeros(NumberStates,NumberStates);
%         G(7:9,13:15) = -C_BN;  %% attitude to gyro
%         G(4:6,10:12) = C_BN; %%  velocity to acceleration
% 
%         G(10:12,10:12) = eye(3,3);
%         G(13:15,13:15) = eye(3,3);
% 
%         G(16,16) = 1;
%         G(16,17) = 1;
%         G(17,16) = 0;
%         G(17,17) = 1;
% 
%         Q2d = phi2 * (G * Q * G') * phi2' * ins_dt;
% 
%         Q2d(16,16) = 1*Sc*ins_dt+Sf*(ins_dt^3)/3;
%         Q2d(16,17) = 1*Sf*(ins_dt^2)/2;
%         Q2d(17,16) = 1*Sf*(ins_dt^2)/2;
%         Q2d(17,17) = 1*Sf*ins_dt;
%         
%         x_hat_minus = phi2 * x_hat_minus;
%         % propagate the covariance matrix
%         P_in_full = phi2 * P_in_full * phi2' + Q2d;
        
        in.x = x_out;
        in.P = P_in_full;
        in.Acc = Acc_in;
        in.Omega = Omega_in;
        in.UserClock = UserClock(Epoch_hi-1,:);
        
        in.x_accel_beta = x_accel_beta;
        in.y_accel_beta = y_accel_beta;
        in.z_accel_beta = z_accel_beta;
        in.x_accel_Q = x_accel_Q;
        in.y_accel_Q = y_accel_Q;
        in.z_accel_Q = z_accel_Q;
        
        in.x_gyro_beta = x_gyro_beta;
        in.y_gyro_beta = y_gyro_beta;
        in.z_gyro_beta = z_gyro_beta;
        in.x_gyro_Q = x_gyro_Q;
        in.y_gyro_Q = y_gyro_Q;
        in.z_gyro_Q = z_gyro_Q;
        

        out = GARD_PropagateINSEKF(in,NumberStates,ins_dt);


        Pos_LLH(Epoch_hi,:) = out.x(1:3);
        Vel_NED(Epoch_hi,:) = out.x(4:6);
        phi_q(Epoch_hi) = out.euler(1);
        theta_q(Epoch_hi) = out.euler(2);
        psi_q(Epoch_hi) = out.euler(3);
        UserClock(Epoch_hi,:) = out.UserClock; 
        x_out = out.x;
        
        P_out_diag(:,Epoch_hi) = diag(out.P);
        %% calculate error bounds
        Pos_var(Epoch_hi,:) = P_out_diag(1:3,Epoch_hi).*[RM^2;RP^2*cos(Pos_LLH(Epoch_hi,1))^2;1];
        
        x_hat_minus = out.x;
        P_in_full = out.P;

        
%         if CalculateSubSolutions            
%             %% for each sub filter calculate the INS mechanisation
%             for SubSolution = 1:NumberSubFilters
%             
%             end
%         end % Calculate SubSolutions


    end % 
        % low-speed loop - kalman filter
        if(mod(Epoch_hi-1,100) == 0)% && Epoch_hi ~= 1)
        %if(0)  %% uncomment this line to turn off the KF update
            Epoch_lo = Epoch_lo + 1;
        
            if Epoch_lo < TruthOffset
                warning('no gps data.. continuing');
                continue;
            end
        
            if(Epoch_lo == StartTime + 60)  % enable fault detection
                FDEnabled = 1;
            end
            
            clear PR_Vec PRR_Vec H R K V z_Vec NumberGPSMeasurements PR_Vec_minus PRR_Vec_minus Tecef2ned Tned2ecef SV_Vec SVPos SVVel;

            % generate current rotation matrix from earth to navigation frame
            % coordinates
            Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch_hi,1),Pos_LLH(Epoch_hi,2));
            Tned2ecef = Tecef2ned';
        
            % get the gps time - used for the iono correction
            GPSTime  = 259199 + Epoch_lo;
            
            % simulate GPS outage
%             if(Epoch_lo > 150 && Epoch_lo < 200)
%                 continue;
%             end

            % UserPos and UserVel are the current state vector
            % which the EKF will linearise on
            UserPos(1:3) =  LLH2ECEF(Pos_LLH(Epoch_hi,1),Pos_LLH(Epoch_hi,2),Pos_LLH(Epoch_hi,3));
            UserPos(4) = UserClock(Epoch_hi,1);
            
            UserVel(1:3) = Tned2ecef * Vel_NED(Epoch_hi,:)';
            UserVel(4) = UserClock(Epoch_hi,2);
            
            
            % format the input vectors
            SVIndex = 0;
            for SV=1:NumberSVs
                if((PR_Sim(Epoch_lo,SV) ~= 0) && SVDontUse(SV) == 0)  %%SV ~= ExcludeSVPRN
                        % add to PR vector
                        SVIndex = SVIndex + 1;

                       SV_TrackTime(SV) = SV_TrackTime(SV) + 1;

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
                    [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV), URA(Epoch_lo,SV)] = ...
                        GPSOrbitPropagator(GPSTime_Week(Epoch_lo-TruthOffset+1), GPSTime_Sec(Epoch_lo-TruthOffset+1)-PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris, 7500);
    else
                    [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV), URA(Epoch_lo,SV)] = ...
                        GPSOrbitPropagator(GPSTime_Week(Epoch_lo), GPSTime_Sec(Epoch_lo), SV, SV_Ephemeris, 7500);    
    end
                    if ValidPosData(Epoch_lo,SV) ~= 1
                        %SVDontUse(SV) = 1;
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
                        %SVDontUse(SV) = 1;
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
                        TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),Pos_LLH(Epoch_hi,3));
                    end

                    PR_Vec(SVIndex) = PR_Vec(SVIndex) - TropoDelay(Epoch_lo,SV);

                end
                
          

            end

    %%%%% SOLUTION SEPARATION GPS-INS %%%%%%
            
            % get the number of GPS measurements -
            %   < 4 - can possibly do a GPS-INS solution using clock
            %   coasting, however could become unstable
            %   = 4 - GPS solution - no integrity
            %   = 5 - Solution Separation - special case
            %   > 5 - Solution Separation

            NumberGPSMeasurements = length(PR_Vec);
            if NumberGPSMeasurements < 8
                warning(sprintf('less than 8 measurements at epoch %d',Epoch_lo));
                continue;
            end
        
            if NumberGPSMeasurements > 8
                NumberGPSMeasurements = 8;
            end
             
            SatellitesUsed(Epoch_lo) = NumberGPSMeasurements;
            SV_Vec_save(Epoch_lo,1:NumberGPSMeasurements) = SV_Vec(1:NumberGPSMeasurements);
            
            %x_hat_in = x_hat(:,Epoch_hi);
            %x_hat_in = x_hat_out_full;
                        
                %%% FIRST DO THE FULL SOLUTION

                % R matrix
                R = eye(2*NumberGPSMeasurements);
                for k = 1:NumberGPSMeasurements
                    R(k,k) = sigma_pr^2;
                    R(k+NumberGPSMeasurements,k+NumberGPSMeasurements) = sigma_prr^2;
                end

                % determine the H matrix
                H = zeros(2*NumberGPSMeasurements,NumberStates);
                for k = 1:NumberGPSMeasurements
                    %Calculated slant ranges

                    for m = 1:3
                         ele(m) =  SVPos(k,m) - UserPos(m);
                    end    

                    r_VecCalc(k) =  norm(ele);   

                    H(k,1) =  -ele(1)/r_VecCalc(k);
                    H(k,2) =  -ele(2)/r_VecCalc(k);
                    H(k,3) =  -ele(3)/r_VecCalc(k);
                    H(k,16) = 1.0;   


                    H(k+NumberGPSMeasurements,4) = -ele(1)/r_VecCalc(k);
                    H(k+NumberGPSMeasurements,5) = -ele(2)/r_VecCalc(k);
                    H(k+NumberGPSMeasurements,6) = -ele(3)/r_VecCalc(k);
                    H(k+NumberGPSMeasurements,17) = 1.0;
                    
                    
                    % calculate hte earth rotation correction as per Kayton pg 228
                    % eq 5.67

                    delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) *UserPos(2) - SVPos(k,2) * UserPos(1));
                
                    % find apriori estimate of pseudorange
                    PR_Vec_minus(k) = r_VecCalc(k) + UserPos(4) - SVPos(k,4) - delta_pr_omegaedot(k);  % geometric range + c * delta_T
                    
                    %predicted relative velocity of sv and receiver
                    r_VecCalcVel(k) = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + ...
                                      (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + ...
                                      (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
                    Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);

                    PRR_Vec_minus(k) = Relative_Velocity(k) + UserVel(4) - SVVel(k,4);
                end
                
                % find measurement vector - delta between predicted and measured
                % pseudorange
                z_Vec = ([PR_Vec(1:NumberGPSMeasurements) PRR_Vec(1:NumberGPSMeasurements)]' - [PR_Vec_minus PRR_Vec_minus]');                
                
                
                H(1:NumberGPSMeasurements,1:3) = H(1:NumberGPSMeasurements,1:3)*Tecef2ned';
                H(NumberGPSMeasurements+1:2*NumberGPSMeasurements,4:6) = H(NumberGPSMeasurements+1:2*NumberGPSMeasurements,4:6)*Tecef2ned';
                
                %% add compass measurement
                Heading_error = CompassSensor(Epoch_hi) - psi_q(Epoch_hi);
                if(Heading_error > pi)
                    Heading_error = Heading_error - 2*pi;
                end
                if(Heading_error < -pi)
                    Heading_error = Heading_error + 2*pi;
                end
                
                HeadingIndex = size(H,1)+1;
                H(HeadingIndex,9) = -1;
                z_Vec(HeadingIndex) = Heading_error;
                R(HeadingIndex,HeadingIndex) = CompassVariance*10;
                
                


                % calculate DOPS 
                H4(1:NumberGPSMeasurements,1:3) = H(1:NumberGPSMeasurements,1:3);
                H4(1:NumberGPSMeasurements,4) = H(1:NumberGPSMeasurements,16);

                 
                H(:,1) = H(:,1)*RM;    
                H(:,2) = H(:,2)*RP*cos(Pos_LLH(Epoch_hi,1)); %I think this is right, Reh*cos(lat), not Rnh*cos(lat)    
                H(:,3) = H(:,3)*-1;
                
                % observability matrix
                F_TOTAL = zeros(NumberStates,NumberStates);
                F_TOTAL(1:9,1:9) = F_INS(1:9,1:9);
                F_TOTAL(10,10) = -x_accel_beta;
                F_TOTAL(11,11) = -y_accel_beta;
                F_TOTAL(12,12) = -z_accel_beta;
                F_TOTAL(13,13) = -x_gyro_beta;
                F_TOTAL(14,14) = -y_gyro_beta;
                F_TOTAL(15,15) = -z_gyro_beta;
                F_TOTAL(16,17) = 1.0;
                
                HSize = size(H,1);
                
                clear ObsMatrix;
                for i=1:NumberStates
                   ObsMatrix((i-1)*HSize+1:i*HSize,1:NumberStates) = H *  phi2^(i-1);
                end
                
                ObsRank(Epoch_lo) = rank(ObsMatrix);    
            
                
                
                phi = eye(NumberStates,NumberStates);
                phi = expm(F_INS);
                phi(7:9,13:15) = -C_BN; %% attitude to gyro
                phi(4:6,10:12) = C_BN; %%  velocity to acceleration
                phi(10,10) = exp( - x_accel_beta );
                phi(11,11) = exp( - y_accel_beta );
                phi(12,12) = exp( - z_accel_beta );
                phi(13,13) = exp( - x_gyro_beta );
                phi(14,14) = exp( - y_gyro_beta );
                phi(15,15) = exp( - z_gyro_beta );
                phi(16,17) = 1.0;
                
                Qd = phi*(G*Q*G')*phi';

                [EKF_GDOP(Epoch_lo) EKF_PDOP(Epoch_lo) EKF_HDOP(Epoch_lo) EKF_VDOP(Epoch_lo) EKF_TDOP(Epoch_lo)] = GARD_CalculateDOPS(H4);

                % calculate the elevations of each SV
                for k = 1:NumberGPSMeasurements
                    Elevation(Epoch_lo,SV_Vec(k)) = asin(H4(k,3));
                end
                
                % save for debugging
                %phi2_save(Epoch_lo,:,:) = phi2;
                %Q2d_save(Epoch_lo,:,:) = Q2d;
                
               


                % evaluate the KF
                %[x_hat_out, P_out] = GARD_EvaluateKF(dt, x_hat_in, P_in, phi, H, z_Vec, Q, R);
                % Prediction step
                %x_hat_minus = phi * x_hat_in;
                %x_hat_minus = x_hat_in;  %%% note that the INS mechanisation already propogates the state vector with a trapezoidal integration
                %P_minus = phi * P_in_full * phi' + Q;
                P_minus = P_in_full;
                
                % correction step

                % calculate the kalman gain
                V = H * P_minus * H' + R;
                K = P_minus * H' * inv(V);

                % update equations
                %v = z - H * x_hat_minus;
                %x_hat_out_full = x_hat_minus + K * (z_Vec);
                x_hat_out_full = K * (z_Vec);
                %x_hat_out_full = x_hat_in +  K * (z_Vec);

                %P_out_full = P_minus - K * H * P_minus;
                P_out_full = (eye(NumberStates,NumberStates) - K*H)*P_minus*(eye(NumberStates,NumberStates)-K*H)' + K*R*K';
    
                P_in_full = P_out_full;
                
                P_out_diag(:,Epoch_hi) = diag(P_out_full);
             
                
                z_save(1:HSize,Epoch_lo) = z_Vec;
                x_save(:,Epoch_lo) = x_hat_out_full;

                %% apply correction
                Vel_NED(Epoch_hi,:) = Vel_NED(Epoch_hi,:) + x_hat_out_full(4:6)';
                Pos_LLH(Epoch_hi,:) = Pos_LLH(Epoch_hi,:) + x_hat_out_full(1:3)';
                
                UserClock(Epoch_hi,1) =  x_hat_out_full(16) + UserClock(Epoch_hi,1);
                UserClock(Epoch_hi,2) = x_hat_out_full(17) + UserClock(Epoch_hi,2);
                %UserClock(Epoch_hi,1) = UserPos_T(Epoch_lo);
                %UserClock(Epoch_hi,2) = UserVel_T(Epoch_lo);
                
                % convert tilt error to DCM Update
                del_alpha = x_hat_out_full(7);
                del_beta = x_hat_out_full(8);
                del_gamma = x_hat_out_full(9);

                % save bias correction
               if(Epoch_lo == 1)
                    AccelBias(Epoch_lo,:) = AccelBias(Epoch_lo,:) - x_hat_out_full(10:12)';   
                    GyroBias(Epoch_lo,:) = GyroBias(Epoch_lo,:) - x_hat_out_full(13:15)';           
                 else
                    AccelBias(Epoch_lo,:) = AccelBias(Epoch_lo-1,:) - x_hat_out_full(10:12)';
                    GyroBias(Epoch_lo,:) = GyroBias(Epoch_lo-1,:) - x_hat_out_full(13:15)';
               end
               
                
                % reset state vector                
                x_hat_out_full(1:9) = 0;
                x_hat_out_full(10:15) = 0;
                x_hat_out_full(16:17) = 0;
                
                x_hat_minus = x_hat_out_full;

                % correct atttitude
                del_att_skew = [0         -del_gamma   del_beta; ...
                               del_gamma  0          -del_alpha; ...
                               -del_beta  del_alpha   0];

                
                %C_BN = (eye(3,3) - del_att_skew) * C_BN; % DCM correction
                C_BN = C_BN * (eye(3,3) + del_att_skew); % DCM correction
                
                C_BN = GARD_OrthogonaliseDCM(C_BN);

%                 q0_INS = q_INS(Epoch_hi,1);
%                 q1_INS = q_INS(Epoch_hi,2);
%                 q2_INS = q_INS(Epoch_hi,3);
%                 q3_INS = q_INS(Epoch_hi,4);
%                 
%                 T_alpha_INS = 0.5*[q1_INS, q2_INS, q3_INS;
%                     -q0_INS -q3_INS q2_INS;
%                     q3_INS, -q0_INS, -q1_INS;
%                     -q2_INS, q1_INS, -q0_INS;];
% 
% 
%                 Quat_errors_INS(Epoch_lo,:) = T_alpha_INS*[ del_alpha,  del_beta,  del_gamma]';                
%                 
%                 q_INS(Epoch_hi,:) = q_INS(Epoch_hi,:) + Quat_errors_INS(Epoch_lo,:);
%                 q_INS(Epoch_hi,:) = q_INS(Epoch_hi,:) / norm(q_INS(Epoch_hi,:));
                
%                 q_INS(Epoch_hi,:) = quatupdate(q_INS(Epoch_hi,:)',-del_alpha,-del_beta,-del_gamma)';

%                C_BN = GARDSim_DCMfromQuat(q_INS(Epoch_hi,:));
                q_INS(Epoch_hi,:) = GARD_DCMToQuat(C_BN);
                
                % calculate euler angles for reference
                phi_q(Epoch_hi) = atan2(C_BN(3,2),C_BN(3,3));
                theta_q(Epoch_hi) = asin(-C_BN(3,1));
                psi_q(Epoch_hi) = atan2(C_BN(2,1),C_BN(1,1));
                
                
                
                x_out = [Pos_LLH(Epoch_hi,:)';Vel_NED(Epoch_hi,:)';q_INS(Epoch_hi,:)'];
                
                
                %%% END FULL SOLUTION
                
                %% start sub solutions
                NumberSubMeasurements = NumberGPSMeasurements-1;
                NumberSubFilters = nchoosek(NumberGPSMeasurements,NumberSubMeasurements);

                R_Sub = eye(2*NumberSubMeasurements);
                for k = 1:NumberSubMeasurements
                    R_Sub(k,k) = RangeNoiseVariance;
                    R_Sub(k+NumberGPSMeasurements,k+NumberGPSMeasurements) = RangeRateNoiseVariance;
                end
                
                % initialise sub filters
                if Epoch_lo == StartTime
                    for i=1:NumberSubFilters
                        x_hat_minus_Sub(:,i) = x_hat_minus;
                        P_minus_Sub(:,:,i) = P_minus;
                        
                        
                        
                        EKF_Px_dual_minus(:,:,i) = [P_out_full, zeros(NumberStates,NumberStates);
                                        zeros(NumberStates,NumberStates),zeros(NumberStates,NumberStates)];
                    end
                end
                 
                for SubSolution = 1:NumberSubFilters
                     
                    
                    %% propagate sub-filter variance (important for HPL)
                    %x_hat_minus_Sub(:,SubSolution) = phi * x_hat_minus_Sub(:,SubSolution);
                    x_hat_minus_Sub(:,SubSolution) = zeros(NumberStates,1);
                    P_minus_Sub(:,:,SubSolution) = phi * P_minus_Sub(:,:,SubSolution) * phi' + Qd;
                    
%                     PR_SubVec(1:SubSolution-1) = PR_Vec(1:SubSolution-1);
%                     PR_SubVec(SubSolution:NumberSubMeasurements) = PR_Vec(SubSolution+1:NumberGPSMeasurements);
% 
%                     PRR_SubVec(1:SubSolution-1) = PRR_Vec(1:SubSolution-1);
%                     PRR_SubVec(SubSolution:NumberSubMeasurements) = PRR_Vec(SubSolution+1:NumberGPSMeasurements);
% 
%                     
%                     
%                     % formulate the measurment vector
%                     EKF_y_k_Sub = [PR_SubVec(1:NumberSubMeasurements)';PRR_SubVec(1:NumberSubMeasurements)'];
% 
% 
%                     PR_Vec_minus_Sub(1:SubSolution-1) = PR_Vec_minus(1:SubSolution-1);
%                     PR_Vec_minus_Sub(SubSolution:NumberSubMeasurements) = PR_Vec_minus(SubSolution+1:NumberGPSMeasurements);
% 
%                     PRR_Vec_minus_Sub(1:SubSolution-1) = PRR_Vec_minus(1:SubSolution-1);
%                     PRR_Vec_minus_Sub(SubSolution:NumberSubMeasurements) = PRR_Vec_minus(SubSolution+1:NumberGPSMeasurements);
% 
% 
%                     ys_kminus_Sub = [PR_Vec_minus_Sub'; PRR_Vec_minus_Sub'];
%                     
%                     z_Vec_Sub = (EKF_y_k_Sub - ys_kminus_Sub);  
%                     
% 
%                 
%                     EKF_H_Sub(1:SubSolution-1,:) = H(1:SubSolution-1,:);
%                     EKF_H_Sub(SubSolution:NumberSubMeasurements,:) = H(SubSolution+1:NumberGPSMeasurements,:);
%                     EKF_H_Sub(NumberSubMeasurements+1:NumberSubMeasurements+SubSolution-1,:) = ...
%                             H(NumberSubMeasurements+1:NumberSubMeasurements+SubSolution-1,:);
%                     EKF_H_Sub(NumberSubMeasurements+SubSolution:NumberSubMeasurements+NumberSubMeasurements,:) = ...
%                             H(NumberGPSMeasurements+SubSolution+1:NumberGPSMeasurements+NumberGPSMeasurements,:);
%                         
%                     EKF_H_Sub(2*NumberSubMeasurements+1,:) = H(2*NumberGPSMeasurements+1,:);
%                     
%                     
%                     
%                     
%                     HeadingIndex = size(EKF_H_Sub,1);
%                     z_Vec_Sub(HeadingIndex) = Heading_error;
%                     R_Sub(HeadingIndex,HeadingIndex) = CompassVariance;
%                 
%                     % calculate the kalman gain
%                     V_Sub = EKF_H_Sub * P_minus_Sub(:,:,SubSolution) * EKF_H_Sub' + R_Sub;
%                     K_Sub = P_minus_Sub(:,:,SubSolution) * EKF_H_Sub' * inv(V_Sub);
% 
%                     % update equations
%                     EKF_x_hat_kplus_Sub(:,SubSolution) = x_hat_minus_Sub(:,SubSolution)  + K_Sub * (z_Vec_Sub);
% 
%                     EKF_Px_kplus_Sub(:,:,SubSolution) = (eye(NumberStates,NumberStates) - K_Sub*EKF_H_Sub) * P_minus_Sub(:,:,SubSolution) * ...
%                                 (eye(NumberStates,NumberStates)-K_Sub*EKF_H_Sub)'...
%                                 + K_Sub*R_Sub*K_Sub';


                    EKF_H_Sub = H;
                    EKF_H_Sub(SubSolution,:) = zeros(1,NumberStates);
                    EKF_H_Sub(NumberGPSMeasurements+SubSolution,:) = zeros(1,NumberStates);

                    R_Sub = R;
                    %R_Sub(SubSolution,:) = zeros(1,NumberStates);
                    %R_Sub(NumberGPSMeasurements+SubSolution,:) = zeros(1,NumberStates);
                    
                    
                    % calculate the kalman gain
                    V_Sub = EKF_H_Sub * P_minus_Sub(:,:,SubSolution) * EKF_H_Sub' + R_Sub;
                    K_Sub = P_minus_Sub(:,:,SubSolution) * EKF_H_Sub' * inv(V_Sub);

                    % update equations
                    EKF_x_hat_kplus_Sub(:,SubSolution) = x_hat_minus_Sub(:,SubSolution)  + K_Sub * (z_Vec);

                    EKF_Px_kplus_Sub(:,:,SubSolution) = (eye(NumberStates,NumberStates) - K_Sub*EKF_H_Sub) * P_minus_Sub(:,:,SubSolution) * ...
                                (eye(NumberStates,NumberStates)-K_Sub*EKF_H_Sub)'...
                                + K_Sub*R_Sub*K_Sub';

                    % save result
                    x_hat_minus_Sub(:,SubSolution) = EKF_x_hat_kplus_Sub(:,SubSolution);
                    % update for next epoch
                    P_minus_Sub(:,:,SubSolution) = EKF_Px_kplus_Sub(:,:,SubSolution);
                    
                    
                    % propagate dual covariance 
                    EKF_Px_dual_plus(:,:,SubSolution) = GARD_DualCovariancePropagator(EKF_Px_dual_minus(:,:,SubSolution),phi,Qd,K,K_Sub,H',R,SubSolution);
                    EKF_Px_dual_save(Epoch_lo,:,SubSolution) = diag(EKF_Px_dual_plus(:,:,SubSolution));
                end % end sub solutions

                
                
                EKF_Px_dual_minus = EKF_Px_dual_plus;


                for j = 1:NumberSubFilters
                    % note that only the position estimates are used - not
                    % clock or velocity
                    P_ned = P_out_full(1:3,1:3);
                    P_ned_Sub = EKF_Px_kplus_Sub(1:3,1:3,j);

                    NED_ss = (x_hat_out_full(1:3) - EKF_x_hat_kplus_Sub(1:3,j));
                    
                    EKF_Beta_ss_H(:,j) = NED_ss(1:2);
                    EKF_B_ss_H(:,:,j) = (P_ned_Sub(1:2,1:2) - P_ned(1:2,1:2));

                    EKF_lambda_ss_H(Epoch_lo,j) = (EKF_Beta_ss_H(:,j)' * EKF_Beta_ss_H(:,j));

                    EKF_B_lambda_H = eigs(EKF_B_ss_H(:,:,j));
                    EKF_TD_H(Epoch_lo,j) = sqrt(max(EKF_B_lambda_H)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));

                    if(EKF_lambda_ss_H(Epoch_lo,j) > EKF_TD_H(Epoch_lo,j))
                        disp(sprintf('%s H-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
                    end

                    EKF_Beta_ss_V(:,j) = NED_ss(3);
                    EKF_B_ss_V(:,:,j) = abs(P_ned_Sub(3,3) - P_ned(3,3));

                    EKF_lambda_ss_V(Epoch_lo,j) = (EKF_Beta_ss_V(:,j)' * EKF_Beta_ss_V(:,j));

                    EKF_B_lambda_V = eigs(EKF_B_ss_V(:,:,j));
                    EKF_TD_V(Epoch_lo,j) = sqrt(max(EKF_B_lambda_V)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));

                    if(EKF_lambda_ss_V(Epoch_lo,j) > EKF_TD_V(Epoch_lo,j))
                        disp(sprintf('%s V-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
                    end

                    %% this HPL represents the fault-free (H0) hypothesis. 
                    Sigma_MD_H = abs(norminv(P_MD_H,0,1));
                    Sigma_MD_V = abs(norminv(P_MD_V,0,1));
                    EKF_HPL_sub(Epoch_lo,j) = Sigma_MD_H * sqrt(P_ned(1,1) + P_ned(2,2)) + EKF_TD_H(Epoch_lo,j);
                    EKF_VPL_sub(Epoch_lo,j) = Sigma_MD_V * sqrt(P_ned(3,3)) + EKF_TD_V(Epoch_lo,j);
                end

            
                clear z_Vec_Sub EKF_H_Sub R_Sub K_Sub V_Sub ys_kminus_Sub EKF_y_k_Sub PR_Vec_minus_Sub PRR_Vec_minus_Sub PR_SubVec PRR_SubVec; 
                
            EKF_HPL(Epoch_lo) = max(EKF_HPL_sub(Epoch_lo,:));
            EKF_VPL(Epoch_lo) = max(EKF_VPL_sub(Epoch_lo,:));
                
        %% calculate Least Squares solution for comparisson
        [LSQ_SolutionVec(Epoch_lo,:), LSQ_VarSolutionVec(Epoch_lo,:), LSQ_NumIterations(Epoch_lo,:), ...
        LSQ_ResidualVector(Epoch_lo,1:NumberGPSMeasurements), LSQ_M, LSQ_Fail(Epoch_lo), LSQ_limit(Epoch_lo), ...
        LSQ_DOP(Epoch_lo,:)] = GARD_LSQ(UserPos,NumberGPSMeasurements,PR_Vec(1:NumberGPSMeasurements),SVPos);

        [LSQVel_SolutionVec(Epoch_lo,:), LSQVel_VarSolutionVec(Epoch_lo,:), LSQVel_NumIterations(Epoch_lo,:),...
            LSQVel_ResidualVector(Epoch_lo,1:NumberGPSMeasurements), LSQVel_M, LSQVel_Fail(Epoch_lo,:), ...
            LSQVel_limit(Epoch_lo,:)] = GARD_LSQVel(UserPos,UserVel,NumberGPSMeasurements,PRR_Vec(1:NumberGPSMeasurements),SVPos, SVVel);
    
        LSQVel_NED(Epoch_lo,:) = (Tecef2ned*LSQVel_SolutionVec(Epoch_lo,1:3)');
        
        [LSQ_LLH(Epoch_lo,1) LSQ_LLH(Epoch_lo,2) LSQ_LLH(Epoch_lo,3)] = ECEF2LLH(LSQ_SolutionVec(Epoch_lo,1:3));

        LSQ_LLH_error(Epoch_lo,:) = LSQ_LLH(Epoch_lo,:) - pos_truth_llh(2:4,Epoch_hi)';
        LSQVel_error(Epoch_lo,:) = LSQVel_NED(Epoch_lo,:) - vel_truth(2:4,Epoch_hi)';

if Sim_Data == 0
        LSQ_Clock_error(Epoch_lo,1) = LSQ_SolutionVec(Epoch_lo,4) - UserPos_T(Epoch_lo);
else
        LSQ_Clock_error(Epoch_lo,1) = LSQ_SolutionVec(Epoch_lo,4);
end

        
        LSQ_Clock_error(Epoch_lo,2) = LSQVel_SolutionVec(Epoch_lo,4);

if Sim_Data == 0
        Clock_error(Epoch_lo) = UserClock(Epoch_hi,1) - UserPos_T(Epoch_lo);
        Freq_error(Epoch_lo) = UserClock(Epoch_hi,2) - UserVel_T(Epoch_lo);
else
        Clock_error(Epoch_lo) = UserClock(Epoch_hi,1);
        Freq_error(Epoch_lo) = UserClock(Epoch_hi,2);
end
        
        %end % this is the end of the if statement which decides if this is the first epoch or not
        
        

        
    end % end low-speed loop
    
    hgt = Pos_LLH(Epoch_hi,3);
    Lat_error(Epoch_hi) = (Pos_LLH(Epoch_hi,1) - pos_truth_llh(2,Epoch_hi)) * (RM+hgt);
    Long_error(Epoch_hi) = (Pos_LLH(Epoch_hi,2) - pos_truth_llh(3,Epoch_hi)) * ((RP+hgt) * cos(Pos_LLH(Epoch_hi,1)));
    Height_error(Epoch_hi) = (Pos_LLH(Epoch_hi,3) - pos_truth_llh(4,Epoch_hi));

    Vel_error(Epoch_hi,1:3) = Vel_NED(Epoch_hi,1:3) - vel_truth(2:4,Epoch_hi)';
    
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
    
    if(mod(Epoch_hi,1000) == 0)
        disp(strcat(ModuleName,sprintf('Completed Epoch %d',Epoch_hi)));
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
        save(ekf_results_filename);
    end
    
end % end high-speed loop



% save results
save(ekf_results_filename);
logline = sprintf('%s Saving Results to %s\n',ModuleName,ekf_results_filename);
fid = fopen(progress_filename,'a');
fprintf(fid,logline);
fclose(fid);
        
disp(logline);
disp(strcat(ModuleName,'Done!'));







%% plot results
% 

% % pos
% figure();
% plot(pos_truth_llh(3,:)*180/pi,pos_truth_llh(2,:)*180/pi,'g')
% hold on; grid on;
% plot(Pos_LLH(:,2)*180/pi,Pos_LLH(:,1)*180/pi,'r')
% hold off;
% xlabel('Longitude');
% ylabel('Latitude');

figure();
hold on;
plot(TimeINS,(Lat_error),'b');
grid on;
hold on;
plot(TimeINS,2*sqrt(Pos_var(:,1)),'g');
plot(TimeINS,-2*sqrt(Pos_var(:,1)),'g');
%plot(TimeGPS(1:Epoch_lo),HPL,'r');
xlabel('Simulation Time (sec)');
ylabel('North Error (m)');
legend('North Error','2-\sigma Error Bound');
axis([TimeINS(1) TimeINS(NumberEpochsINS) -20 20]);
%axis([0 3000 0 50]);
hold off;

figure();
hold on;
plot(TimeINS,(Long_error),'b');
grid on;
hold on;
plot(TimeINS,2*sqrt(Pos_var(:,2)),'g');
plot(TimeINS,-2*sqrt(Pos_var(:,2)),'g');
%plot(TimeGPS(1:Epoch_lo),HPL,'r');
xlabel('Simulation Time (sec)');
ylabel('East Error (m)');
legend('East Error','2-\sigma Error Bound');
axis([TimeINS(1) TimeINS(NumberEpochsINS) -20 20]);
%axis([0 3000 0 50]);
hold off;




%% vertical error with VPL
figure();
plot(TimeINS,(Height_error),'b');
hold on;
plot(TimeINS,2*sqrt(Pos_var(:,3)),'g');
plot(TimeINS,-2*sqrt(Pos_var(:,3)),'g');
%plot(TimeGPS(1:Epoch_lo),VPL,'r');
xlabel('Simulation Time (sec)');
ylabel('Vertical Error (m)');
legend('Vertical Error','2-\sigma Error Bound','VPL');
%axis([TimeINS(1) TimeINS(NumberEpochsINS) 0 20]);
%axis([0 3000 0 50]);
hold off;
grid on;

%% plot satellite usage
% figure();
% plot(TimeGPS,SatellitesUsed);
% xlabel('Test Time (sec)');
% ylabel('Satellites Used');



%% attitude errors
figure();
subplot(3,1,1),plot(TimeINS,phi_err*r2d,'r');grid on;hold on;ylabel('Roll (deg)');
subplot(3,1,1),plot(TimeINS,2*sqrt(P_out_diag(7,:)*180/pi),'g');
subplot(3,1,1),plot(TimeINS,-2*sqrt(P_out_diag(7,:)*180/pi),'g');
%axis([TimeINS(1) TimeINS(NumberEpochsINS) -5 5]);
subplot(3,1,2),plot(TimeINS,theta_err*r2d,'r');grid on;hold on;ylabel('Pitch (deg)');
subplot(3,1,2),plot(TimeINS,2*sqrt(P_out_diag(8,:)*180/pi),'g');
subplot(3,1,2),plot(TimeINS,-2*sqrt(P_out_diag(8,:)*180/pi),'g');
%axis([TimeINS(1) TimeINS(NumberEpochsINS) -5 5]);
subplot(3,1,3),plot(TimeINS,psi_err*r2d,'r');grid on;hold on;ylabel('Yaw (deg)');
subplot(3,1,3),plot(TimeINS,2*sqrt(P_out_diag(9,:)*180/pi),'g');
subplot(3,1,3),plot(TimeINS,-2*sqrt(P_out_diag(9,:)*180/pi),'g');
%axis([TimeINS(1) TimeINS(NumberEpochsINS) -10 10]);
xlabel('Time (Sec)');

%legend('Roll','Pitch','Yaw');
%axis([TimeINS(1) TimeINS(NumberEpochsINS) -5 5]);

% figure();
% plot(lamda_ss_out);
% grid on;
% hold on;
% plot(TD);
% axis([0 TimeGPS(Epoch_lo) -2 120]);
% xlabel('Simulation Time (Seconds)');
% ylabel('Test Statistic');
% title('Solution Separation Statistics - PRN1 15m 1m/s ramp fault Fault at 500s');


%% plot velocity error
figure(); 
subplot(3,1,1),plot(TimeINS,Vel_error(:,1));hold on; grid on;ylabel('North (m/s)');title('Velocity Error');
subplot(3,1,1),plot(TimeINS,sqrt(P_out_diag(4,:)),'g');
subplot(3,1,1),plot(TimeINS,-sqrt(P_out_diag(4,:)),'g');
subplot(3,1,2),plot(TimeINS,Vel_error(:,2));hold on; grid on;ylabel('East (m/s)');
subplot(3,1,2),plot(TimeINS,sqrt(P_out_diag(5,:)),'g');
subplot(3,1,2),plot(TimeINS,-sqrt(P_out_diag(5,:)),'g');
subplot(3,1,3),plot(TimeINS,Vel_error(:,3));hold on; grid on;ylabel('Down (m/s)');
subplot(3,1,3),plot(TimeINS,sqrt(P_out_diag(6,:)),'g');
subplot(3,1,3),plot(TimeINS,-sqrt(P_out_diag(6,:)),'g');
xlabel('Simulation Time (seconds)');
%axis([TimeINS(1) TimeINS(Epoch_hi) -5 5]);

%% plot gyro and accelerometer bias
figure(); 
subplot(3,1,1),plot(TimeGPS,(GyroBias(:,1) - GyroBiasTruth(1))*r2d); grid on; hold on; ylabel('Roll (deg/sec)');title('Gyroscope Bias Estimate Error');
subplot(3,1,1),plot(TimeINS,2*sqrt(P_out_diag(13,:)*r2d),'g');
subplot(3,1,1),plot(TimeINS,-2*sqrt(P_out_diag(13,:)*r2d),'g');
subplot(3,1,2);plot(TimeGPS,(GyroBias(:,2) - GyroBiasTruth(2))*r2d); grid on; hold on; ylabel('Pitch (deg/sec)');
subplot(3,1,2),plot(TimeINS,2*sqrt(P_out_diag(14,:)*r2d),'g');
subplot(3,1,2),plot(TimeINS,-2*sqrt(P_out_diag(14,:)*r2d),'g');
subplot(3,1,3);plot(TimeGPS,(GyroBias(:,3) - GyroBiasTruth(3))*r2d); grid on; hold on; ylabel('Yaw (deg/sec)');
subplot(3,1,3),plot(TimeINS,2*sqrt(P_out_diag(15,:)*r2d),'g');
subplot(3,1,3),plot(TimeINS,-2*sqrt(P_out_diag(15,:)*r2d),'g');
xlabel('Simulation Time (sec)');

figure(); 
subplot(3,1,1),plot(TimeGPS,AccelBias(:,1) - AccelBiasTruth(1)); grid on; hold on; ylabel('X (m/s/s)');title('Accelerometer Bias Estimate Error');
subplot(3,1,1),plot(TimeINS,2*sqrt(P_out_diag(10,:)),'g');
subplot(3,1,1),plot(TimeINS,-2*sqrt(P_out_diag(10,:)),'g');
subplot(3,1,2);plot(TimeGPS,AccelBias(:,2) - AccelBiasTruth(2)); grid on; hold on; ylabel('Y (m/s/s)');
subplot(3,1,2),plot(TimeINS,2*sqrt(P_out_diag(11,:)),'g');
subplot(3,1,2),plot(TimeINS,-2*sqrt(P_out_diag(11,:)),'g');
subplot(3,1,3);plot(TimeGPS,AccelBias(:,3) - AccelBiasTruth(3)); grid on; hold on; ylabel('Z (m/s/s)');
subplot(3,1,3),plot(TimeINS,2*sqrt(P_out_diag(12,:)),'g');
subplot(3,1,3),plot(TimeINS,-2*sqrt(P_out_diag(12,:)),'g');
xlabel('Simulation Time (sec)');

% figure();hold on; 
% subplot(2,1,1);plot(TimeGPS,GyroBias); grid on; ylabel('Gyro bias');
% subplot(2,1,2);plot(TimeGPS,AccelBias); grid on; ylabel('Accel bias'); 
% xlabel('Simulation Time (sec)');
% hold off;

%% plot clock
% figure(); hold on;
% plot(TimeINS,UserClock(:,1));
% plot(TimeINS,UserClock(:,2));
% plot(TimeGPS,UserPos_T,'r');
% plot(LSQ_SolutionVec(:,4),'r*');
% grid on;
% xlabel('Simulation Time (sec)');
% legend('Clock Estimate','Clock Truth');
figure();
plot(TimeGPS,Clock_error); hold on; grid on;
plot(TimeINS,2*sqrt(P_out_diag(16,:)),'g');
plot(TimeINS,-2*sqrt(P_out_diag(16,:)),'g');
axis([0 TimeGPS(Epoch_lo) -15 15]);
xlabel('Simulation Time (sec)');
ylabel('Error (metres)');
title('Clock Estimation Error');

figure();
plot(TimeGPS,Freq_error); hold on; grid on;
plot(TimeINS,2*sqrt(P_out_diag(17,:)),'g');
plot(TimeINS,-2*sqrt(P_out_diag(17,:)),'g');
axis([StartTime StopTime -1 1]);
xlabel('Simulation Time (sec)');
ylabel('Error (metres/sec)');
title('Frequency Estimation Error');



%% 
% figure(); grid on; hold on;
% plot(LSQVel_NED,'*');
% plot(vel_truth(1,:),vel_truth(2:4,:))
% plot(TimeINS,Vel_NED);

%% calculate RMS errors
Height_RMS = sqrt(sum(Height_error.^2)/NumberEpochsINS)
Height_RMS_sub = sqrt(sum(Height_error(3000:NumberEpochsINS).^2)/(NumberEpochsINS-3000))

Lat_RMS = sqrt(sum(Lat_error.^2)/NumberEpochsINS)
Lat_RMS_sub = sqrt(sum(Lat_error(3000:NumberEpochsINS).^2)/(NumberEpochsINS-3000))

Long_RMS = sqrt(sum(Long_error.^2)/NumberEpochsINS)
Long_RMS_sub = sqrt(sum(Long_error(3000:NumberEpochsINS).^2)/(NumberEpochsINS-3000))




figure();hold on;grid on;
plot(TimeINS(1:Epoch_hi),sqrt(Lat_error.^2+Long_error.^2));
plot(TimeGPS(1:Epoch_lo),EKF_HPL*sqrt(RMh*RPh),'r','LineWidth',2);
legend('Horizontal Error','HPL');
xlabel('Simulation Time (sec)');
ylabel('Error (metres)');


figure();hold on;grid on;
plot(TimeINS(1:Epoch_hi),sqrt(Height_error.^2));
plot(TimeGPS(1:Epoch_lo),EKF_VPL,'r','LineWidth',2);
legend('Vertical Error','VPL');
xlabel('Simulation Time (sec)');
ylabel('Error (metres)');

figure();
hold on; grid on;
plot(TimeGPS,EKF_HPL*sqrt(RMh*RPh),'LineWidth',2);
plot(TimeGPS,EKF_VPL,'r','LineWidth',2);
legend('HPL','VPL');
axis([StartTime StopTime 0 20]);
xlabel('Simulation Time (sec)');
ylabel('Protection Level (m)');
title('EKF GPS-INS Protection Level - RNAV Approach');

% plot lat error vs lsq
figure();
hold on; grid on;
plot(TimeINS,abs(Lat_error),'b');
plot(TimeGPS,abs(LSQ_LLH_error(:,1))*RM,'r');
legend('UKF INS','LSQ');
axis([StartTime StopTime 0 1]);
title('Latitude Error');

figure();
hold on; grid on;
plot(TimeINS,abs(Long_error),'b');
plot(TimeGPS,abs(LSQ_LLH_error(:,2))*RP*cos(-0.47),'r');
legend('UKF INS','LSQ');
axis([StartTime StopTime 0 1]);
title('Longitude Error');

figure();
hold on; grid on;
plot(TimeINS,Height_error,'b');
plot(TimeGPS,(LSQ_LLH_error(:,3)),'r');
legend('EKF INS','LSQ');
%axis([StartTime StopTime 0 2]);
title('Height Error');

