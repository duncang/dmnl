% GPS-INS position solution using the UKF
% Written by Duncan Greer 29 May 2007
% $Id: GARDSim_GPSINS_TightlyCoupledUKF.m 1879 2008-07-15 05:20:21Z n2523710 $
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


ModuleName = '[GARDSim_GPSINS_TightlyCoupledUKF] ';


% setup constants
WGS84Constants;
d2r = pi/180;
r2d = 180/pi;
g = -9.79;  % m/s/s
GPSConstants;

NumberStates = 17;


% load data path
DataPath = 'data/rnav_approach/';

StartTime = 000; % start at 100 seconds into the data set
StopTime = 3000;

disp(strcat(ModuleName,'Loading data from: "',DataPath,'"'));

% load the GPS and INS measurements
if(~exist('PR_Sim','var'))
    GPSMeasurementData = strcat(DataPath,'PR_Simulation.mat');
    load(GPSMeasurementData);
end

% load IMU measurements data
UseNoisy = 0;

GRASOn = 0;  % if GRAS is enabled then there is no iono bias

if(~exist('sensors','var'))

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
if(~exist('gravity','var'))
    load(strcat(DataPath,'gravity.mat'));
end

% load truth data
if(~exist('pos_truth_llh','var'))
    load(strcat(DataPath,'pos_truth_llh.mat'));
    load(strcat(DataPath,'vel_truth.mat'));
    load(strcat(DataPath,'att_truth.mat'));
end

disp(strcat(ModuleName,'Finished loading data'));

gps_dt = 1;
ins_dt = 0.01;

% initialise epochs
Epoch_lo = (StartTime) / gps_dt;
Epoch_hi = 0;

% position
%InitialPosition = [-27.4*pi/180 153.2*pi/180 4500/3.28]';  % somewhere near brisbane in LLH, rads and meters
InitialPosition = [pos_truth_llh(2,StartTime/ins_dt+1),pos_truth_llh(3,StartTime/ins_dt+1),pos_truth_llh(4,StartTime/ins_dt+1)]';
% velocity
%InitialVelocity = [0 -60 0]';  % NED velocities in m/s
InitialVelocity = [vel_truth(2,StartTime/ins_dt+1),vel_truth(3,StartTime/ins_dt+1),vel_truth(4,StartTime/ins_dt+1)]';
% attitude
%InitialAttitude = [0 0 270*pi/180]';  % roll, pitch, yaw in radians
InitialAttitude = [att_truth(2,StartTime/ins_dt+1),att_truth(3,StartTime/ins_dt+1),att_truth(4,StartTime/ins_dt+1)]';

C_BN = GARDSim_DCMfromEuler(InitialAttitude(1),InitialAttitude(2),InitialAttitude(3))
% 

%UserPos = [InitialPosition',0];
%UserVel = [InitialVelocity',0];


% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));


NumberEpochsGPS=min(length(PR_Sim),StopTime/gps_dt+1);
NumberSVs = 32;

SVDontUse = zeros(1,NumberSVs);

SVDontUse(11) = 1;

%% artificially inflate the PR noise at a certain time
NoiseOnTime = StartTime + 500;
NoiseOffTime = StopTime;
NoiseAmplification = 1;

%PR_Sim(NoiseOnTime:NoiseOffTime,:) = PR_Sim(NoiseOnTime:NoiseOffTime,:) + randn(NoiseOffTime-NoiseOnTime+1,NumberSVs) * NoiseAmplification;


LOCAL_GRAVITY = 9.80;

TimeINS = 0:ins_dt:NumberEpochsGPS-1;
TimeGPS = 0:gps_dt:NumberEpochsGPS-1;

NumberEpochsINS = min(length(TimeINS),StopTime/ins_dt + 1);

% initialise storage vars to svae time
Lat_error = zeros(1,NumberEpochsINS);
Long_error = zeros(1,NumberEpochsINS);
Height_error = zeros(1,NumberEpochsINS);
% P_xn = zeros(1,NumberEpochsINS);
% P_yn = zeros(1,NumberEpochsINS);
% P_zn = zeros(1,NumberEpochsINS);
% V_xn = zeros(1,NumberEpochsINS);
% V_yn = zeros(1,NumberEpochsINS);
% V_zn = zeros(1,NumberEpochsINS);

Pos_LLH = zeros(NumberEpochsINS,3);
Vel_NED = zeros(NumberEpochsINS,3);
Acc_NED = zeros(NumberEpochsINS,3);

UserClock = [0 0];

% a = zeros(1,NumberEpochsINS);
% b = zeros(1,NumberEpochsINS);
% c = zeros(1,NumberEpochsINS);
% d = zeros(1,NumberEpochsINS);
% 
% a_dot = zeros(1,NumberEpochsINS);
% b_dot = zeros(1,NumberEpochsINS);
% c_dot = zeros(1,NumberEpochsINS);
% d_dot = zeros(1,NumberEpochsINS);


phi_q = zeros(NumberEpochsINS,1);
theta_q = zeros(NumberEpochsINS,1);
psi_q = zeros(NumberEpochsINS,1);

phi_uc = zeros(NumberEpochsINS,1);
theta_uc = zeros(NumberEpochsINS,1);
psi_uc = zeros(NumberEpochsINS,1);

gvec_phi = zeros(NumberEpochsINS,1);
gvec_theta = zeros(NumberEpochsINS,1);

% A_xn = zeros(1,NumberEpochsINS);
% A_yn = zeros(1,NumberEpochsINS);
% A_zn = zeros(1,NumberEpochsINS);
%Acc_NED = zeros(NumberEpochsINS,3);

lat_dot = zeros(NumberEpochsINS,1);
long_dot = zeros(NumberEpochsINS,1);
h_dot = zeros(NumberEpochsINS,1);


Coriolis = zeros(NumberEpochsINS,3);

P_out_diag = zeros(NumberStates,NumberEpochsINS);
Pos_var = zeros(NumberEpochsINS,3);

SatellitesUsed = zeros(NumberEpochsGPS,1);

% setup error values
if(GRASOn == 1)
    PositionError = 5.0;
    VelocityError = 2.0;
    ClockBiasError = 1E5;
    ClockDriftError = 10.0;
    RangeNoiseVariance = 3.0^2; % m
    RangeRateNoiseVariance = 2.0^2; % m/s
else
    PositionError = 5.0;
    VelocityError = 2.0;
    ClockBiasError = 1E5;
    ClockDriftError = 10.0;
    RangeNoiseVariance = 5.0^2; % m
    RangeRateNoiseVariance = 0.5^2; % m/s
end



% setup state transition - time invariant if dt is constant
% power spectral densitites of the process noise - 
x_gyro_beta = 1/100;
y_gyro_beta = 1/100;
z_gyro_beta = 1/100;
x_accel_beta = 1/100;
y_accel_beta = 1/100;
z_accel_beta = 1/100;
x_gyro_Q = 0.05; % PSD of gyro noise
y_gyro_Q = 0.05; % PSD of gyro noise
z_gyro_Q = 0.05; % PSD of gyro noise
x_accel_Q = 0.01; % PSD of accelerometer noise
y_accel_Q = 0.01; % PSD of accelerometer noise
z_accel_Q = 0.01; % PSD of accelerometer noise
% phi = eye(NumberStates,NumberStates);
% phi(1,5) = gps_dt;
% phi(2,6) = gps_dt;
% phi(3,7) = gps_dt;
% phi(4,8) = gps_dt;
% 
% phi2 = eye(NumberStates,NumberStates);
% phi2(1,5) = ins_dt;
% phi2(2,6) = ins_dt;
% phi2(3,7) = ins_dt;
% phi2(4,8) = ins_dt;

% generate initial PHI matrix
F = GARD_GenerateINSFMatrix(InitialPosition,InitialVelocity,[0 0 LOCAL_GRAVITY]);
phi = eye(NumberStates,NumberStates);
phi(1:9,1:9) = phi(1:9,1:9) + F * gps_dt;
% % augmented state matrix
phi(7:9,10:12) = -C_BN * gps_dt;
phi(4:6,13:15) = C_BN * gps_dt;

phi(10,10) = 1 - x_accel_beta * gps_dt;
phi(11,11) = 1 - y_accel_beta * gps_dt;
phi(12,12) = 1 - z_accel_beta * gps_dt;
phi(13,13) = 1 - x_gyro_beta * gps_dt;
phi(14,14) = 1 - y_gyro_beta * gps_dt;
phi(15,15) = 1 - z_gyro_beta * gps_dt;
phi(16,17) = 1.0*gps_dt;

phi2 = eye(NumberStates,NumberStates);
phi2(1:9,1:9) = phi2(1:9,1:9) + F * ins_dt;
% % augmented state matrix
phi2(7:9,10:12) = -C_BN * ins_dt;
phi2(4:6,13:15) = C_BN * ins_dt;

phi2(10,10) = 1 - x_accel_beta * ins_dt;
phi2(11,11) = 1 - y_accel_beta * ins_dt;
phi2(12,12) = 1 - z_accel_beta * ins_dt;
phi2(13,13) = 1 - x_gyro_beta * ins_dt;
phi2(14,14) = 1 - y_gyro_beta * ins_dt;
phi2(15,15) = 1 - z_gyro_beta * ins_dt;
phi2(16,17) = 1.0*ins_dt;



Sf = 8*pi^2*2e-20*Speedoflight^2;  % frequency noise psd
Sc = 2*2e-19*Speedoflight^2;  % bias noise psd

Q = zeros(NumberStates,NumberStates);

Q(10,10) = 2*x_gyro_beta*x_gyro_Q;
Q(11,11) = 2*y_gyro_beta*y_gyro_Q;
Q(12,12) = 2*z_gyro_beta*z_gyro_Q;
Q(13,13) = 2*x_accel_beta*x_accel_Q;
Q(14,14) = 2*y_accel_beta*y_accel_Q;
Q(15,15) = 2*z_accel_beta*z_accel_Q;

Q(16,16) = Sc * gps_dt + Sf * (gps_dt^2)/2;
Q(16,17) = Sc * (gps_dt^2) / 2;
Q(17,16) = Sc * (gps_dt^2) / 2;
Q(17,17) = Sf * gps_dt;

G = zeros(NumberStates,NumberStates);
G(7:9,10:12) = -C_BN;
G(4:6,13:15) = C_BN;
G(10:12,10:12) = eye(3,3);
G(13:15,13:15) = eye(3,3);

% convert to discrete time form
Qd = phi * (G * Q * G') * phi';
Q2d = phi2 * (G * Q * G') * phi2';

% Sp_dt3 = Sp * (gps_dt ^ 3) / 3;
% Sp_dt2 = Sp * (gps_dt ^ 2) / 2;
% Sp_dt = Sp * gps_dt;
% Sf_dt = Sf * gps_dt;
% Sg_dt3 = Sg * (gps_dt ^ 3) / 3;
% Sg_dt2 = Sg * (gps_dt ^ 2) / 2;
% Sg_dt = Sg * gps_dt;
% 
% Q(1,1) = Sp_dt3;
% Q(2,2) = Sp_dt3;
% Q(3,3) = Sp_dt3;
% 
% Q(1,5) = Sp_dt2;
% Q(2,6) = Sp_dt2;
% Q(3,7) = Sp_dt2;
% Q(4,8) = Sg_dt2;
% Q(5,1) = Sp_dt2;
% Q(6,2) = Sp_dt2;
% Q(7,3) = Sp_dt2;
% Q(8,4) = Sg_dt2;
% 
% Q(4,4) = Sf_dt + Sg_dt3;
% Q(5,5) = Sp_dt;
% Q(6,6) = Sp_dt;
% Q(7,7) = Sp_dt;
% Q(8,8) = Sg_dt;



% Sp_dt3 = Sp * (ins_dt ^ 3) / 3;
% Sp_dt2 = Sp * (ins_dt ^ 2) / 2;
% Sp_dt = Sp * ins_dt;
% Sf_dt = Sf * ins_dt;
% Sg_dt3 = Sg * (ins_dt ^ 3) / 3;
% Sg_dt2 = Sg * (ins_dt ^ 2) / 2;
% Sg_dt = Sg * ins_dt;
% 
% Q2(1,1) = Sp_dt3;
% Q2(2,2) = Sp_dt3;
% Q2(3,3) = Sp_dt3;
% 
% Q2(1,5) = Sp_dt2;
% Q2(2,6) = Sp_dt2;
% Q2(3,7) = Sp_dt2;
% Q2(4,8) = Sg_dt2;
% Q2(5,1) = Sp_dt2;
% Q2(6,2) = Sp_dt2;
% Q2(7,3) = Sp_dt2;
% Q2(8,4) = Sg_dt2;
% 
% Q2(4,4) = Sf_dt + Sg_dt3;
% Q2(5,5) = Sp_dt;
% Q2(6,6) = Sp_dt;
% Q2(7,7) = Sp_dt;
% Q2(8,8) = Sg_dt;

%x_hat_in = [UserPos, UserVel]';
%x_hat = zeros(NumberStates,NumberEpochsINS);
x_hat_out_full = zeros(17,1);

PositionVariance = PositionError ^ 2;
VelocityVariance = VelocityError ^ 2;
ClockBiasVariance = ClockBiasError ^ 2;
ClockDriftVariance = ClockDriftError ^2;
%
P_in = zeros(NumberStates,NumberStates);
% P_in_dash = [PositionVariance,
%         PositionVariance,
%         PositionVariance,
%         ClockBiasVariance,
%         VelocityVariance,
%         VelocityVariance,
%         VelocityVariance,
%         ClockDriftVariance];
% for index = 1:8
%     P_in(index,index) = P_in_dash(index);
% end

P_in(1,1) = PositionVariance/RM;
P_in(2,2) = PositionVariance/RP*cos(-0.47);
P_in(3,3) = PositionVariance;

P_in(4:6,4:6) = eye(3,3)*VelocityVariance;
P_in(16,16) = ClockBiasVariance;
P_in(17,17) = ClockDriftVariance;

P_in_full = P_in;

%% initialise sub filter
for j=1:10
    P_in_sub(:,:,j) = P_in;
end



% iono model parameters for the above Nav file
ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA           
BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA  

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

% two loops will run - a high speed loop calculating teh INS solution at
% 50 Hz - The data is collected at 100Hz, and a Runge-Kutta integration
% scheme is used 

Pfa = 2.22e-8;
SV_TrackTime = zeros(NumberSVs,1);

disp(strcat(ModuleName,sprintf('Beginning navigation loop for %d Epochs',NumberEpochsINS)));

for Epoch_hi = StartTime/ins_dt+1:NumberEpochsINS
    
    g = -gravity(2,Epoch_hi);
    %g = -9.79; 
    
    % implement the inertial navigation solution here
    if Epoch_hi == StartTime/ins_dt+1
%        % first time step
%        a(Epoch_hi) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
%        b(Epoch_hi) = sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) - cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
%        c(Epoch_hi) = cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
%        d(Epoch_hi) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2);
%        
%        % initialise geodetic position

        %dt = 0.02;
        
%         P_xn(Epoch_hi) = InitialPosition(1);
%         P_yn(Epoch_hi) = InitialPosition(2);
%         P_zn(Epoch_hi) = InitialPosition(3);
%         V_xn(Epoch_hi) = InitialVelocity(1);
%         V_yn(Epoch_hi) = InitialVelocity(2);
%         V_zn(Epoch_hi) = InitialVelocity(3);
%         V_n = [V_xn(Epoch_hi);V_yn(Epoch_hi);V_zn(Epoch_hi)];
        
        Pos_LLH(Epoch_hi,:) = InitialPosition;       
        Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch_hi,2),Pos_LLH(Epoch_hi,1));
        %Tecef2ned= T_ECEF2NED(P_yn(Epoch_hi),P_xn(Epoch_hi));
        Tned2ecef = Tecef2ned';
        
        % initialise state vector in ECEF coordinates
        %x_hat(1:3,Epoch_hi) = LLH2ECEF(InitialPosition(1),InitialPosition(2),InitialPosition(3));  % position
        %x_hat(5:7,Epoch_hi) = Tned2ecef * [InitialVelocity(1);InitialVelocity(2);InitialVelocity(3)];  % velocity;

%         % initialise accelerations
%         A_xn(Epoch_hi) = 0;
%         A_yn(Epoch_hi) = 0;
%         A_zn(Epoch_hi) = 0;
        
   else
        ins_dt = TimeINS(Epoch_hi) - TimeINS(Epoch_hi-1);
        
        % get sensor measurements
        omega_x = sensors(5,Epoch_hi);
        omega_y = sensors(6,Epoch_hi);
        omega_z = sensors(7,Epoch_hi);

        A_xb = sensors(2,Epoch_hi);
        A_yb = sensors(3,Epoch_hi);
        A_zb = sensors(4,Epoch_hi);
        
        
%         a_dot(Epoch_hi) = -0.5 * (b(Epoch_hi-1) * omega_x + c(Epoch_hi-1) * omega_y + d(Epoch_hi-1) * omega_z);
%         b_dot(Epoch_hi) = 0.5 * (a(Epoch_hi-1) * omega_x - d(Epoch_hi-1) * omega_y + c(Epoch_hi-1) * omega_z);
%         c_dot(Epoch_hi) = 0.5 * (d(Epoch_hi-1) * omega_x + a(Epoch_hi-1) * omega_y - b(Epoch_hi-1) * omega_z);
%         d_dot(Epoch_hi) = -0.5 * (c(Epoch_hi-1) * omega_x - b(Epoch_hi-1) * omega_y - a(Epoch_hi-1) * omega_z);
%         
%         
%         
%         a(Epoch_hi) = a(Epoch_hi-1) + a_dot(Epoch_hi-1) * ins_dt + (a_dot(Epoch_hi)- a_dot(Epoch_hi-1))*ins_dt/2;
%         b(Epoch_hi) = b(Epoch_hi-1) + b_dot(Epoch_hi-1) * ins_dt + (b_dot(Epoch_hi)- b_dot(Epoch_hi-1))*ins_dt/2;
%         c(Epoch_hi) = c(Epoch_hi-1) + c_dot(Epoch_hi-1) * ins_dt + (c_dot(Epoch_hi)- c_dot(Epoch_hi-1))*ins_dt/2;
%         d(Epoch_hi) = d(Epoch_hi-1) + d_dot(Epoch_hi-1) * ins_dt + (d_dot(Epoch_hi)- d_dot(Epoch_hi-1))*ins_dt/2;
        
        % perform Runge-Kutta itegration for attitude
        
%         a_k1 = a_dot(Epoch_hi-1);
%         a_k2 = a_dot(
%         a(Epoch_hi) = a(Epoch_hi-1) + (ins_dt/6) * (a_k1 + 2*a_k2 + 2 * a_k3 + a_k4); 
        
        
        %Tecef2ned= T_ECEF2NED(P_yn(Epoch_hi-1),P_xn(Epoch_hi-1));
        Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch_hi-1,2),Pos_LLH(Epoch_hi-1,1));
        Tned2ecef = Tecef2ned';
        
        %V_n = [V_xn(Epoch_hi-1);V_yn(Epoch_hi-1);V_zn(Epoch_hi-1)];
    
        


        % find gravity vector for attitude estimate
        gvec_phi(Epoch_hi) =  atan2(-A_yb,sqrt(A_xb^2 + A_zb^2));% roll
        gvec_theta(Epoch_hi) = atan2(A_xb,-A_zb); % pitch

%         c11 = (a(Epoch_hi)^2 + b(Epoch_hi)^2 - c(Epoch_hi)^2 - d(Epoch_hi)^2);
%         c12 = 2 * (b(Epoch_hi)*c(Epoch_hi) - a(Epoch_hi)*d(Epoch_hi));
%         c13 = 2 * (b(Epoch_hi)*d(Epoch_hi) + a(Epoch_hi)*c(Epoch_hi));
%         c21 = 2 * (b(Epoch_hi)*c(Epoch_hi) + a(Epoch_hi)*d(Epoch_hi));
%         c22 = (a(Epoch_hi)^2 - b(Epoch_hi)^2 + c(Epoch_hi)^2 - d(Epoch_hi)^2);
%         c23 = 2 * (c(Epoch_hi)*d(Epoch_hi) - a(Epoch_hi)*b(Epoch_hi));
%         c31 = 2 * (b(Epoch_hi)*d(Epoch_hi) - a(Epoch_hi)*c(Epoch_hi));
%         c32 = 2 * (c(Epoch_hi)*d(Epoch_hi) + a(Epoch_hi)*b(Epoch_hi));
%         c33 = (a(Epoch_hi)^2 - b(Epoch_hi)^2 - c(Epoch_hi)^2 + d(Epoch_hi)^2);

       % C_bn = [c11, c12, c13; c21, c22, c23; c31, c32, c33];
       %
       %C_bn = GARDSim_DCMfromEuler(phi_true(Epoch_hi),theta_true(Epoch_hi),psi_true(Epoch_hi));
       %C_bn = GARDSim_DCMfromQuat([a(Epoch_hi),b(Epoch_hi),c(Epoch_hi),d(Epoch_hi)]);
       
       %     % form the skew symmetric form
        OMEGA_b = [ 0         -omega_z     omega_y; ...
                omega_z    0          -omega_x; ...
               -omega_y    omega_x     0        ];
           
        OMEGA_in = [0, -OMEGAedot * sin(-27*pi/180), 0; ...
                    OMEGAedot * sin(-27*pi/180), 0, -OMEGAedot * cos(-27*pi/180);
                    0, OMEGAedot * cos(-27*pi/180), 0];

        %C_BN = C_BN + (C_BN * OMEGA_b - OMEGA_in * C_BN)*ins_dt;
        C_BN = C_BN + (C_BN * OMEGA_b)*ins_dt;
        
        c32 = C_BN(3,2);
        c33 = C_BN(3,3);
        c31 = C_BN(3,1);
        c21 = C_BN(2,1);
        c11 = C_BN(1,1);
    
        phi_q(Epoch_hi) = atan2(C_BN(3,2),C_BN(3,3));
        theta_q(Epoch_hi) = asin(-C_BN(3,1));
        psi_q(Epoch_hi) = atan2(C_BN(2,1),C_BN(1,1));

        % save the uncorrected attitude
        phi_uc(Epoch_hi) = phi_q(Epoch_hi);
        theta_uc(Epoch_hi) = theta_q(Epoch_hi);
        psi_uc(Epoch_hi) = psi(Epoch_hi);

        % rotate accelerometer measurements to nav frame using A_n = C_bn *
        % A_b

        A_n = C_BN * [A_xb;A_yb;A_zb];
        
%         A_xn(Epoch_hi) = A_n(1);
%         A_yn(Epoch_hi) = A_n(2);
%         A_zn(Epoch_hi) = A_n(3) - LOCAL_GRAVITY;  % gravity is local gravity as computed by aerosim
        Acc_NED(Epoch_hi,:) = A_n - [0; 0; -LOCAL_GRAVITY];
        
        % perform coriolis correction
        %OMEGA_e_n = Tecef2ned * [0;0;7.292115e-5]; % rotation of the earth
        %OMEGA_n = Tecef2ned * [lat_dot(Epoch_hi-1);long_dot(Epoch_hi-1);0]; % rotation of the reference frame    
        %Coriolis(:,Epoch_hi) = cross(2 * OMEGA_e_n + OMEGA_n,V_n);
        
        % calculate coriolis correction as per Aerosim
        OMEGA_e = 7.292115e-5; % eearth rotation rate (inertial frame)
        
        
        CorMat(1,1) = 0;
        CorMat(1,2) = (2*OMEGA_e + long_dot(Epoch_hi-1))*sin(Pos_LLH(Epoch_hi-1,1));
        CorMat(1,3) = -lat_dot(Epoch_hi-1);
        CorMat(2,1) = -(2*OMEGA_e + long_dot(Epoch_hi-1))*sin(Pos_LLH(Epoch_hi-1,1));
        CorMat(2,2) = 0;
        CorMat(2,3) = -(2*OMEGA_e + long_dot(Epoch_hi-1))*cos(Pos_LLH(Epoch_hi-1,1));
        CorMat(3,1) = lat_dot(Epoch_hi-1);
        CorMat(3,2) = (2*OMEGA_e + long_dot(Epoch_hi-1))*cos(Pos_LLH(Epoch_hi-1,1));
        CorMat(3,3) = 0;
        
        Coriolis(Epoch_hi,:) = CorMat * Vel_NED(Epoch_hi-1,:)';
        
%         A_xn(Epoch_hi) = A_xn(Epoch_hi) + Coriolis(1);
%         A_yn(Epoch_hi) = A_yn(Epoch_hi) + Coriolis(2);
%         A_zn(Epoch_hi) = A_zn(Epoch_hi) + Coriolis(3);

        Acc_NED(Epoch_hi,:) = Acc_NED(Epoch_hi,:) + Coriolis(Epoch_hi,:);
        
%         % rotate acceleration measurement to ECEF FRAME
%         A_n = [A_xn(Epoch_hi);A_yn(Epoch_hi);A_zn(Epoch_hi)];
%         A_e = Tned2ecef * A_n;
        
        % TODO perform Runge-Kutta 4th order integration here
        
        % Trapezoidal INTEGRATION POSITION UPDATE       
        % propogate velocity estimate
        %V_xn(Epoch_hi) = V_xn(Epoch_hi-1) + A_xn(Epoch_hi-1) * ins_dt + (A_xn(Epoch_hi) - A_xn(Epoch_hi-1))*ins_dt/2;
        %V_yn(Epoch_hi) = V_yn(Epoch_hi-1) + A_yn(Epoch_hi-1) * ins_dt + (A_yn(Epoch_hi) - A_yn(Epoch_hi-1))*ins_dt/2;
        %V_zn(Epoch_hi) = V_zn(Epoch_hi-1) + A_zn(Epoch_hi-1) * ins_dt + (A_zn(Epoch_hi) - A_zn(Epoch_hi-1))*ins_dt/2;

        % determine latitude and longitude rates

%         lat_dot(Epoch_hi) = V_xn(Epoch_hi-1) / (RM + P_zn(Epoch_hi-1));
%         long_dot(Epoch_hi) = V_yn(Epoch_hi-1) / (cos(P_xn(Epoch_hi-1))*(RP + P_zn(Epoch_hi-1)));

        % propogate position estimate
        %P_xn(Epoch_hi) = P_xn(Epoch_hi-1) + lat_dot(Epoch_hi-1) * ins_dt + (lat_dot(Epoch_hi) - lat_dot(Epoch_hi-1))*ins_dt/2;
        %P_yn(Epoch_hi) = P_yn(Epoch_hi-1) + long_dot(Epoch_hi-1) * ins_dt + (long_dot(Epoch_hi) - long_dot(Epoch_hi-1))*ins_dt/2;
        %P_zn(Epoch_hi) = P_zn(Epoch_hi-1) + V_zn(Epoch_hi-1) * ins_dt  + (V_zn(Epoch_hi) - V_zn(Epoch_hi-1))*ins_dt/2;

        
        
%         
%         x_hat(4:6,Epoch_hi) = x_hat(4:6,Epoch_hi-1) + A_e * ins_dt;
%         x_hat(17,Epoch_hi) = x_hat(17,Epoch_hi-1);
%         
%         x_hat(1:3,Epoch_hi) = x_hat(1:3,Epoch_hi-1) + x_hat(4:6,Epoch_hi-1) * ins_dt + (x_hat(4:6,Epoch_hi) - x_hat(4:6,Epoch_hi-1))*ins_dt/2;
%         x_hat(16,Epoch_hi) = x_hat(16,Epoch_hi-1) + x_hat(17,Epoch_hi-1) * ins_dt;
%         
%         [P_xn(Epoch_hi),P_yn(Epoch_hi),P_zn(Epoch_hi)] = ECEF2LLH(x_hat(1:3,Epoch_hi));
%         V_n = Tecef2ned * x_hat(4:6,Epoch_hi);
%         V_xn(Epoch_hi) = V_n(1); V_yn(Epoch_hi) = V_n(2);V_zn(Epoch_hi) = V_n(3);
%         


        % increment the velocity and position 

        if(Epoch_hi == 1)
            Vel_NED(Epoch_hi,:) = InitialVelocity; %Acc_NED(Epoch_hi,:)*dt*0;
            Pos_LLH(Epoch_hi,:) = InitialPosition;
        else
            % trapezoidal velocity update
            Vel_NED(Epoch_hi,:) = Vel_NED(Epoch_hi-1,:) + Acc_NED(Epoch_hi,:)*ins_dt;% + (Acc_NED(Epoch_hi,:)-Acc_NED(Epoch_hi-1,:))/2;
            
            % position update
            lat_dot(Epoch_hi) = Vel_NED(Epoch_hi,1) / (RM + Pos_LLH(Epoch_hi-1,3));
            long_dot(Epoch_hi) = Vel_NED(Epoch_hi,2) / (cos(Pos_LLH(Epoch_hi-1,1))*(RP + Pos_LLH(Epoch_hi-1,3)));
            h_dot(Epoch_hi) =  -Vel_NED(Epoch_hi,3);
            
            Pos_LLH(Epoch_hi,1) = Pos_LLH(Epoch_hi-1,1) + lat_dot(Epoch_hi)*ins_dt;
            Pos_LLH(Epoch_hi,2) = Pos_LLH(Epoch_hi-1,2) + long_dot(Epoch_hi)*ins_dt;
            Pos_LLH(Epoch_hi,3) = Pos_LLH(Epoch_hi-1,3) + h_dot(Epoch_hi)*ins_dt;

            UserClock(Epoch_hi,2) = UserClock(Epoch_hi-1,2);
            UserClock(Epoch_hi,1) = UserClock(Epoch_hi-1,1) + UserClock(Epoch_hi,2)*ins_dt;
            
        end
    
        % propagate the state error
        x_hat_minus = phi2 * x_hat_minus;
        % propagate the covariance matrix
        P_in_full = phi2 * P_in_full * phi2' + Q2d;
        
        P_out_diag(:,Epoch_hi) = diag(P_in_full);
        %% calculate error bounds
        Pos_var(Epoch_hi,:) = Tecef2ned * P_out_diag(1:3,Epoch_hi);
        Pos_var(Epoch_hi,1) = sqrt(abs(Pos_var(Epoch_hi,1)));
        Pos_var(Epoch_hi,2) = sqrt(abs(Pos_var(Epoch_hi,2)));
        Pos_var(Epoch_hi,3) = sqrt(abs(Pos_var(Epoch_hi,3)));
        
        % low-speed loop - kalman filter
        if(mod(Epoch_hi-1,100) == 0 && Epoch_hi ~= 1)
        %if(0)  %% uncomment this line to turn off the KF update
            Epoch_lo = Epoch_lo + 1;
            
            if(Epoch_lo == StartTime + 60)  % enable fault detection
                FDEnabled = 1;
            end
            
            
            clear PR_Vec PRR_Vec H R K V z_Vec NumberGPSMeasurements PR_Vec_minus PRR_Vec_minus;
            
%             if(Epoch_lo == NoiseOnTime)
%                 RangeNoiseVariance = RangeNoiseVariance + NoiseAmplification^2;
%             end

            % get the gps time - used for the iono correction
            GPSTime  = 259199 + Epoch_lo;

%             % UserPos and UserVel are the current state vector
%             % which the EKF will linearise on
%             UserPos(1:3) = x_hat(1:3,Epoch_hi);
%             UserPos(4) = x_hat(16,Epoch_hi);
%             UserVel(1:3) = x_hat(4:6,Epoch_hi);
%             UserVel(4) = x_hat(17,Epoch_hi);
%             

            UserPos(1:3) =  LLH2ECEF(Pos_LLH(Epoch_hi,1),Pos_LLH(Epoch_hi,2),Pos_LLH(Epoch_hi,3));
            UserPos(4) = UserClock(Epoch_hi,1);
            
            UserVel(1:3) = Tned2ecef * Vel_NED(Epoch_hi,:)';
            UserVel(4) = UserClock(Epoch_hi,2);
            
            
            %% introduce a ramp fault onto PRN 1
%             if(Epoch_lo > 500)
%                 PR_Error = PR_Error + PR_Error_Rate * gps_dt;
%                 %PR_Error = 40;
%                 PR_Sim(Epoch_lo,3) = PR_Sim(Epoch_lo,3) + PR_Error;
%             end
            
            % format the input vectors
            SVIndex = 0;
            for SV=1:NumberSVs
                if((PR_Sim(Epoch_lo,SV) > 100) && SVDontUse(SV) == 0)  %%SV ~= ExcludeSVPRN
                    % add to PR vector
                    SVIndex = SVIndex + 1;
                    
                   SV_TrackTime(SV) = SV_TrackTime(SV) + 1;
                    
                    SV_Vec(SVIndex) = SV;
                    PR_Vec(SVIndex) = PR_Sim(Epoch_lo,SV);
                    
                    
                    
                    %PRR_Vec(SVIndex) = PRR_Sim(Epoch_lo,SV);
                    
                    if(Epoch_lo == StartTime+1)
                        PRR_Vec(SVIndex) = 0;
                        PR_csc(SVIndex) = PR_Vec(SVIndex);
                    else
                        if(CP_Sim(Epoch_lo-1,SV) == 0)
                            PRR_Vec(SVIndex) = (CP_Sim(Epoch_lo+1,SV) - CP_Sim(Epoch_lo,SV));
                            PR_csc(SVIndex) = PR_Vec(SVIndex);
                        else
                            PRR_Vec(SVIndex) = (CP_Sim(Epoch_lo,SV) - CP_Sim(Epoch_lo-1,SV));
                            
                         %%% PERFORM CARRIER PHASE SMOOTHING of PR
                         %%% HERE
                            
                            
                            alpha = gps_dt/100;
                            P_proj = PR_prev + (L1_Wavelength/(2*pi))*(CP_Sim(Epoch_lo) - CP_Sim(Epoch_lo-1));
                            PR_csc(SVIndex) = alpha * PR_Vec(SVIndex) + (1-alpha)*P_proj;                            
                            
                        end
                    end
%                    save previous value
                   PR_prev = PR_csc(SVIndex);
                   
                    % dont use the carrier smoothed code
                    %PR_Vec = PR_csc;                  
                    
                    
                    
                    SVPos(SVIndex,:) = [SV_X_Data(Epoch_lo,SV) SV_Y_Data(Epoch_lo,SV) SV_Z_Data(Epoch_lo,SV) SV_T_Data(Epoch_lo,SV)];
                    SVVel(SVIndex,:) = [SV_Xvel_Data(Epoch_lo,SV) SV_Yvel_Data(Epoch_lo,SV) SV_Zvel_Data(Epoch_lo,SV) SV_Tvel_Data(Epoch_lo,SV)];
                    SVAcc(SVIndex,:) = [SV_Xacc_Data(Epoch_lo,SV) SV_Yacc_Data(Epoch_lo,SV) SV_Zacc_Data(Epoch_lo,SV) SV_Tacc_Data(Epoch_lo,SV)];

                   
                    % calculate hte earth rotation correction as per Kayton pg 228
                    % eq 5.67

                    delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) *UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
                    PR_Vec(SVIndex) = PR_Vec(SVIndex) + delta_pr_omegaedot;
                    
                    % correct for satellite clock and drift
                    PR_Vec(SVIndex) = PR_Vec(SVIndex)  + SVPos(SVIndex,4);
                    PRR_Vec(SVIndex) = PRR_Vec(SVIndex) + SVVel(SVIndex,4);
            
                    % calculate the iono delay correction - single frequency user
                    % model from ICD 200
                    if(GRASOn == 0)
                        ionodelay = ionomodel(GPSTime, UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
                        PR_Vec(SVIndex) = PR_Vec(SVIndex) - ionodelay;
                    end
                   
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

            SatellitesUsed(Epoch_lo) = NumberGPSMeasurements;
            
            %x_hat_in = x_hat(:,Epoch_hi);
            %x_hat_in = x_hat_out_full;
                        
                %%% FIRST DO THE FULL SOLUTION

                % R matrix
                R = eye(2*NumberGPSMeasurements);
                for k = 1:NumberGPSMeasurements
                    R(k,k) = RangeNoiseVariance;
                    R(k+NumberGPSMeasurements,k+NumberGPSMeasurements) = RangeRateNoiseVariance;
                end

                % determine the H matrix
                H = zeros(2*NumberGPSMeasurements,NumberStates);
                for k = 1:NumberGPSMeasurements
                    %Calculated slant ranges

                    for m = 1:3
                         ele(m) =  SVPos(k,m) - UserPos(m);
                    end    

                    r_VecCalc(k) =  norm(ele);   

%                     H(k,1) =  -ele(1)/r_VecCalc(k);
%                     H(k,2) =  -ele(2)/r_VecCalc(k);
%                     H(k,3) =  -ele(3)/r_VecCalc(k);
%                     H(k,16) = 1.0;   
% %                     H(k,5) = 0;
% %                     H(k,6) = 0;
% %                     H(k,7) = 0;
% %                     H(k,8) = 0;
% 
% %                     H(k+NumberGPSMeasurements,1) = 0;
% %                     H(k+NumberGPSMeasurements,2) = 0;
% %                     H(k+NumberGPSMeasurements,3) = 0;
% %                     H(k+NumberGPSMeasurements,4) = 0;   
%                     H(k+NumberGPSMeasurements,4) = -ele(1)/r_VecCalc(k);
%                     H(k+NumberGPSMeasurements,5) = -ele(2)/r_VecCalc(k);
%                     H(k+NumberGPSMeasurements,6) = -ele(3)/r_VecCalc(k);
%                     H(k+NumberGPSMeasurements,17) = 1.0;

                    
                    
%%%%%%%% TROYS PARTIAL DERIVATIVES  %%%%%%%%%
                    

        %H matrix elements for the range rates:

        LatitudeGPS=Pos_LLH(Epoch_hi,1);
        LongitudeGPS=Pos_LLH(Epoch_hi,2);
        HeightGPS=Pos_LLH(Epoch_hi,3);

        dxdot_dVn = -cos(LongitudeGPS)*sin(LatitudeGPS);
        dxdot_dVe = -sin(LongitudeGPS);
        dxdot_dVd = -cos(LongitudeGPS)*cos(LatitudeGPS);

        dydot_dVn = -sin(LongitudeGPS)*sin(LatitudeGPS);
        dydot_dVe = cos(LongitudeGPS);
        dydot_dVd = -sin(LongitudeGPS)*cos(LatitudeGPS);

        dzdot_dVn = cos(LatitudeGPS);
        dzdot_dVd = -sin(LatitudeGPS);    

        drhodot_dxdot = -ele(1)/r_VecCalc(k);           
        drhodot_dydot = -ele(2)/r_VecCalc(k);                       
        drhodot_dzdot = -ele(3)/r_VecCalc(k);
           

         %partial derivatives

        drhodot_dVn = drhodot_dxdot*dxdot_dVn + drhodot_dydot*dydot_dVn + drhodot_dzdot*dzdot_dVn;
        drhodot_dVe = drhodot_dxdot*dxdot_dVe + drhodot_dydot*dydot_dVe;
        drhodot_dVd = drhodot_dxdot*dxdot_dVd + drhodot_dydot*dydot_dVd + drhodot_dzdot*dzdot_dVd;

        drho_dx = -ele(1)/r_VecCalc(k); 
        drho_dy = -ele(2)/r_VecCalc(k); 
        drho_dz = -ele(3)/r_VecCalc(k); 

        %H matrix elements for the ranges:

        %note that dx_dPn is the same as dxdot_dVn

        %      %this is from the chain rule of course.
        %      drho_dPn = drho_dx*dxdot_dVn + drho_dy*dydot_dVn + drho_dz*dzdot_dVn;
        %      drho_dPe = drho_dx*dxdot_dVe + drho_dy*dydot_dVe;
        %      drho_dPd = drho_dx*dxdot_dVd + drho_dy*dydot_dVd + drho_dz*dzdot_dVd;
        %       


        %this is the lat lon hgt error version
        %this is from the chain rule of course.



        dxdot_dVn = sin(LongitudeGPS)*cos(LatitudeGPS);
        dxdot_dVe = -cos(LongitudeGPS);
        dxdot_dVd = -sin(LongitudeGPS)*sin(LatitudeGPS);

        dydot_dVn = -cos(LongitudeGPS)*cos(LatitudeGPS);
        dydot_dVe = -sin(LongitudeGPS);
        dydot_dVd = cos(LongitudeGPS)*sin(LatitudeGPS);

        dzdot_dVn = -sin(LatitudeGPS);
        dzdot_dVd = -cos(LatitudeGPS);    



        lon = LongitudeGPS;
        lat = LatitudeGPS;
        hgt = HeightGPS;

        %         dxdlat = -cos(lon)*cos(lat)*lat*(Re+hgt)-cos(lon)*sin(lat)*(Re+hgt)+sin(lon)*lon*(Rn+hgt)*sin(lat)-cos(lon)*sin(lat)*hgt;
        %         dxdlon = sin(lon)*sin(lat)*lat*(Re+hgt)-cos(lon)*lon*(Rn+hgt)*cos(lat)-sin(lon)*(Rn+hgt)*cos(lat)-sin(lon)*cos(lat)*hgt;
        %         dxdhgt = -cos(lon)*sin(lat)*lat-sin(lon)*lon*cos(lat)+cos(lon)*cos(lat);
        %         
        %         dydlat = -sin(lon)*cos(lat)*lat*(Re+hgt)-sin(lon)*sin(lat)*(Re+hgt)-cos(lon)*lon*(Rn+hgt)*sin(lat)-sin(lon)*sin(lat)*hgt;
        %         dydlon = -cos(lon)*sin(lat)*lat*(Re+hgt)-sin(lon)*lon*(Rn+hgt)*cos(lat)+cos(lon)*(Rn+hgt)*cos(lat)+cos(lon)*cos(lat)*hgt;  
        %         dydhgt = -sin(lon)*sin(lat)*lat+cos(lon)*lon*cos(lat)+sin(lon)*cos(lat);


        %     dzdlat = -sin(lat)*lat*(Re+hgt)+cos(lat)*(Re+hgt)+cos(lat)*hgt;
        %     dzdlon = 0;
        %     dzdhgt = cos(lat)*lat+sin(lat);


        NN = a/(sqrt(1-e^2*sin(lat)^2));

        dxdlat = -(NN + hgt)*sin(lat)*cos(lon);
        dxdlon = -(NN + hgt)*cos(lat)*sin(lon);
        dxdhgt = cos(lat)*cos(lon);

        dydlat = -(NN + hgt)*sin(lat)*sin(lon);
        dydlon = (NN + hgt)*cos(lat)*cos(lon);
        dydhgt = cos(lat)*sin(lon);


        dzdlat = (NN+hgt)*cos(lat);
        dzdlon = 0;
        dzdhgt = sin(lat);




        drho_dPn = drho_dx*dxdlat + drho_dy*dydlat + drho_dz*dzdlat;
        drho_dPe = drho_dx*dxdlon + drho_dy*dydlon;
        drho_dPd = drho_dx*dxdhgt + drho_dy*dxdhgt + drho_dz*dzdhgt;






        %       %FOR THE INS
        %this is for the ranges
        H(k,1) =  drho_dx;
        H(k,2) =  drho_dy;
        H(k,3) =  drho_dz;
        H(k,4) = 0;
        H(k,5) = 0;
        H(k,6) = 0;
        H(k,7) = 0;  
        H(k,8) = 0;
        H(k,9) = 0;
        H(k,16) = 1;  %clock bias



        %this is for the range rates

        H(k+NumberGPSMeasurements,1) = 0;
        H(k+NumberGPSMeasurements,2) = 0;
        H(k+NumberGPSMeasurements,3) = 0;
        H(k+NumberGPSMeasurements,4) = drhodot_dxdot;
        H(k+NumberGPSMeasurements,5) = drhodot_dydot;
        H(k+NumberGPSMeasurements,6) = drhodot_dzdot;
        H(k+NumberGPSMeasurements,7) = 0;
        H(k+NumberGPSMeasurements,8) = 0;
        H(k+NumberGPSMeasurements,9) = 0;     
        H(k+NumberGPSMeasurements,17) = 1; %clock drift        

                    
%%%%%%%%% END OF TROYS PARTIAL DERIVATIVES %%%%%%%%%%%
                    
                    

                    
                    % find apriori estimate of pseudorange
                    %PR_Vec_minus(k) = r_VecCalc(k) + UserPos(4) + UserVel(4) * gps_dt;  % geometric range + c * delta_T
                    PR_Vec_minus(k) = r_VecCalc(k) + UserPos(4);  % geometric range + c * delta_T
                    
                    %predicted relative velocity of sv and receiver
                    r_VecCalcVel(k) = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + ...
                                      (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + ...
                                      (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
                    Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);

                    PRR_Vec_minus(k) = Relative_Velocity(k) + UserVel(4);
                end
                
                H(1:NumberGPSMeasurements,1:3) = H(1:NumberGPSMeasurements,1:3)*Tecef2ned';
                H(NumberGPSMeasurements+1:2*NumberGPSMeasurements,4:6) = H(NumberGPSMeasurements+1:2*NumberGPSMeasurements,4:6)*Tecef2ned';
                


                % calculate DOPS 
                H4(1:NumberGPSMeasurements,1:3) = H(1:NumberGPSMeasurements,1:3);
                H4(1:NumberGPSMeasurements,4) = H(1:NumberGPSMeasurements,16);
                H3 = H(1:NumberGPSMeasurements,1:3);
                Hh = H3;
                 
                H(:,1) = H(:,1)*RM;    
                H(:,2) = H(:,2)*RP*cos(LatitudeGPS); %I think this is right, Reh*cos(lat), not Rnh*cos(lat)    
                H(:,3) = H(:,3)*-1;
    
            
                %GDOP(Epoch_lo) = sqrt(trace(inv(H4' * H4)));
                %PDOP(Epoch_lo) = sqrt(trace(inv(H3' * H3)));
                
                AA = inv(H4' * H4);
                %HDOP(Epoch_lo) = sqrt(trace(Hh_H(1:2,1:2)));
                %VDOP(Epoch_lo) = sqrt(Hh_H(3,3));

                var_x = AA(1,1);
                var_y = AA(2,2);
                var_z = AA(3,3);
                var_dt = AA(4,4);

                GDOP(Epoch_lo) = sqrt(var_x + var_y + var_z + var_dt);

                PDOP(Epoch_lo) = sqrt(var_x + var_y + var_z);

                HDOP(Epoch_lo) = sqrt(var_x + var_y);

                VDOP(Epoch_lo) = sqrt(var_z);

                TDOP(Epoch_lo) = sqrt(var_dt);
                
                % calculate the elevations of each SV
                for k = 1:NumberGPSMeasurements
                    Elevation(Epoch_lo,SV_Vec(k)) = asin(Hh(k,3));
           
                end
                

                % generate PHI matrix
                F = GARD_GenerateINSFMatrix(Pos_LLH(Epoch_hi,:),Vel_NED(Epoch_hi,:),Acc_NED(Epoch_hi,:)-[0 0 LOCAL_GRAVITY]);
                phi = eye(NumberStates,NumberStates);
                phi(1:9,1:9) = phi(1:9,1:9) + F * gps_dt;
                % % augmented state matrix
                phi(7:9,10:12) = -C_BN * gps_dt;
                phi(4:6,13:15) = C_BN * gps_dt;

                phi(10,10) = 1 - x_accel_beta * gps_dt;
                phi(11,11) = 1 - y_accel_beta * gps_dt;
                phi(12,12) = 1 - z_accel_beta * gps_dt;
                phi(13,13) = 1 - x_gyro_beta * gps_dt;
                phi(14,14) = 1 - y_gyro_beta * gps_dt;
                phi(15,15) = 1 - z_gyro_beta * gps_dt;
                phi(16,17) = 1.0*gps_dt;

                phi2 = eye(NumberStates,NumberStates);
                phi2(1:9,1:9) = phi2(1:9,1:9) + F * ins_dt;
                % % augmented state matrix
                phi2(7:9,10:12) = -C_BN * ins_dt;
                phi2(4:6,13:15) = C_BN * ins_dt;

                phi2(10,10) = 1 - x_accel_beta * ins_dt;
                phi2(11,11) = 1 - y_accel_beta * ins_dt;
                phi2(12,12) = 1 - z_accel_beta * ins_dt;
                phi2(13,13) = 1 - x_gyro_beta * ins_dt;
                phi2(14,14) = 1 - y_gyro_beta * ins_dt;
                phi2(15,15) = 1 - z_gyro_beta * ins_dt;
                phi2(16,17) = 1.0*ins_dt;



                Sf = 8*pi^2*2e-20*Speedoflight^2;  % frequency noise psd
                Sc = 2*2e-19*Speedoflight^2;  % bias noise psd

                Q = zeros(NumberStates,NumberStates);

                Q(10,10) = 2*x_gyro_beta*x_gyro_Q;
                Q(11,11) = 2*y_gyro_beta*y_gyro_Q;
                Q(12,12) = 2*z_gyro_beta*z_gyro_Q;
                Q(13,13) = 2*x_accel_beta*x_accel_Q;
                Q(14,14) = 2*y_accel_beta*y_accel_Q;
                Q(15,15) = 2*z_accel_beta*z_accel_Q;
                
                Q(16,16) = 10*Sc * gps_dt + Sf * (gps_dt^2)/2;
                Q(16,17) = 10*Sc * (gps_dt^2) / 2;
                Q(17,16) = 10*Sc * (gps_dt^2) / 2;
                Q(17,17) = 10*Sf * gps_dt;

                G = zeros(NumberStates,NumberStates);
                G(7:9,10:12) = -C_BN;
                G(4:6,13:15) = C_BN;
                G(10:12,10:12) = eye(3,3);
                G(13:15,13:15) = eye(3,3);

                % convert to discrete time form
                Qd = phi * (G * Q * G') * phi';
                Q2d = phi2 * (G * Q * G') * phi2';

                % find measurement vector - delta between predicted and measured
                % pseudorange
                z_Vec = ([PR_Vec PRR_Vec]' - [PR_Vec_minus PRR_Vec_minus]');

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
                x_hat_out_full = x_hat_minus + K * (z_Vec);
                %x_hat_out_full = x_hat_in +  K * (z_Vec);

                %P_out_full = P_minus - K * H * P_minus;
                P_out_full = (eye(NumberStates,NumberStates) - K*H)*P_minus*(eye(NumberStates,NumberStates)-K*H)' + K*R*K';
    
                P_in_full = P_out_full;
                
                P_out_diag(:,Epoch_hi) = diag(P_out_full);
             
                
                z_save(1:2*NumberGPSMeasurements,Epoch_lo) = z_Vec;
                x_save(:,Epoch_lo) = x_hat_out_full;

                %% apply correction
                Vel_NED(Epoch_hi,:) = Vel_NED(Epoch_hi,:) + x_hat_out_full(4:6)';
                Pos_LLH(Epoch_hi,:) = Pos_LLH(Epoch_hi,:) + x_hat_out_full(1:3)';
                UserClock(Epoch_hi,1) = UserClock(Epoch_hi,1) + x_hat_out_full(16);
                UserClock(Epoch_hi,2) = UserClock(Epoch_hi,2) + x_hat_out_full(17);

                % convert tilt error to DCM Update
                del_alpha = x_hat_out_full(7);
                del_beta = x_hat_out_full(8);
                del_gamma = x_hat_out_full(9);

                % save bias correction
                GyroBias(Epoch_lo,:) = x_hat_out_full(13:15)';
                AccelBias(Epoch_lo,:) = x_hat_out_full(10:12)';
                
                x_hat_out_full(1:9) = 0;
                x_hat_out_full(16:17) = 0;
                
                x_hat_minus = x_hat_out_full;

                % reset state vector
                %x_hat(1:9) = 0;


%                 % correct atttitude
%                 del_att_skew = [0         -del_gamma   del_beta; ...
%                                del_gamma  0          -del_alpha; ...
%                                -del_beta  del_alpha   0];
% 
%                 C_BN = (eye(3,3) + del_att_skew) * C_BN;
% 
%                 C_BN = C_BN / norm(C_BN);
                
                
                %%% END FULL SOLUTION
%                 
%             if(NumberGPSMeasurements < 4)
%                 % do something
%             elseif (NumberGPSMeasurements == 4)
%                 % position solution - no integrity
%             elseif (NumberGPSMeasurements > 4)
%                 % special case of solution separation
%             %elseif (NumberGPSMeasurements > 5)
%                 % full solution separation
%                 
%                 
%                 
%                 
%                 
%                 %%% START SUB SOLUTIONS
%                 % do a position solution for each measurement
%                 NumberSubMeasurements = NumberGPSMeasurements-1;
%                 
% 
%                 for (SubSolution = 1:NumberGPSMeasurements)
% 
%                     
%                     % form the PR and CP vectors with one measurement left
%                     % out
%                     
%                     PR_SubVec(1:SubSolution-1) = PR_Vec(1:SubSolution-1);
%                     PR_SubVec(SubSolution:NumberSubMeasurements) = PR_Vec(SubSolution+1:NumberGPSMeasurements);
%                     
%                     PRR_SubVec(1:SubSolution-1) = PRR_Vec(1:SubSolution-1);
%                     PRR_SubVec(SubSolution:NumberSubMeasurements) = PRR_Vec(SubSolution+1:NumberGPSMeasurements);
%                     
%                     SVPos_Sub(1:SubSolution-1,:) = SVPos(1:SubSolution-1,:);
%                     SVPos_Sub(SubSolution:NumberSubMeasurements,:) = SVPos(SubSolution+1:NumberGPSMeasurements,:);
%                     
%                     SVVel_Sub(1:SubSolution-1,:) = SVVel(1:SubSolution-1,:);
%                     SVVel_Sub(SubSolution:NumberSubMeasurements,:) = SVVel(SubSolution+1:NumberGPSMeasurements,:);
%                     
%                     SVAcc_Sub(1:SubSolution-1,:) = SVAcc(1:SubSolution-1,:);
%                     SVAcc_Sub(SubSolution:NumberSubMeasurements,:) = SVAcc(SubSolution+1:NumberGPSMeasurements,:);
%                     
%                     % R matrix
%                     R = eye(2*NumberSubMeasurements);
%                     for k = 1:NumberSubMeasurements
%                         R(k,k) = RangeNoiseVariance;
%                         R(k+NumberSubMeasurements,k+NumberSubMeasurements) = RangeRateNoiseVariance;
%                     end
% 
%                     % determine the H matrix
%                     H = zeros(NumberSubMeasurements,NumberStates);
%                     for k = 1:NumberSubMeasurements
%                         %Calculated slant ranges
% 
%                         for m = 1:3
%                              ele(m) =  SVPos_Sub(k,m) - UserPos(m);
%                         end    
% 
%                         r_VecCalc(k) =  norm(ele);   
% 
%                         H(k,1) =  -ele(1)/r_VecCalc(k);
%                         H(k,2) =  -ele(2)/r_VecCalc(k);
%                         H(k,3) =  -ele(3)/r_VecCalc(k);
%                         H(k,4) = 1;   
%                         H(k,5) = 0;
%                         H(k,6) = 0;
%                         H(k,7) = 0;
%                         H(k,8) = 0;
% 
%                         H(k+NumberSubMeasurements,1) = 0;
%                         H(k+NumberSubMeasurements,2) = 0;
%                         H(k+NumberSubMeasurements,3) = 0;
%                         H(k+NumberSubMeasurements,4) = 0;   
%                         H(k+NumberSubMeasurements,5) = -ele(1)/r_VecCalc(k);
%                         H(k+NumberSubMeasurements,6) = -ele(2)/r_VecCalc(k);
%                         H(k+NumberSubMeasurements,7) = -ele(3)/r_VecCalc(k);
%                         H(k+NumberSubMeasurements,8) = 1;
% 
%                         % find apriori estimate of pseudorange
%                         PR_SubVec_minus(k) = r_VecCalc(k) + UserPos(4) + UserVel(4) * gps_dt;  % geometric range + c * delta_T
% 
%                         %predicted relative velocity of sv and receiver
%                         r_VecCalcVel(k) = (SVVel_Sub(k,1) - UserVel(1))*(SVPos_Sub(k,1)-UserPos(1)) + ...
%                                           (SVVel_Sub(k,2) - UserVel(2))*(SVPos_Sub(k,2)-UserPos(2)) + ...
%                                           (SVVel_Sub(k,3) - UserVel(3))*(SVPos_Sub(k,3)-UserPos(3));
%                                       
%                         Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);
% 
%                         PRR_SubVec_minus(k) = Relative_Velocity(k) + UserVel(4);
%                     end
% 
%                     % find measurement vector - delta between predicted and measured
%                     % pseudorange
%                     z_Vec = ([PR_SubVec PRR_SubVec]' - [PR_SubVec_minus PRR_SubVec_minus]');
% 
%                     % evaluate the KF
%                     %[x_hat_out, P_out] = GARD_EvaluateKF(dt, x_hat_in, P_in, phi, H, z_Vec, Q, R);
%                     % Prediction step
%                     x_hat_minus = x_hat_in;
%                     P_minus = phi * P_in_sub(:,:,SubSolution) * phi' + Q;
%                     
%                     
%                     % correction step
% 
%                     % calculate the kalman gain
%                     V = H * P_minus * H' + R;
%                     K = P_minus * H' * inv(V);
% 
%                     % update equations
%                     %v = z - H * x_hat_minus;
%                     x_hat_out(:,SubSolution) = x_hat_minus + K * (z_Vec);
%                     P_out(:,:,SubSolution) = P_minus - K * H * P_minus;
% 
%                     P_in_sub(:,:,SubSolution) = P_out(:,:,SubSolution);
%                     
%                     % clear variables
%                     clear PR_SubVec PRR_SubVec SVPos_Sub SVVel_Sub SVAcc_Sub PR_SubVec_minus PRR_SubVec_minus R H;  
% 
%                 end % end sub models
                
                %%%% END SUB SOLUTIONS
                
                %%%% START FAULT DETECTION
                
%                 
%                 
%                 % form the solution separation vectors
%                 for j = 1:NumberGPSMeasurements
%                     % note that only the position estimates are used - not
%                     % clock or velocity
%                     Beta_ss(:,j) = x_hat_out_full(1:3) - x_hat_out(1:3,j);
%                     B_ss(:,:,j) = P_out(1:3,1:3,j) - P_out_full(1:3,1:3);
%                     
%                     lamda_ss(j) = Beta_ss(:,j)' * pinv(B_ss(:,:,j)) * Beta_ss(:,j);
%                     
%                     B_lambda = eigs(B_ss(:,:,j));
%                     TD(Epoch_lo,j) = sqrt(max(B_lambda)) * abs(norminv(Pfa/NumberGPSMeasurements,0,1));
%                     
%                     PPL_sub(Epoch_lo,j) = sqrt(P_out(1,1,j)^2 + P_out(2,2,j)^2 + P_out(3,3,j)^2) + TD(Epoch_lo,j);
%                     
%                     
%                     P_ned = Tecef2ned * P_out(1:3,1:3,j);
%                     
%                     HPL_sub(Epoch_lo,j) = sqrt(P_ned(1,1)^2 + P_ned(2,2)^2) + TD(Epoch_lo,j);
%                     VPL_sub(Epoch_lo,j) = sqrt(P_ned(3,3)^2) + TD(Epoch_lo,j);
%                     
%                     
%                 end
%                 PPL(Epoch_lo) = max(PPL_sub(Epoch_lo,:));
%                 HPL(Epoch_lo) = max(HPL_sub(Epoch_lo,:));
%                 VPL(Epoch_lo) = max(VPL_sub(Epoch_lo,:));
%                 
%                 Beta_ss_out = Beta_ss;
%                 B_ss_out = B_ss;
%                 lamda_ss_out(Epoch_lo,1:length(lamda_ss)) = lamda_ss;
%                 
%                 %% we are really only interested in the position solution
%                 %% here - should hte clock estimate be included in the test
%                 %% statistic??
%                 
%                 % rotate solution separation vector from ECEF to navigation
%                 % frame (NED)
%                 
%                 % rotate solution separation covariance matrix 
%                 
%                 %% for each satellite used, if hte minimum track time is
%                 %% less than 10 seconds, do not do fault detection
%                 for i=1:NumberGPSMeasurements
%                    TrackTime(i) = SV_TrackTime(SV_Vec(i)); 
%                 end
% 
%                 MinTrackTime = min(TrackTime);
%                 
%                 if MinTrackTime == 1
%                     disp(sprintf('New Satellite Set at time %d',Epoch_lo));
%                 end
%                 
%                 %% test the lamda vector for faults
%                 faults = find(lamda_ss > TD(Epoch_lo));
%                 
%                 NumFaults = length(faults);
%                 
%                 if(FDEnabled && MinTrackTime > 10)
%                     if(NumFaults == 0)
% 
%                         % choose which solution to use
%                         UserPos = x_hat_out_full(1:4)';
%                         UserVel = x_hat_out_full(5:8)';
%                         x_hat(:,Epoch_hi) = x_hat_out_full;
%                     elseif (NumFaults > 1)
%                         %% perform multiple fault FDE
%                         disp(sprintf(strcat(ModuleName,'Multiple Faults Detected at time %d'),Epoch_lo));
%                         
%                         %% do not incorporate correction
%                         
%                     else
%                         ExcludeSV = faults;
%                         ExcludeSVPRN = SV_Vec(ExcludeSV);
%                         SVDontUse(ExcludeSVPRN) = 1;  % set don't use flag
%                         
%                         %% Single Fault FDE - 
%                         x_hat(:,Epoch_hi) = x_hat_out(:,ExcludeSV);
%                         UserPos = x_hat_out(1:4,ExcludeSV)';
%                         UserVel = x_hat_out(5:8,ExcludeSV)';
%                         
%                         disp(strcat(ModuleName,sprintf('Fault Detected at Time %d on PRN: %d',Epoch_lo,ExcludeSVPRN)))
%                     end
%                 else
%                         %% if FD Off, use Full Solution
%                         UserPos = x_hat_out_full(1:4)';
%                         UserVel = x_hat_out_full(5:8)';
%                         x_hat(:,Epoch_hi) = x_hat_out_full;
%                 end
%                     
%                 
%                 %%%%% END FAULT DETECTION
%                 
%                 clear PR_Vec PRR_Vec SVPos PR_Vec_minus PRR_Vec_minus PR_csc Beta_ss B_ss lamda_ss SV_Vec TrackTime
%             end % end GPS Meas > 5 multiple model
%             
% 
% 
% 


            % save position history
            %x_out(Epoch_lo,:) = [UserPos, UserVel];

            
            
            
            
%             %%%%%%%% update P_n and V_n %%%%%%%%%
%             [P_xn(Epoch_hi),P_yn(Epoch_hi),P_zn(Epoch_hi)] = ECEF2LLH(x_hat(1:3,Epoch_hi));
%             V_n = Tecef2ned * x_hat(5:7,Epoch_hi);
%             V_xn(Epoch_hi) = V_n(1); V_yn(Epoch_hi) = V_n(2);V_zn(Epoch_hi) = V_n(3);
%             
%                     
%             %% calculate error bounds
%             Pos_var(:,Epoch_hi) = Tecef2ned * P_out_diag(1:3,Epoch_hi);
%             Pos_var(1,Epoch_hi) = sqrt(abs(Pos_var(1,Epoch_hi)));
%             Pos_var(2,Epoch_hi) = sqrt(abs(Pos_var(2,Epoch_hi)));
%             Pos_var(3,Epoch_hi) = sqrt(abs(Pos_var(3,Epoch_hi)));
            
            
            %%%%%%%% LEAST-SQUARES GPS + FDI for comparisson %%%%%%%%%%%%

            % Perform Least-Squares position solution
            %[LSQSolutionVec(Epoch_lo,:) VarLSQSolutionVec(Epoch_lo,:) NumIterations ResidualVector(Epoch_lo,1:SVIndex) ...
            %    M LSQ_Fail(Epoch,:) limit LSQDOP(Epoch_lo,:)] = GARD_LSQ(UserPos,SVIndex,PR_Vec,SVPos);
            % perform RAIM Check
            %[BadGeometry(Epoch_lo), RAIM_ALARM(Epoch_lo), SLOPE_Max(Epoch_lo), r, rd, ARP(Epoch)] = GARD_RAIM(SVIndex,1/15000,7.5,300,ResidualVector(Epoch,:),M);

            % perform RAIM Parity check
            %[BadGeometry, RAIM_ALERT, SLOPE_Max, r, Td, HPL,VPL, FaultySatFDI] = GARDSim_RAIMParity(a, lambda, N,PFalseAlarm,SigmaS,Alert_Limit,ResVec,M, PR);
        end % this is the end of the if statement which decides if this is the first epoch or not
        

    end % end low-speed loop
    
%     Lat_error(Epoch_hi) = (P_xn(Epoch_hi) - pos_truth_llh(2,Epoch_hi))*RM;
%     Long_error(Epoch_hi) = (P_yn(Epoch_hi) - pos_truth_llh(3,Epoch_hi))*RP*cos(P_xn(Epoch_hi));
%     Height_error(Epoch_hi) = (P_zn(Epoch_hi) - pos_truth_llh(4,Epoch_hi));
    
    hgt = Pos_LLH(Epoch_hi,3);
    Lat_error(Epoch_hi) = (Pos_LLH(Epoch_hi,1) - pos_truth_llh(2,Epoch_hi)) * (RM+hgt);
    Long_error(Epoch_hi) = (Pos_LLH(Epoch_hi,2) - pos_truth_llh(3,Epoch_hi)) * ((RP+hgt) * cos(Pos_LLH(Epoch_hi,1)));
    Height_error(Epoch_hi) = (Pos_LLH(Epoch_hi,3) - pos_truth_llh(4,Epoch_hi));

    Vel_error(Epoch_hi,1:3) = Vel_NED(Epoch_hi,1:3) - vel_truth(2:4,Epoch_hi)';
    

    if(mod(Epoch_hi,1000) == 0)
        disp(strcat(ModuleName,sprintf('Completed Epoch %d',Epoch_hi)));
    end

end % end high-speed loop



disp(strcat(ModuleName,'Done!'));



% calculate attitde error
phi_err = (phi_true - phi_q');
theta_err = (theta_true - theta_q');
psi_err = (psi_true - psi_q');

for Epoch=1:NumberEpochsINS
    if(psi_err(Epoch) > pi)
        psi_err(Epoch) = psi_err(Epoch) - 2*pi;
    end
    if(psi_err(Epoch) < -pi)
        psi_err(Epoch) = psi_err(Epoch),+ 2*pi;
    end
end


% plot results
% 
% % % plot velocity
% figure();
% plot(TimeINS,V_xn,'b')
% hold on;
% grid on;
% %plot(TimeGPS(1:Epoch_lo),gps_vel_ned(:,1),'r.');
% plot(vel_truth(1,1:NumberEpochsINS),vel_truth(2,1:NumberEpochsINS),'g');
% xlabel('Time (sec)');
% ylabel('North Velocity (m/s)');
% 
% figure();
% plot(TimeINS,V_yn,'b')
% hold on;
% grid on;
% %plot(TimeGPS(1:Epoch_lo),gps_vel_ned(:,2),'r.');
% plot(vel_truth(1,1:NumberEpochsINS),vel_truth(3,1:NumberEpochsINS),'g')
% xlabel('Time (sec)');
% ylabel('East Velocity (m/s)');
% 
% figure();
% plot(TimeINS,V_zn,'b')
% hold on;
% grid on;
% %plot(TimeGPS(1:Epoch_lo),gps_vel_ned(:,3),'r.');
% plot(vel_truth(1,1:NumberEpochsINS),vel_truth(4,1:NumberEpochsINS),'g');
% xlabel('Time (sec)');
% ylabel('Down Velocity (m/s)');

% %latitude
% figure()
% hold on;
% plot(TimeINS,abs(Lat_error),'b');
% grid on;
% hold on;
% plot(TimeINS,Pos_var(1,:),'r');
% plot(TimeGPS(1:Epoch_lo),HPL,'g');
% xlabel('Simulation Time (sec)');
% ylabel('North-South Latitude Error (m)');
% hold off;
% 
% % longitude
% figure()
% hold on;
% plot(TimeINS,Long_error,'b');
% grid on;
% hold on;
% plot(TimeINS,Pos_var(2,:),'r');
% plot(TimeGPS(1:Epoch_lo),HPL,'g');
% xlabel('Simulation Time (sec)');
% ylabel('East-West Longitude Error (m)');
% hold off;

%% horizontal error with HPL
figure();
hold on;
plot(TimeINS,sqrt(Lat_error.^2 + Long_error.^2),'b');
grid on;
hold on;
plot(TimeINS,2*sqrt(Pos_var(:,1).^2 + Pos_var(:,2).^2),'g');
%plot(TimeGPS(1:Epoch_lo),HPL,'r');
xlabel('Simulation Time (sec)');
ylabel('Horizontal Error (m)');
legend('Horizontal Error','2-\sigma Error Bound','HPL');
%axis([0 3000 0 50]);
hold off;

%% vertical error with VPL
figure();
plot(TimeINS,abs(Height_error),'b');
hold on;
plot(TimeINS,2*sqrt(Pos_var(:,3)),'g');
%plot(TimeGPS(1:Epoch_lo),VPL,'r');
xlabel('Simulation Time (sec)');
ylabel('Vertical Error (m)');
legend('Vertical Error','2-\sigma Error Bound','HPL');
%axis([0 3000 0 50]);
hold off;
grid on;

%% plot satellite usage
figure();
plot(TimeGPS,SatellitesUsed);
xlabel('Test Time (sec)');
ylabel('Satellites Used');


% figure();
% hold on;
% grid on;
% plot(TimeINS,Pos_var(2,:),'r');
% xlabel('Simulation Time (sec)');
% ylabel('East-West Position Variance (m)');
% axis([480 620 0 10]);
% title('Loss of GRAS at T=500s');

% figure();
% plot(TimeINS,Pos_var_nav(1,:),'b');
% hold on;
% plot(TimeINS,Pos_var_tac(1,:),'r');
% xlabel('Simulation Time (sec)');
% ylabel('North-South Position Variance (m)');
% axis([480 620 0 5]);
% title('Loss of GRAS at T=500s');
% legend('High Quality INS','Low Quality INS');
% grid on;

% figure();
% plot(TimeINS,Pos_var_nav(2,:),'b');
% hold on;
% plot(TimeINS,Pos_var_tac(2,:),'r');
% xlabel('Simulation Time (sec)');
% ylabel('East-West Position Variance (m)');
% axis([480 620 0 10]);
% title('Loss of GRAS at T=500s');
% legend('High Quality INS','Low Quality INS');
% grid on;

% % height
% figure()
% hold on;
% plot(TimeINS,Height_error,'b');
% grid on;
% hold on;
% plot(TimeINS,Pos_var(3,:),'r');
% plot(TimeINS,-Pos_var(3,:),'r');
% xlabel('Simulation Time (sec)');
% ylabel('Height Error (m)');
% hold off;

% figure();
% grid on;
% hold on;
% plot(P_yn*r2d,P_xn*r2d,'.');
% %plot(gps_out_llh(:,2),gps_out_llh(:,1),'r.');
% plot(pos_truth_llh(3,:)*r2d,pos_truth_llh(2,:)*r2d,'g');
% hold off;
% xlabel('Longitude (deg)');
% ylabel('Latitude (deg)');
% legend('EKF Estimate','True Trajectory');


% attitude errors
% figure();
% grid on;
% hold on;
% plot(TimeINS,phi_q * r2d,'b');
% plot(TimeINS,phi_uc * r2d,'r');
% plot(att_truth(1,1:NumberEpochsINS),att_truth(2,1:NumberEpochsINS) * r2d,'g');
% hold off;
% xlabel('Time (sec)');
% ylabel('Roll Angle (deg)');
% 
figure();
grid on;
hold on;
plot(TimeINS,phi_err*r2d,'r');
plot(TimeINS,theta_err*r2d,'g');
plot(TimeINS,psi_err*r2d,'b');
hold off;
xlabel('Time (Sec)');
ylabel('Attitude Error (deg)');
legend('Roll','Pitch','Yaw');


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
figure(); hold on; grid on;
plot(TimeINS,Vel_error(:,1),'r');
plot(TimeINS,Vel_error(:,2),'g');
plot(TimeINS,Vel_error(:,3),'b');
xlabel('Simulation Time (seconds)');
ylabel('Velocity Error (m/s)');
title('Velocity Error');
legend('North','East','Down');


%% calculate RMS errors
Height_RMS = sqrt(sum(Height_error.^2)/NumberEpochsINS)
Height_RMS_sub = sqrt(sum(Height_error(3000:NumberEpochsINS).^2)/(NumberEpochsINS-3000))

Lat_RMS = sqrt(sum(Lat_error.^2)/NumberEpochsINS)
Lat_RMS_sub = sqrt(sum(Lat_error(3000:NumberEpochsINS).^2)/(NumberEpochsINS-3000))

Long_RMS = sqrt(sum(Long_error.^2)/NumberEpochsINS)
Long_RMS_sub = sqrt(sum(Long_error(3000:NumberEpochsINS).^2)/(NumberEpochsINS-3000))
