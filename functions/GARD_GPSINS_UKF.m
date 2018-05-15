function results = GARD_GPSINS_UKF(imudata, gps,barodata,magdata, ...
    IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE, Sensor_Params, CompassCalibration, StartTime, StopTime, InitialPosition, InitialVelocity,InitialAttitude)

% Unscented Kalman Filter Implementation for integrated GPS-INS
%
% Implements an integrated GPS-Inertial position solution
%
% Written by Duncan Greer
%
% $Id: GARDSim_GPSINS_UKF.m 1879 2008-07-15 05:20:21Z n2523710 $
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


ModuleName = '[GARDSim_GPSINS_UKF] ';





% initialise epochs




gps_dt = 1/GPS_RATE;
ins_dt = 1/IMU_RATE;

NumberEpochs_lo = (StopTime - StartTime ) * GPS_RATE;
NumberEpochs_hi = NumberEpochs_lo * IMU_RATE;

NumberStates = 18;

% process noise states - 6 for IMU sensors, 2 for gps clock
ProcessNoiseStates = 8;  

% measurement noise states - pseudorange and pseuorange rate for each
% satellite
MeasurementNoiseStates=16;


% setup constants

WGS84Constants;
d2r = pi/180;
r2d = 180/pi;


% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));

% gravity model - Kayton eqn 2.6
%g = -9.79;  % m/s/s
g = -0.01 * (978.049 * (1 + .00529 * sin(InitialPosition(1))^2));


% for the purpose of this simulation, the number of measurements is
% fixed at 16 = 2*8 - the best 8 measurements based on GDOP(??) should be chosen, plus there are PRR measurments for each SV.
NumberMeasurements = 16;
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
x_gyro_beta = Sensor_Params.gyro_beta(1);
y_gyro_beta = Sensor_Params.gyro_beta(2);
z_gyro_beta = Sensor_Params.gyro_beta(3);
x_accel_beta = Sensor_Params.accel_beta(1);
y_accel_beta = Sensor_Params.accel_beta(2);
z_accel_beta = Sensor_Params.accel_beta(3);

x_gyro_Q = Sensor_Params.gyro_Q(1);
y_gyro_Q = Sensor_Params.gyro_Q(2);
z_gyro_Q = Sensor_Params.gyro_Q(3);
x_accel_Q = Sensor_Params.accel_Q(1);
y_accel_Q = Sensor_Params.accel_Q(2);
z_accel_Q = Sensor_Params.accel_Q(3);


% gps clock bias and frequency error PSD
Sc = 1*0.14;
Sf = 1*0.0359;

% initialise INS process noise
UKF_Q = zeros(ProcessNoiseStates,ProcessNoiseStates);

UKF_Q(1,1) = 2*x_accel_beta*x_accel_Q;
UKF_Q(2,2) = 2*y_accel_beta*y_accel_Q;
UKF_Q(3,3) = 2*z_accel_beta*z_accel_Q;
UKF_Q(4,4) = 2*x_gyro_beta*x_gyro_Q;
UKF_Q(5,5) = 2*y_gyro_beta*y_gyro_Q;
UKF_Q(6,6) = 2*z_gyro_beta*z_gyro_Q;

%UKF_Q = UKF_Q * ins_dt;
% 
% UKF_Q(7,7) = (Sc)*ins_dt + Sf*(ins_dt^3)/3;
% UKF_Q(8,8) = Sf*ins_dt;
% UKF_Q(7,8) = Sf*(ins_dt^2)/2;
% UKF_Q(8,7) = Sf*(ins_dt^2)/2;

UKF_Q(7,7) = (Sc)*gps_dt + Sf*(gps_dt^3)/3;
UKF_Q(8,8) = Sf*gps_dt;
UKF_Q(7,8) = Sf*(gps_dt^2)/2;
UKF_Q(8,7) = Sf*(gps_dt^2)/2;



% initialise measurement niose - pseudorange, then doppler 
sigma_pr_gnd = 5.0;
sigma_prr = 0.5;

GRASOn = 0;

if GRASOn == 0
    sigma_iono = 2.0;
else
    sigma_iono = 0.1;
end

sigma_pr = sqrt(sigma_pr_gnd^2 + sigma_iono^2);

UKF_R = zeros(NumberMeasurements,NumberMeasurements);
UKF_R(1:8,1:8) = sigma_pr^2*eye(8,8);
UKF_R(9:16,9:16) = sigma_prr^2*eye(8,8);


Na = NumberStates+ProcessNoiseStates+MeasurementNoiseStates;

% UKF scaling parameters
alpha = 1; % range: 1e-3 < alpha <= 1
beta = 2;  % 2 is optimal for gaussian priors
kapa = 0;  %

LAMBDA = alpha^2 * (Na+kapa) - Na;

% nana is a variable which describes the convergence speed of the
% quaternion to unity
nana=0.5;


NumberEpochsGPS = NumberEpochs_lo;

TimeINS = 0:ins_dt:NumberEpochsGPS;
TimeGPS = 0:gps_dt:NumberEpochsGPS-1;

NumberSVs = 32;
SVDontUse = zeros(1,NumberSVs);
SV_TrackTime = zeros(NumberSVs,1);


%dt = gps_dt;

NumberEpochsINS = NumberEpochs_hi;


% initialise epochs
Epoch_lo = 0;
Epoch_hi = 0;


FDEnabled = 0;

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



    lat = InitialPosition(1);
    lon = InitialPosition(2);
    Tecef2ned2 = T_ECEF2NED(lat,lon);
    Tecef2ned2(4,1:3) = [0 0 0];
    Tecef2ned2(1:4,4) = [0 0 0 1];


NumberEpochs_lo = (StopTime - StartTime ) * GPS_RATE;
NumberEpochs_hi = NumberEpochs_lo * IMU_RATE;


gpstime = zeros(size(gps.RangeData,2),3);
for i=1:size(gps.RangeData,2)
   gpstime(i,1) = gps.RangeData(i).rtTimeStamp*1e-9;
   gpstime(i,2) = gps.RangeData(i).GPSWeek;
   gpstime(i,3) = gps.RangeData(i).GPSSec;
end

% initialise imudataindex_i - this value is updated every gps update
imudataindex_i = find(imudata(:,1) < gpstime(60,1), 1,'last');
imudataindex = imudataindex_i;

%% start processing
for Epoch_hi = 1:NumberEpochs_hi
    
    imudataindex_save = imudataindex;
    
    if mod(Epoch_hi, IMU_RATE/GPS_RATE) == 0
        imudataindex = imudataindex_i + IMU_RATE/GPS_RATE;
    else
        imudataindex = imudataindex_i + mod(Epoch_hi, IMU_RATE/GPS_RATE);
    end
    
    imudatatime(Epoch_hi,1) = imudata(imudataindex,1);
    imudatatime(Epoch_hi,2) = imudataindex;
    imudatatime(Epoch_hi,3) = imudataindex_i;
    
    if imudataindex - imudataindex_save > 1
       disp(sprintf('imu data index jump was greater than 1: %d %d', imudataindex, imudataindex_save)); 
    end
    
   if Epoch_hi == 1
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
        
        for i=1:8
            xs_Sub(:,:,i) = xs;
        end
   else
                % get sensor measurements
          if(Epoch_lo < 2)
             omega_x_0 = imudata(imudataindex,8);
             omega_y_0 = imudata(imudataindex,9);
             omega_z_0 = imudata(imudataindex,10);
 
             A_xb_0 = imudata(imudataindex,5);
             A_yb_0 = imudata(imudataindex,6);
             A_zb_0 = imudata(imudataindex,7);
         else
            omega_x_0 = imudata(imudataindex,8) - GyroBias(Epoch_lo,1);
            omega_y_0 = imudata(imudataindex,9) - GyroBias(Epoch_lo,2);
            omega_z_0 = imudata(imudataindex,10) - GyroBias(Epoch_lo,3);

            A_xb_0 = imudata(imudataindex,5) - AccelBias(Epoch_lo,1);
            A_yb_0 = imudata(imudataindex,6) - AccelBias(Epoch_lo,2);
            A_zb_0 = imudata(imudataindex,7) - AccelBias(Epoch_lo,3);
          end
         
          t_INS(Epoch_hi,:) = imudata(imudataindex,1);
          
          % save sensor inputs
          omega_b(Epoch_hi,:) = [omega_x_0,omega_y_0,omega_z_0];
          A_b(Epoch_hi,:) = [A_xb_0,A_yb_0,A_zb_0];
        
        % if we have just come out of a GPS loop, recalculate the sigma
        % points
        if mod(Epoch_hi-1,IMU_RATE/GPS_RATE) == 0
             %disp(sprintf('Resetting x_hat_kminus at Time %f',Epoch_hi/100));            
            % generate augmented state vector

            [xs,Ws] = GARD_GenerateSigmaPoints(UKF_x_hat_kminus,UKF_Px_kminus,UKF_Q,UKF_R,alpha,beta,kapa);
            xs_0 = xs(:,1);
            xs_i = xs(:,2:size(xs,2));
            W_0_m = Ws(1);
            W_0_c = Ws(2);
            W_i_m = Ws(3);
            W_i_c = Ws(4);
            if exist('NumberSubFilters','var')
                for i=1:NumberSubFilters
                    xs_Sub(:,:,i) = GARD_GenerateSigmaPoints(UKF_x_hat_kminus_Sub(:,i),UKF_Px_kminus_Sub(:,:,i),UKF_Q,UKF_R,alpha,beta,kapa);
                end
            end
            
        end

        
        A_b_in = [A_xb_0;A_yb_0;A_zb_0];
        Omega_b_in = [omega_x_0; omega_y_0;omega_z_0];
        
        [xs_0_k,xs_i_k] = GARD_PropagateSigmaPointsINS(Na,xs_0,xs_i,A_b_in,Omega_b_in,ins_dt,NumberStates);

        
        
        
        if exist('NumberSubFilters','var')
            for i=1:NumberSubFilters
        
                A_b_in_Sub = imudata(imudataindex,5:7) - UKF_x_hat_kminus_Sub(11:13,i)';
                Omega_b_in_Sub = imudata(imudataindex,8:10) - UKF_x_hat_kminus_Sub(14:16,i)';
                [xs_0_k_Sub(:,i),xs_i_k_Sub(:,:,i)] = GARD_PropagateSigmaPointsINS(Na,xs_Sub(:,1,i),xs_Sub(:,2:size(xs_Sub,2),i),A_b_in_Sub,Omega_b_in_Sub,ins_dt,NumberStates);
                xs_Sub(:,:,i) = [xs_0_k_Sub(:,i),xs_i_k_Sub(:,:,i)];
            end
        end
            
        % update xs_i
        xs_0 = xs_0_k;
        xs_i = xs_i_k;
        
        % calculate the predicted 
        % mean and state covariance
         xs = [xs_0_k,xs_i_k];
         [UKF_x_hat_kminus,UKF_Px_kminus] = GARD_GenerateMeanFromSigmaPoints(xs,Ws,NumberStates);

        if exist('NumberSubFilters','var')
            for i=1:NumberSubFilters
                [UKF_x_hat_kminus_Sub(:,i),UKF_Px_kminus_Sub(:,:,i)] = GARD_GenerateMeanFromSigmaPoints(xs_Sub(:,:,i),Ws,NumberStates);
            end
        end
        
         
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
    if(mod(Epoch_hi,IMU_RATE/GPS_RATE) == 0)

        Epoch_lo = Epoch_lo + 1;
        
            
            % find the gps data index of this epoch
            gpsdataindex = find(gpstime(:,3) == gpstime(60,3) + StartTime + (Epoch_lo-1)*GPS_RATE,1);
            
            % update the current imu data index
            imudataindex_i = find(imudata(:,1) < gpstime(gpsdataindex,1), 1,'last');
            
           
%             
            %if(Epoch_lo == 60)  % enable fault detection
                FDEnabled = 1;
            %end

            if Epoch_lo == 1
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
            NumberGPSMeasurements = gps.RangeData(gpsdataindex).lNumberObservations;
            for Obs=1:NumberGPSMeasurements
                SV = gps.RangeData(gpsdataindex).Obs(Obs).usPRN;
                if( SVDontUse(SV) == 0)  %%SV ~= ExcludeSVPRN
                        % add to PR vector
                        SVIndex = SVIndex + 1;

                       SV_TrackTime(SV) = SV_TrackTime(SV) + 1;

                        SV_Vec(SVIndex) = SV;

                        PR_Vec(SVIndex) = gps.RangeData(gpsdataindex).Obs(Obs).dPseudorange;
                        PRR_Vec(SVIndex) = gps.RangeData(gpsdataindex).Obs(Obs).fDoppler * -L1_Wavelength;



                        [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV)] = ...
                           GPSOrbitPropagator(gps.RangeData(gpsdataindex).GPSWeek, gps.RangeData(gpsdataindex).GPSSec - PR_Vec(SVIndex)/Speedoflight, SV, gps.SV_Ephemeris, 7500);

                       if ValidPosData(Epoch_lo,SV) == 0
                           SVIndex = SVIndex - 1;
                           continue;
                       end

                         [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
                          SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(Epoch_lo,SV)] = ...
                          GPSOrbitPropagatorVelocities(gps.RangeData(gpsdataindex).GPSWeek, gps.RangeData(gpsdataindex).GPSSec-PR_Vec(SVIndex)/Speedoflight, SV, gps.SV_Ephemeris,7500);

                      if ValidVelData(Epoch_lo,SV) == 0
                          SVIndex = SVIndex - 1;
                          continue;
                      end

                        SVPos(SVIndex,4) = SVPos(SVIndex,4) * Speedoflight;
                        SVVel(SVIndex,4) = SVVel(SVIndex,4) * Speedoflight;
                        SVAcc(SVIndex,4) = SVAcc(SVIndex,4) * Speedoflight;

                        [SV_Azimuth(Epoch_lo,SV), SV_Elevation(Epoch_lo,SV)] = AzEl(UserPos(1:3), SVPos(SVIndex,1:3));

                        % calculate the iono delay correction - single frequency user
                        % model from ICD 200

                        ALPHA = [gps.IONUTCData(1).a0 gps.IONUTCData(1).a1 gps.IONUTCData(1).a2 gps.IONUTCData(1).a3];
                        BETA = [gps.IONUTCData(1).b0 gps.IONUTCData(1).b1 gps.IONUTCData(1).b2 gps.IONUTCData(1).b3];

                        ionodelay = ionomodel(gps.RangeData(gpsdataindex).GPSSec, UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
                        PR_Vec(SVIndex) = PR_Vec(SVIndex) - ionodelay;


                        %if ~exist('x_save_llh','var')
                        %    TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),InitialPosition(3));
                        %else    
                            TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),UKF_x_hat_kminus(3));
                        %end
                        PR_Vec(SVIndex) = PR_Vec(SVIndex) - TropoDelay(Epoch_lo,SV);

                        % calculate hte earth rotation correction as per Kayton pg 228
                        % eq 5.67


                        PR_Vec_raw(SVIndex) = PR_Vec(SVIndex);  % save a raw (uncorrected copy) of the PR vector for use in the LSQ algorithm later.
                        PRR_Vec_raw(SVIndex) = PRR_Vec(SVIndex);
                        delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) *UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
                        PR_Vec(SVIndex) = PR_Vec(SVIndex) + delta_pr_omegaedot + SVPos(SVIndex,4);
                        PRR_Vec(SVIndex) = PRR_Vec(SVIndex) + SVVel(SVIndex,4);



                end
                
          

            end
            
            
        NumberGPSMeasurements = SVIndex;

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
        UKF_y_k = [PR_Vec_raw(1:8)';PRR_Vec_raw(1:8)'];
        
%         PR_Vec_minus_0 = zeros(NumberGPSMeasurements,1);
%         PRR_Vec_minus_0 = zeros(NumberGPSMeasurements,1);
%         
%         PR_Vec_minus_i = zeros(NumberGPSMeasurements,2*Na);
%         PRR_Vec_minus_i = zeros(NumberGPSMeasurements,2*Na);
%         
%         Relative_Velocity = zeros(NumberGPSMeasurements,1);
%         
%         
%         % get the a-proiri measuremnt prediction (y_minus)
%         for k = 1:NumberGPSMeasurements
%             Tecef2ned= T_ECEF2NED(xs_0_k(1),xs_0_k(2)); 
%             Tned2ecef = Tecef2ned';
%             % find apriori estimate of pseudorange for each sigma point
%             
%             % zero-th sigma point
%             %Get User position in ECEF
%             UserPos = LLH2ECEF(xs_0_k(1),xs_0_k(2),xs_0_k(3));
%             UserPos(4) = xs_0_k(17);
%             
%             UserVel =  Tned2ecef * xs_0_k(4:6);
%             UserVel(4) = xs_0_k(18);
%             
%             geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
%             geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + ...
%                              (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + ...
%                              (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
%             delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) *UserPos(2) - SVPos(k,2) * UserPos(1));
%             
%             %% calculate measurement predicition
%             % pseudorange prediction
%             PR_Vec_minus_0(k) = geo_range_to_sat + UserPos(4) - delta_pr_omegaedot(k) - SVPos(k,4) + xs_0(NumberStates+ProcessNoiseStates+k);
%             %predicted relative velocity of sv and receiver
%             Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
%             PRR_Vec_minus_0(k) = Relative_Velocity(k) + UserVel(4) + xs_0(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+k) - SVVel(k,4);
%             
%             for i=1:2*Na
%                 Tecef2ned= T_ECEF2NED(xs_i_k(1,i),xs_i_k(2,i));
%                 Tned2ecef = Tecef2ned';
%                 
%                % Get User position in ECEF
%                 UserPos = LLH2ECEF(xs_i_k(1,i),xs_i_k(2,i),xs_i_k(3,i));
%                 UserPos(4) = xs_i_k(17,i);
%                 
%                 UserVel =  Tned2ecef * xs_i_k(4:6,i);
%                 UserVel(4) = xs_i_k(18,i);
%                 
%                 for m = 1:3
%                     ele(m) =  SVPos(k,m) - UserPos(m);
%                 end
% 
%                 geo_range_to_sat =  norm(ele);
%         
%                 %geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
%                 geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + ...
%                                  (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + ...
%                                  (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
%                 delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) *UserPos(2) - SVPos(k,2) * UserPos(1));
%                 
%                 PR_Vec_minus_i(k,i) = geo_range_to_sat + UserPos(4) - delta_pr_omegaedot(k)  - SVPos(k,4) + xs_i(NumberStates+ProcessNoiseStates+k,i);
% 
%                 %predicted relative velocity of sv and receiver
%                 Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
%                 PRR_Vec_minus_i(k,i) = Relative_Velocity(k) + UserVel(4) - SVVel(k,4) + xs_i(NumberStates+ProcessNoiseStates+NumberGPSMeasurements+k,i) ;
%             end % for i=1:2*Na
% 
% 
% 
%         end  % for k = 1:NumberGPSMeasurements
%         
%         ys_kminus_0 = [PR_Vec_minus_0; PRR_Vec_minus_0];
%         ys_kminus_i = [PR_Vec_minus_i(:,:); PRR_Vec_minus_i(:,:)];
% 

        [ys_kminus_0, ys_kminus_i] = GARD_GenerateMeasurementPredictionSigmaPoints(Na,xs_0_k,xs_i_k,NumberGPSMeasurements,NumberStates,ProcessNoiseStates,SVPos,SVVel);


        % calculate measurement covariance
         ys = [ys_kminus_0,ys_kminus_i];
         [UKF_y_hat_kminus,UKF_Py_kminus,UKF_Pxy_kminus] = GARD_GenerateMeasurementMeanFromSigmaPoints(xs,ys,Ws,NumberStates,NumberMeasurements);

        
        % calculate kalman gain
        UKF_K_k = UKF_Pxy_kminus * inv(UKF_Py_kminus);

        
        % apply correction
        UKF_z = (UKF_y_k - UKF_y_hat_kminus);
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

         
         
        %save results
        t_save(Epoch_hi) = Epoch_hi/IMU_RATE/GPS_RATE;
        UKF_x_hat_save(Epoch_hi,:) = UKF_x_hat_kplus;
        UKF_P_save(Epoch_hi,:) = diag(UKF_Px_kplus);
    
        % save the post-update residual (shoudl be whiteish)
        [PR_Vec_plus, PRR_Vec_plus] = GARD_GetMeasurementPrediction(UKF_x_hat_kplus,NumberGPSMeasurements,SVPos,SVVel);
        UKF_y_hat_kplus = [PR_Vec_plus; PRR_Vec_plus];
        UKF_z_post_save(Epoch_lo,:) = UKF_y_k - UKF_y_hat_kplus;
        
        % copy updated state to high-speed loop state for next round
        UKF_x_hat_kminus = UKF_x_hat_kplus;
        UKF_Px_kminus = UKF_Px_kplus;
        
        GyroBias(Epoch_lo,:) = [UKF_x_hat_kplus(14),UKF_x_hat_kplus(15),UKF_x_hat_kplus(16)];
        AccelBias(Epoch_lo,:) = [UKF_x_hat_kplus(11),UKF_x_hat_kplus(12),UKF_x_hat_kplus(13)];

        
%% calculate sub solutions

    if FDEnabled == 1
        NumberSubMeasurements = NumberGPSMeasurements-1;
        NumberSubFilters = nchoosek(NumberGPSMeasurements,NumberSubMeasurements);
        
        % if we're on the first run, initialise xs_Sub from xs
        if ~exist('xs_Sub','var')
           for i=1:NumberSubFilters
                xs_Sub(:,:,i) = xs;
           end
        end

        for SubSolution = 1:NumberSubFilters

%% ===== start of original ======
            
%             PR_SubVec(1:SubSolution-1) = PR_Vec(1:SubSolution-1);
%             PR_SubVec(SubSolution:NumberSubMeasurements) = PR_Vec(SubSolution+1:NumberGPSMeasurements);
% 
%             PRR_SubVec(1:SubSolution-1) = PRR_Vec(1:SubSolution-1);
%             PRR_SubVec(SubSolution:NumberSubMeasurements) = PRR_Vec(SubSolution+1:NumberGPSMeasurements);
%             
%             % formulate the measurment vector
%             UKF_y_k_Sub = [PR_SubVec(1:NumberSubMeasurements)';PRR_SubVec(1:NumberSubMeasurements)'];
% 
%             
% %             PR_Vec_minus_0_Sub(1:SubSolution-1) = PR_Vec_minus_0(1:SubSolution-1);
% %             PR_Vec_minus_0_Sub(SubSolution:NumberSubMeasurements) = PR_Vec_minus_0(SubSolution+1:NumberGPSMeasurements);
% %             
% %             PRR_Vec_minus_0_Sub(1:SubSolution-1) = PRR_Vec_minus_0(1:SubSolution-1);
% %             PRR_Vec_minus_0_Sub(SubSolution:NumberSubMeasurements) = PRR_Vec_minus_0(SubSolution+1:NumberGPSMeasurements);
% %             
% %             PR_Vec_minus_i_Sub(1:SubSolution-1,:) = PR_Vec_minus_i(1:SubSolution-1,:);
% %             PR_Vec_minus_i_Sub(SubSolution:NumberSubMeasurements,:) = PR_Vec_minus_i(SubSolution+1:NumberGPSMeasurements,:);
% %             
% %             PRR_Vec_minus_i_Sub(1:SubSolution-1,:) = PRR_Vec_minus_i(1:SubSolution-1,:);
% %             PRR_Vec_minus_i_Sub(SubSolution:NumberSubMeasurements,:) = PRR_Vec_minus_i(SubSolution+1:NumberGPSMeasurements,:);
% %              
% %             ys_kminus_0_Sub = [PR_Vec_minus_0_Sub'; PRR_Vec_minus_0_Sub'];
% %             ys_kminus_i_Sub = [PR_Vec_minus_i_Sub(:,:); PRR_Vec_minus_i_Sub(:,:)];
% 
%             SVPos_Sub(1:SubSolution-1,:) = SVPos(1:SubSolution-1,:); 
%             SVPos_Sub(SubSolution:NumberSubMeasurements,:) = SVPos(SubSolution+1:NumberGPSMeasurements,:);
%             
%             SVVel_Sub(1:SubSolution-1,:) = SVVel(1:SubSolution-1,:); 
%             SVVel_Sub(SubSolution:NumberSubMeasurements,:) = SVVel(SubSolution+1:NumberGPSMeasurements,:);
%             
%             %SVVel_Sub = 
%             
%             [ys_kminus_0_Sub, ys_kminus_i_Sub] = GARD_GenerateMeasurementPredictionSigmaPoints(Na,xs_Sub(:,1,SubSolution),xs_Sub(:,2:size(xs_Sub,2),SubSolution),NumberGPSMeasurements,NumberStates,ProcessNoiseStates,SVPos_Sub,SVVel_Sub,NumberSubMeasurements);
% 
%             ys_Sub = [ys_kminus_0_Sub,ys_kminus_i_Sub];
%             
%             
%             [UKF_y_hat_kminus_Sub,UKF_Py_kminus_Sub,UKF_Pxy_kminus_Sub] = GARD_GenerateMeasurementMeanFromSigmaPoints(xs_Sub(:,:,SubSolution),ys_Sub,Ws,NumberStates,NumberSubMeasurements);
%             
%             % calculate kalman gain
%             UKF_K_k_Sub = UKF_Pxy_kminus_Sub * inv(UKF_Py_kminus_Sub);
% 
%         
%             % apply correction
%             UKF_z_Sub = (UKF_y_k_Sub - UKF_y_hat_kminus_Sub);
%             
%            [UKF_x_hat_kminus_Sub(:,SubSolution),UKF_Px_kminus_Sub(:,:,SubSolution)] = GARD_GenerateMeanFromSigmaPoints(xs_Sub(:,:,SubSolution),Ws,NumberStates);
%             
%             UKF_x_hat_kplus_Sub(:,SubSolution) = UKF_x_hat_kminus_Sub(:,SubSolution) + UKF_K_k_Sub * UKF_z_Sub;
% 
%             UKF_Px_kplus_Sub(:,:,SubSolution) = UKF_Px_kminus_Sub(:,:,SubSolution) - UKF_K_k_Sub * UKF_Py_kminus_Sub * UKF_K_k_Sub';

%% ==== END OF ORIGINAL ======
% new stuff
            
            [ys_kminus_0_Sub, ys_kminus_i_Sub] = GARD_GenerateMeasurementPredictionSigmaPoints(Na,xs_Sub(:,1,SubSolution),xs_Sub(:,2:size(xs_Sub,2),SubSolution),NumberGPSMeasurements,NumberStates,ProcessNoiseStates,SVPos,SVVel);

            ys_Sub = [ys_kminus_0_Sub,ys_kminus_i_Sub];
            
            
            [UKF_y_hat_kminus_Sub,UKF_Py_kminus_Sub,UKF_Pxy_kminus_Sub] = GARD_GenerateMeasurementMeanFromSigmaPoints(xs_Sub(:,:,SubSolution),ys_Sub,Ws,NumberStates,NumberMeasurements);
            
            % calculate kalman gain
            UKF_K_k_Sub = real(UKF_Pxy_kminus_Sub * inv(real(UKF_Py_kminus_Sub)));

            
            % zero the Kalman gain for hte satellites we want to ignore
            UKF_K_k_Sub(:,SubSolution) = zeros(NumberStates,1);
            UKF_K_k_Sub(:,NumberSubFilters+SubSolution) = zeros(NumberStates,1);
            
        
            % apply correction
            UKF_z_Sub = (UKF_y_k - real(UKF_y_hat_kminus_Sub));
            
           [UKF_x_hat_kminus_Sub(:,SubSolution),UKF_Px_kminus_Sub(:,:,SubSolution)] = GARD_GenerateMeanFromSigmaPoints(xs_Sub(:,:,SubSolution),Ws,NumberStates);
            
            UKF_x_hat_kplus_Sub(:,SubSolution) = real(UKF_x_hat_kminus_Sub(:,SubSolution) + UKF_K_k_Sub * UKF_z_Sub);

            UKF_Px_kplus_Sub(:,:,SubSolution) = real(UKF_Px_kminus_Sub(:,:,SubSolution) - UKF_K_k_Sub * UKF_Py_kminus_Sub * UKF_K_k_Sub');
            
            




%% ==== END ==== 
            % save for next iteration
            UKF_x_hat_kminus_Sub(:,SubSolution) = real(UKF_x_hat_kplus_Sub(:,SubSolution));
            UKF_Px_kminus_Sub(:,:,SubSolution) = real(UKF_Px_kplus_Sub(:,:,SubSolution));
            
            
            % save results for output
            UKF_x_hat_Sub_Save(Epoch_hi,:,SubSolution) = real(UKF_x_hat_kplus_Sub(:,SubSolution));
            UKF_P_Sub_Save(Epoch_hi,:,SubSolution) = diag(real(UKF_Px_kplus_Sub(:,:,SubSolution)));
        end % end sub solutions

        
        
        for j = 1:NumberSubFilters
            % note that only the position estimates are used - not
            % clock or velocity
            P_ned = UKF_Px_kplus(1:3,1:3);
            P_ned_Sub = UKF_Px_kplus_Sub(1:3,1:3,j);

            NED_ss = (UKF_x_hat_kplus(1:3) - UKF_x_hat_kplus_Sub(1:3,j));
            UKF_Beta_ss_H(:,j) = NED_ss(1:2);
            UKF_B_ss_H(:,:,j) = abs(abs(P_ned_Sub(1:2,1:2)) - abs(P_ned(1:2,1:2)));

            UKF_lambda_ss_H(Epoch_lo,j) = UKF_Beta_ss_H(:,j)' * pinv(UKF_B_ss_H(:,:,j)) * UKF_Beta_ss_H(:,j);

            UKF_B_lambda_H = abs(eigs(UKF_B_ss_H(:,:,j))); % calculates eigenvalues of separation covariance
            %UKF_TD_H(Epoch_lo,j) = sqrt(max(UKF_B_lambda_H)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));
            UKF_TD_H(Epoch_lo,j) = chi2inv(1-PFalseAlarm/NumberGPSMeasurements,NumberGPSMeasurements-4);

            if(sqrt(abs(UKF_lambda_ss_H(Epoch_lo,j))) > UKF_TD_H(Epoch_lo,j))
                disp(sprintf('%s H-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
            end

            UKF_Beta_ss_V(:,j) = NED_ss(3);
            UKF_B_ss_V(:,:,j) = abs(P_ned_Sub(3,3) - P_ned(3,3));

            UKF_lambda_ss_V(Epoch_lo,j) = UKF_Beta_ss_V(:,j)' * pinv(UKF_B_ss_V(:,:,j)) * UKF_Beta_ss_V(:,j);

            UKF_B_lambda_V = abs(eigs(UKF_B_ss_V(:,:,j)));
            
            %UKF_TD_V(Epoch_lo,j) = sqrt(max(UKF_B_lambda_V)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));
            UKF_TD_V(Epoch_lo,j) = chi2inv(1-PFalseAlarm/NumberGPSMeasurements,NumberGPSMeasurements-4);  % actually, this will be the same as horizontal.

            if(sqrt(abs(UKF_lambda_ss_V(Epoch_lo,j))) > UKF_TD_V(Epoch_lo,j))
                disp(sprintf('%s V-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
            end

            %% this HPL represents the fault-free (H0) hypothesis. 
            %Sigma_MD_H = abs(norminv(P_MD_H,0,1));
            %Sigma_MD_V = abs(norminv(P_MD_V,0,1));
            %UKF_HPL_sub(Epoch_lo,j) = Sigma_MD_H * sqrt(P_ned(1,1) + P_ned(2,2)) + UKF_TD_H(Epoch_lo,j);
            %UKF_VPL_sub(Epoch_lo,j) = Sigma_MD_V * sqrt(P_ned(3,3)) + UKF_TD_V(Epoch_lo,j);

            UKF_P_lambda_H = abs(eigs(P_ned(1:2,1:2)));
            UKF_P_lambda_V = abs(eigs(P_ned(3,3)));
            
            UKF_K_H0_H = CEP_TableLookup_PmdH(min(UKF_P_lambda_H) / max(UKF_P_lambda_H));
            UKF_K_H0_V = CEP_TableLookup_Pmd_V(1); % only have one eigen value for vert... 
            
            UKF_HPL_H0_sub(Epoch_lo,j) = UKF_K_H0_H * sqrt(max(UKF_P_lambda_H));
            UKF_VPL_H0_sub(Epoch_lo,j) = UKF_K_H0_V * sqrt(max(UKF_P_lambda_V));
            
            %% todo - calculate the Fault-in progress (h1) hypothesis based
            %% HPL_H1.

            % HPL_H1 = max (HPE_B, HPE_NP, HPE_NB)
            Pbias_H = sqrt(RAIM_lambda_H(NumberGPSMeasurements-4));
            Pbias_V = sqrt(RAIM_lambda_V(NumberGPSMeasurements-4));
            
            HPE_B = Pbias_H * sqrt(max(UKF_B_lambda_H));
            eigs_sub = abs(eigs(P_ned_Sub(1:2,1:2)));
            K_NP_H = CEP_TableLookup_PmdH(min(eigs_sub)/max(eigs_sub));
            HPE_NP = K_NP_H * sqrt(max(eigs_sub));
            
            K_NB_H = CEP_TableLookup_PmdH(min(UKF_B_lambda_H) / max(UKF_B_lambda_H));
            HPE_NB = K_NB_H * sqrt(max(UKF_B_lambda_H));
            
            
            UKF_HPE(j) = HPE_B + HPE_NP - HPE_NB;

            VPE_B = Pbias_V * sqrt(max(UKF_B_lambda_V));
            
            eigs_sub = abs(eigs(P_ned_Sub(3,3)));
            K_NP_V = CEP_TableLookup_Pmd_V(1);
            VPE_NP = K_NP_V * sqrt(max(eigs_sub));
            
            K_NB_V = CEP_TableLookup_Pmd_V(1);
            VPE_NB = K_NB_V * sqrt(max(UKF_B_lambda_V));
            
            UKF_VPE(j) = VPE_B + VPE_NP - VPE_NB;
        end

        UKF_HPL_H1(Epoch_lo) = max(UKF_HPE);
        UKF_VPL_H1(Epoch_lo) = max(UKF_VPE);
        
        UKF_HPL_H0(Epoch_lo) = max(UKF_HPL_H0_sub(Epoch_lo,:));
        UKF_VPL_H0(Epoch_lo) = max(UKF_VPL_H0_sub(Epoch_lo,:));
        
        UKF_HPL(Epoch_lo) = max(UKF_HPL_H0(Epoch_lo),UKF_HPL_H1(Epoch_lo));
        UKF_VPL(Epoch_lo) = max(UKF_VPL_H0(Epoch_lo),UKF_VPL_H1(Epoch_lo));
    end

                    
    else
        
        % not a low-speed loop
        t_save(Epoch_hi) = Epoch_hi/100;
        UKF_x_hat_save(Epoch_hi,:) = UKF_x_hat_kminus;
        UKF_P_save(Epoch_hi,:) = diag(UKF_Px_kminus);

        if exist('NumberSubFilters','var')
            for SubSolution=1:NumberSubFilters
                UKF_x_hat_Sub_Save(Epoch_hi,:,SubSolution) = UKF_x_hat_kminus_Sub(:,SubSolution);
                UKF_P_Sub_Save(Epoch_hi,:,SubSolution) = diag(UKF_Px_kminus_Sub(:,:,SubSolution));
            end
        end
        
        
    end% end low speed loop

    
    

    if(mod(Epoch_hi,IMU_RATE/GPS_RATE) == 0)
        logline = sprintf('%s Completed Epoch %d of %d (%3.1f percent)\n',ModuleName,Epoch_lo,NumberEpochs_lo, 100*(Epoch_lo/NumberEpochs_lo));
        disp(logline);
        
        
    end
    

end % end high-speed loop

        

disp('[GARDSim_GPSINS_UKF] Done!');

results.UKF_x_hat_save = UKF_x_hat_save;
results.UKF_P_save = UKF_P_save;

if exist('UKF_HPL','var')
    results.UKF_HPL = UKF_HPL;
    results.UKF_VPL = UKF_VPL;
end
results.GyroBias = GyroBias;
results.AccelBias = AccelBias;

results.imudatatime = imudatatime;
results.phi_q = phi_q;
results.theta_q = theta_q;
results.psi_q = psi_q;


results.UKF_TD_V = UKF_TD_V;
results.UKF_TD_H = UKF_TD_H;
results.UKF_lambda_ss_V = UKF_lambda_ss_V;
results.UKF_lambda_ss_H = UKF_lambda_ss_H;

results.UKF_x_hat_Sub_Save = UKF_x_hat_Sub_Save;
results.UKF_P_Sub_Save = UKF_P_Sub_Save;
