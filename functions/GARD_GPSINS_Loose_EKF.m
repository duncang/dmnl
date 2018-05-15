function results = GARD_GPSINS_Loose_EKF(imudata, gpsdata,barodata,magdata, GPSLeverArm, ...
    IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE, Sensor_Params, CompassCalibration, ...
    GPSPos_E,GPSPos_EE,GPSVel_NED,GPSVel_EE)
% function result = GARD_GPSINS_Loose_EKF(imudata,gpsdata,barodata,magdata);
% Loosely Coupled GPS-INS using an EKF
% Written by Duncan Greer
% $Id: GARD_GPSINS_Loose_EKF.m 4095 2010-10-30 01:38:06Z greerd $
%
% INPUTS
% ======
% for i measurements of each type:
% imudata: ix12 array of IMU measurements
% gpsdata: ix19 array of GPS positions and velocities
% barodata: ix2 array of Barometric altitude 
% magdata: ix4 array of magnetometer measurements
% IMU_RATE: rate of IMU data in Hz (e.g. 200)
% GPS_RATE: rate of GPS data in Hz
% MAG_RATE: rate of Magnetometer data in Hz
% BARO_RATE: rate of barometric altitude data in Hz
%
% gpsdata:
% 1  - timestamp (seconds Unix time)
% 2  - 
% 3  -
% 4  -
% 5  -
% 6  - X Position ECEF (m)
% 7  - Y Position ECEF (m)
% 8  - Z Position ECEF (m)
% 9  - X Position Estimated Error (m, 1-sigma)
% 10 - Y Position Estimated Error (m, 1-sigma)
% 11 - Z Position Estimated Error (m, 1-sigma)
% 12 -
% 13 -
% 14 - X Velocity ECEF (m/s)
% 15 - Y Velocity ECEF (m/s)
% 16 - Z Velocity ECEF (m/s)
% 17 - X Velocity Estimated Error ECEF (m/s, 1-sigma)
% 18 - Y Velocity Estimated Error ECEF (m/s, 1-sigma)
% 19 - Z Velocity Estimated Error ECEF (m/s, 1-sigma)
%
% imudata:
% 1  - timestamp (seconds Unix time)
% 2  - 
% 3  -
% 4  -
% 5  - X Acceleration (m/s/s)
% 6  - Y Acceleration (m/s/s)
% 7  - Z Acceleration (m/s/s)
% 8  - p rate (rad/sec)
% 9  - q rate (rad/sec)
% 10 - r rate (rad/sec)
% 11 - Roll Inclination (rad)
% 12 - Pitch Inclination (rad)
%
% barodata
% 1 - timestamp (seconds unix time)
% 2 - barometric pressure altitude (referenced to 1013.25) (m)
%
% magdata
% 1 - timestamp (seconds unix time)
% 2 - X magnetic field (gauss) 
% 3 - Y magnetic field (gauss)
% 4 - Z magnetic field (gauss)
%
% OUTPUTS
% =======
% results: data structure containing fused data
%
% $Id: GARD_GPSINS_Loose_EKF.m 4095 2010-10-30 01:38:06Z greerd $
%

% get the time base for each sensor
time_imu = imudata(:,1);
time_gps = gpsdata(:,1);
time_baro = barodata(:,1);
time_mag = magdata(:,1);

% expected/ideal data rates of sensors
if ~exist('IMU_RATE','var')
    disp('Setting IMU to default rate of 200Hz');
    IMU_RATE = 200;
end

if ~exist('MAG_RATE','var')
    disp('Setting Magnetometer to default rate of 20Hz');
    MAG_RATE = 20;
end

if ~exist('GPS_RATE')
    disp('Setting GPS to default rate of 10Hz');
    GPS_RATE = 10;
end

if ~exist('BARO_RATE')
    disp('Setting BARO to default rate of 120Hz');
    BARO_RATE = 120;
end



NumberStates = 15;

% average the first 10 seconds of roll and pitch gravity vector
CoarseRoll = mean(imudata(1:10*IMU_RATE,11));
CoarsePitch = mean(imudata(1:10*IMU_RATE,12));

% find the initial dcm
C_BN = GARDSim_DCMfromEuler(CoarseRoll,CoarsePitch,0);



if size(magdata) ~= [1 1]
    disp('Calculating Heading Data');

    MagVariation = 11*pi/180;
    if(~exist('MagHeading'))
        MagHeading = zeros(1,length(magdata));
        for(index=1:1000)
           MagHeading(index) =  (atan2(magdata(index,3),magdata(index,2)));
           if MagHeading(index) > pi
               MagHeading(index) = MagHeading(index) - (2*pi);
           end
           if MagHeading(index) < -pi
               MagHeading(index) = MagHeading(index) + (2*pi);
           end
        end
    end

end

InitialAttitude(1) = CoarseRoll;
InitialAttitude(2) = CoarsePitch;
InitialAttitude(3) = 100*pi/180;% mean(MagHeading(1:10*MAG_RATE));



a = zeros(length(imudata),1);
b = zeros(length(imudata),1);
c = zeros(length(imudata),1);
d = zeros(length(imudata),1);

a(1) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
b(1) = sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) - cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
c(1) = cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
d(1) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2);




C_BN = GARDSim_DCMfromEuler(InitialAttitude(1),InitialAttitude(2),InitialAttitude(3));

% initialise position and attitude vector
Pos_LLH = zeros(length(imudata),3);
Attitude = zeros(length(imudata),3);

phi_q = zeros(length(imudata),1);
theta_q = zeros(length(imudata),1);
psi_q = zeros(length(imudata),1);

phi_q(1) = InitialAttitude(1);
theta_q(1) = InitialAttitude(2);
psi_q(1) = InitialAttitude(3);

disp('Calibrating initial gyro bias');

% calibrate initial gyro bias
GyroBias(1) = mean(imudata(1:30*IMU_RATE,8));
GyroBias(2) = mean(imudata(1:30*IMU_RATE,9));
GyroBias(3) = mean(imudata(1:30*IMU_RATE,10));


AccelBias = zeros(3,1);

%imudata(:,8) = imudata(:,8) - GyroBias(1);
%imudata(:,9) = imudata(:,9) - GyroBias(2);
%imudata(:,10) = imudata(:,10) - GyroBias(3);



dt = 1/IMU_RATE;
gps_dt = 1 / GPS_RATE;
ins_dt = 1/IMU_RATE;

% x_gyro_beta = 1/500;
% y_gyro_beta = 1/500;
% z_gyro_beta = 1/500;
% x_accel_beta = 1/300;
% y_accel_beta = 1/300;
% z_accel_beta = 1/300;
% 
% 
% x_gyro_Q = (1.0e-2)^2;
% y_gyro_Q = (1.0e-2)^2;
% z_gyro_Q = (1.0e-2)^2;
% x_accel_Q = (1.0e-2)^2;
% y_accel_Q = (1.0e-2)^2;
% z_accel_Q = (1.0e-2)^2;

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

% rotate the pos and velocity measurements and errors to NED frame



disp('Calculating GPS Positions');
if(~exist('GPSPos_E','var'))
    
    GPSPos_E = zeros(length(gpsdata),3);
    GPSPos_EE = zeros(length(gpsdata),3);
    GPSVel_NED = zeros(length(gpsdata),3);
    GPSVel_EE = zeros(length(gpsdata),3);
    
    for(i = 1:length(gpsdata))
        [GPSPos_E(i,1), GPSPos_E(i,2), GPSPos_E(i,3)] = ECEF2LLH(gpsdata(i,6:8));
        C_EN = T_ECEF2NED(GPSPos_E(i,1),GPSPos_E(i,2));

        GPSPos_EE(i,:) = abs(C_EN * gpsdata(i,9:11)');
        GPSVel_NED(i,:) = C_EN * gpsdata(i,14:16)';
        GPSVel_EE(i,:) = abs(C_EN * gpsdata(i,17:19)');
        
        if mod(i,1000) == 0
            disp(sprintf('Completed conversion %d of %d',i,length(gpsdata)));
        end
        
    end
    
    GPSAcc_NED(2:length(gpsdata),1) = diff(GPSVel_NED(:,1));
    GPSAcc_NED(2:length(gpsdata),2) = diff(GPSVel_NED(:,2));
    GPSAcc_NED(2:length(gpsdata),3) = diff(GPSVel_NED(:,3));
    
    disp('Calculating GPS Attitudes');
    pseudo_roll = zeros(length(gpsdata),1);
    gps_flight_path_angle = zeros(length(gpsdata),1);
    gps_track = zeros(length(gpsdata),1);
    
    for i=1:length(gpsdata)
    [pseudo_roll(i),gps_flight_path_angle(i),gps_track(i)] = GARD_PseudoAttitude(GPSVel_NED(i,:),GPSAcc_NED(i,:));
    end
end




InitialPosition = [GPSPos_E(1,1),GPSPos_E(1,2), GPSPos_E(1,3)];

LOCAL_GRAVITY = -GravityModel(InitialPosition);


% calculate baro height


if size(barodata) ~= [1 1]
disp('Calculating Baro Height');
% QNH = barodata(1,2) + InitialPosition(3)/(29.7/3.28);
% if(~exist('BaroHeight'))
%     BaroHeight = zeros(1,length(barodata));
%     for(index = 1:length(barodata))
%         BaroHeight(index) = (QNH - barodata(index,2)) * 30/3.28;
%     end
% end
BaroOffset = InitialPosition(3) - barodata(1,2);
if(~exist('BaroHeight'))
    BaroHeight = zeros(1,length(barodata));
    for(index = 1:length(barodata))
        BaroHeight(index) = barodata(index,2)+BaroOffset;
    end
end

end

disp('Initialising Filter variables');

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));

R = sqrt(RM * RP);

% initialise the state vector
x_hat_plus = zeros(15,1);
x_hat_plus(10) = GyroBias(1);
x_hat_plus(11) = GyroBias(2);
x_hat_plus(12) = GyroBias(3);


% initialise state covariance
P_hat_plus = zeros(15,15);
P_hat_plus(1,1) = 0.000001;
P_hat_plus(2,2) = 0.000001;
P_hat_plus(3,3) = 10;
P_hat_plus(4,4) = 1;
P_hat_plus(5,5) = 1;
P_hat_plus(6,6) = 1;
P_hat_plus(7,7) = 0.01;
P_hat_plus(8,8) = 0.01;
P_hat_plus(9,9) = 0.01;
P_hat_plus(10,10) = 0.001;
P_hat_plus(11,11) = 0.001;
P_hat_plus(12,12) = 0.001;
P_hat_plus(13,13) = 0.001;
P_hat_plus(14,14) = 0.001;
P_hat_plus(15,15) = 0.001;


saveindex = 1;
x_hat_save = zeros(length(gpsdata),15);
P_hat_save = zeros(length(gpsdata),15);
%Vel_NED_corrected = zeros(length(gpsdata),3);
%Pos_NED_corrected = zeros(length(gpsdata),3);
Acc_NED = zeros(length(imudata),3);
Vel_NED = zeros(length(imudata),3);
insupdatetime = zeros(length(imudata),1);

pl = zeros(length(imudata),1);

Vel_UVW = zeros(length(imudata),3);

disp('Calculating gravity vector');
accelerometer_filter = ones(1,2000)/2000;
gvec_roll = filter(accelerometer_filter,1,imudata(:,11));
gvec_pitch = filter(accelerometer_filter,1,imudata(:,12));

% Damien - 15th January 2008 - changed because longitude is the first
% argument to this function
Tecef2ned= T_ECEF2NED(GPSPos_E(1,1),GPSPos_E(1,2));


OMEGA_e = 7.292115e-5; % eearth rotation rate (inertial frame)

updatetime = zeros(length(gpsdata),4);
Epoch_lo = 0;
disp('Beginning navigation...');
for index=1:length(imudata)
    
    % TODO perform attitude update
    
    % get sensor measurements
    omega_x = imudata(index,8) - GyroBias(1);
    omega_y = imudata(index,9)  - GyroBias(2);
    omega_z = imudata(index,10)  - GyroBias(3);
    
%     % form the skew symmetric form
    Omega_b_vec = [omega_x;omega_y;omega_z];
    
    OMEGA_b = [ 0         -omega_z     omega_y; ...
                omega_z    0          -omega_x; ...
               -omega_y    omega_x     0        ];


    insupdatetime(index) = imudata(index,1);

    if(index==1)
        % initialise the DCM
        C_BN = GARD_EulerToDCM(InitialAttitude(1),InitialAttitude(2),InitialAttitude(3));
    else



        %OMEGA_in = [0, -OMEGA_e * sin(-27*pi/180), 0; ...
        %            OMEGA_e * sin(-27*pi/180), 0, -OMEGA_e * cos(-27*pi/180);
        %            0, OMEGA_e * cos(-27*pi/180), 0];

        %C_BN = C_BN + (C_BN * OMEGA_b - OMEGA_in * C_BN)*dt;
        C_BN = GARD_DCMUpdate2(C_BN,Omega_b_vec*ins_dt);

    end


    c32 = C_BN(3,2);
    c33 = C_BN(3,3);
    c31 = C_BN(3,1);
    c21 = C_BN(2,1);
    c11 = C_BN(1,1);
    
    phi_q(index) = atan2(c32,c33);
    theta_q(index) = asin(-c31);
    psi_q(index) = atan2(c21,c11);
   

    if(psi_q(index) < -pi)
        psi_q(index) = psi_q(index) + (2*pi);
    end
    if(psi_q(index) > pi)
        psi_q(index) = psi_q(index) - (2*pi);
    end

    % perform velocity update
    % get the navigation frame accelerations
    %Acc_NED(index,:) = (C_BN * (imudata(index,5:7)' - AccelBias))';
    Acc_NED(index,:) = (C_BN * (imudata(index,5:7)'))';
    
    % correct vertical channel for gravity
    Acc_NED(index,3) = Acc_NED(index,3) + LOCAL_GRAVITY;
    
    % TODO: correct for coriolis effects
             
%         
%   

    OMEGA_e_n = Tecef2ned * [0;0;7.292115e-5];
    
    if(index>1)
         Coriolis = cross(2 * OMEGA_e_n,Vel_NED(index-1,:));
    else
        Coriolis = cross(2 * OMEGA_e_n,[0,0,0]);
    end
    
    Acc_NED(index,:) = Acc_NED(index,:) - Coriolis;
    
    % increment the velocity and position 
    
    if(index == 1)
        Vel_NED(index,:) = Acc_NED(index,:)*dt*0;
        Pos_LLH(index,:) = InitialPosition;
    else
        % trapezoidal velocity update
        Vel_NED(index,:) = Vel_NED(index-1,:) + Acc_NED(index,:)*dt;% + (Acc_NED(index,:)-Acc_NED(index-1,:))/2;
        % position update
        RMh = RM + Pos_LLH(index-1,3);
        RPh = RP + Pos_LLH(index-1,3);
        lat_dot = Vel_NED(index,1) / (RMh);
        long_dot = Vel_NED(index,2) / (cos(Pos_LLH(index-1,1))*(RPh));
        Pos_LLH(index,1) = Pos_LLH(index-1,1) + lat_dot*dt;
        Pos_LLH(index,2) = Pos_LLH(index-1,2) + long_dot*dt;
        Pos_LLH(index,3) = Pos_LLH(index-1,3) + -Vel_NED(index,3)*dt;
        
    end
    
    
   
    % save the GPS time for this data
    if index ~= 1
        GPSTime(index) = GPSTime(index-1) + 1/IMU_RATE;
    else
        GPSTime(index) = 0;
    end
    if(mod(index,IMU_RATE/GPS_RATE) == 0)
    %if(0)    
        % run kalman filter update
        Epoch_lo = Epoch_lo + 1;
        % find the index of gps, baro and mag data corresponding to the current
        % time
        gpsindex_prior = find(gpsdata(:,1)<imudata(index,1),1,'last');
        gpsindex_after = find(gpsdata(:,1)>imudata(index,1),1,'first');
        baroindex = find(barodata(:,1)>imudata(index,1),1,'first')-1;
        magindex = find(magdata(:,1)>imudata(index,1),1,'first')-1;
        
        
        if abs(imudata(index,1)-gpsdata(gpsindex_prior,1)) < 0.02
               gpsindex = gpsindex_prior;
        elseif abs(imudata(index,1)-gpsdata(gpsindex_after,1)) < 0.02
               gpsindex = gpsindex_after;
        else
            gpsindex = gpsindex_after;
                
        end
        
        
        
        if(gpsindex<1) gpsindex = 1; end
        
        if isempty(gpsindex)
            disp('gps index out of bound');
            break;
        end
        
        
        
        % check if there is a gap in the gps data
        gps_delay = imudata(index,1)-gpsdata(gpsindex,1);
        if  abs(gps_delay) > 1.0
            disp(sprintf('warning - gps data delay: %d:%d',gpsindex,gps_delay));
            continue;
        end
        
        if(baroindex<1) baroindex = 1; end
        if(magindex<1) magindex = 1; end

       
        if(gpsdata(gpsindex,4) ~= 0)
            disp(sprintf('Warning No GPS Solution at t=%f',gpsdata(gpsindex,1)));
        else
            
            %gps_dt = gpsdata(gpsindex,1) - gpslastupdate;
            
            
            try 
            updatetime(saveindex,1) = imudata(index,1);
            updatetime(saveindex,2) = gpsdata(gpsindex,1);
            %updatetime(saveindex,3) = magdata(magindex,1);
            %updatetime(saveindex,4) = barodata(baroindex,1);
            catch
               disp('[GARD_GPSINS_Loose_EKF] caught exception');
            end
            
            % save the GPS time for this data
            GPSTime(index) = gpsdata(gpsindex,3);
            
            % correct gps measurements for lever arm
            IMUPos_LLH(gpsindex,:) = GPSPos_E(gpsindex,:)' - (C_BN * GPSLeverArm) .* [1/RMh; 1/RPh*cos(Pos_LLH(index-1,2)); 1];
            IMUVel_NED(gpsindex,:) = GPSVel_NED(gpsindex,:)' - C_BN * cross(Omega_b_vec,GPSLeverArm);
            
            % initialise measurements
            z = zeros(7,1);
            z(1:3) = Pos_LLH(index,:) - IMUPos_LLH(gpsindex,:);
            z(4:6) = Vel_NED(index,:) - IMUVel_NED(gpsindex,:);
           
            if size(magdata) == [1 1]
                HeadingError = 0;
                

            else
                MagHeading(magindex) = GARD_CompassHeading(magdata(magindex,2:4),CompassCalibration(1),CompassCalibration(2),...
                                                            CompassCalibration(3),CompassCalibration(4),phi_q(index),theta_q(index));
                HeadingError = psi_q(index) - MagHeading(magindex);

                if(HeadingError > pi)
                    HeadingError = HeadingError - 2*pi;
                end
                if(HeadingError < -pi)
                    HeadingError = HeadingError + 2*pi;
                end
            end
            z(7) = 0; %-HeadingError;
            
            z(8) = 0; %Pos_LLH(index,3) - BaroHeight(baroindex);
            
            
             % propagate covariance

    F_INS = GARD_GenerateINSFMatrix(Pos_LLH(index,:),Vel_NED(index,:),Acc_NED(index,:)-[0 0 LOCAL_GRAVITY]);

    F_INS(7:9,10:12) = -C_BN; % tilt error to gyro bias
    F_INS(4:6,13:15) = C_BN; % velocity error to accel bias

    F_INS(10,10) = - x_gyro_beta;
    F_INS(11,11) = - y_gyro_beta;
    F_INS(12,12) = - z_gyro_beta;
    F_INS(13,13) = - x_accel_beta;
    F_INS(14,14) = - y_accel_beta;
    F_INS(15,15) = - z_accel_beta;

    PHI_k = expm(F_INS*gps_dt);

    Q = zeros(15,15);

    Q(10,10) = 2*x_gyro_beta*x_gyro_Q;
    Q(11,11) = 2*y_gyro_beta*y_gyro_Q;
    Q(12,12) = 2*z_gyro_beta*z_gyro_Q;
    Q(13,13) = 2*x_accel_beta*x_accel_Q;
    Q(14,14) = 2*y_accel_beta*y_accel_Q;
    Q(15,15) = 2*z_accel_beta*z_accel_Q;

    G = zeros(15,15);
    G(7:9,10:12) = -C_BN; % tilt error to gyro bias
    G(4:6,13:15) = C_BN; % velocity error to accel bias
    G(10:12,10:12) = eye(3,3);
    G(13:15,13:15) = eye(3,3);

    % convert to discrete time form
    Q_k = PHI_k * (G * Q * G') * PHI_k' * gps_dt;
    %Q_k = (G * Q * G') * gps_dt;

    
            % setup state transition matrix

            
            
            H_k = zeros(7,15);
            H_k(1,1) = 1;
            H_k(2,2) = 1;
            H_k(3,3) = 1;
            H_k(4,4) = 0;
            H_k(5,5) = 0;
            H_k(6,6) = 0;
            H_k(7,9) = 0;  %% compass to heading
            H_k(8,3) = 0;  %% baro to alt
            %H_k(9,7) = 1;  %% horizon to roll
            %H_k(10,8) = 1;  %% horizon to pitch
            
            
            HSize = size(H_k,1);
                
            clear ObsMatrix;
            for i=1:NumberStates
               ObsMatrix((i-1)*HSize+1:i*HSize,1:NumberStates) = H_k *  PHI_k^(i-1);
            end

            ObsRank(Epoch_lo) = rank(ObsMatrix);   

            % the gps position error estimates are needed here
            %GPS_Pos_EE = C_EN * gpsdata(gpsindex,9:11)';
            %GPS_Vel_EE = C_EN * gpsdata(gpsindex,17:19)';

            R_k = zeros(7,7);
            R_k(1,1) = (5/RM)^2;
            R_k(2,2) = (5/RP*cos(Pos_LLH(index,1)))^2;
            R_k(3,3) = 7.5^2;
            R_k(4,4) = 1^2;
            R_k(5,5) = 1^2;
            R_k(6,6) = 1^2;

            % mag heading
            R_k(7,7) = (10 * pi/180)^2;
            
            
            % baralt
            if(imudata(index,1) < 200)
                R_k(8,8) = 20^2;
            else
                R_k(8,8) = 20^2;
            end
            % convert to discrete time form
            %R_k = R_k / gps_dt;

            % prediction
            x_hat_minus = PHI_k * x_hat_plus;
            P_hat_minus = (PHI_k * P_hat_plus * PHI_k') + Q_k;
    
            % update
            z_hat = H_k * x_hat_minus;
            V_k = H_k * P_hat_minus * H_k' + R_k;

            % matrix inverse via svd
            [U S V] = svd(V_k);
            K_k = P_hat_minus * H_k' * V * diag(1./diag(S)) * U';
            
            %K_k = P_hat_minus * H_k' * inv(V_k);
            

            % get the innovation
            v_k = z - z_hat;
            v_save(saveindex,:) = v_k;

            x_hat_plus = x_hat_minus + K_k * (v_k);
            P_hat_plus = P_hat_minus - K_k * H_k * P_hat_minus;

            % save state vector and P matrix
            x_hat_save(saveindex,:) = x_hat_plus;
            P_hat_save(saveindex,:) = diag(P_hat_plus);
            z_save(saveindex,:) = z;
            

            saveindex = saveindex + 1;
            
            %% todo - check integrity

            %% apply correction
            Vel_NED(index,:) = Vel_NED(index,:) - x_hat_plus(4:6)';
            Pos_LLH(index,:) = Pos_LLH(index,:) - x_hat_plus(1:3)';
            
            % convert tilt error to DCM Update
            del_alpha = x_hat_plus(7);
            del_beta = x_hat_plus(8);
            del_gamma = x_hat_plus(9);

            % save bias correction
            %GyroBias = GyroBias + x_hat_plus(10:12)';
            %AccelBias = AccelBias + x_hat_plus(13:15);
            GyroBias = x_hat_plus(10:12)';
            AccelBias = x_hat_plus(13:15);
            
            % reset state vector
            x_hat_plus(1:9) = 0;


            % correct atttitude
            del_att_skew = [0         -del_gamma   del_beta; ...
                           del_gamma  0          -del_alpha; ...
                           -del_beta  del_alpha   0];
             
            C_BN = (eye(3,3) + del_att_skew) * C_BN;
            
            C_BN = GARD_OrthogonaliseDCM(C_BN);
             
            
            c32 = C_BN(3,2);
            c33 = C_BN(3,3);
            c31 = C_BN(3,1);
            c21 = C_BN(2,1);
            c11 = C_BN(1,1);

            phi_q(index) = atan2(c32,c33);
            theta_q(index) = asin(-c31);
            psi_q(index) = atan2(c21,c11);


            if(psi_q(index) < -pi)
                psi_q(index) = psi_q(index) + (2*pi);
            end
            if(psi_q(index) > pi)
                psi_q(index) = psi_q(index) - (2*pi);
            end
            
            
            GyroBias_save(Epoch_lo,:) = GyroBias;
            AccelBias_save(Epoch_lo,:) = AccelBias;
            
            LOCAL_GRAVITY = -GravityModel(Pos_LLH(index,:));
        end
        
    end
    
    
    % update body velocity
    Vel_UVW(index,:) = (C_BN' * Vel_NED(index,:)')';
    
    if(mod(index,1000) == 0)
       disp(sprintf('Completed epoch %d',index));        
    end
end

% create Pos NED from Pos_LLH
disp('Calculating NED Positions');
Pos_NED = zeros(size(Pos_LLH));
for(index=1:length(Pos_LLH))
   Pos_NED(index,1) = (Pos_LLH(index,1) - InitialPosition(1)) * R; 
   Pos_NED(index,2) = (Pos_LLH(index,2) - InitialPosition(2)) * R * cos(InitialPosition(1));
   Pos_NED(index,3) = -(Pos_LLH(index,3) - InitialPosition(3));
   if(Pos_LLH(index,1) == 0)
        Pos_NED(index,1) = 0;
    end
    if(Pos_LLH(index,2) == 0)
        Pos_NED(index,2) = 0;
    end
    if(Pos_LLH(index,3) == 0)
        Pos_NED(index,3) = 0;
    end
end
       
results.GPSTime = GPSTime;
results.GPSPos_E = GPSPos_E;
results.GPSPos_EE = GPSPos_EE;
results.GPSVel_NED = GPSVel_NED;
results.GPSVel_EE = GPSVel_EE;
results.Pos_LLH = Pos_LLH;
results.Pos_NED = Pos_NED;
results.Acc_NED = Acc_NED;
results.Vel_NED = Vel_NED;
results.Vel_UVW = Vel_UVW;
results.phi = phi_q;
results.theta = theta_q;
results.psi = psi_q;
results.x_hat = x_hat_save;
results.P_hat = P_hat_save;
results.z = z_save;
results.C_BN = C_BN;
results.ObsRank = ObsRank;
results.GyroBias = GyroBias_save;
results.AccelBias = AccelBias_save;
results.PseudoAtt =  [pseudo_roll,gps_flight_path_angle,gps_track];
results.updatetime = updatetime;
results.insupdatetime = insupdatetime;
results.v = v_save;
results.IMUPos_LLH = IMUPos_LLH;
results.IMUVel_NED = IMUVel_NED;

%results.MagHeading = MagHeading;

