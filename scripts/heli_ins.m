%
% CSIRO Helicopter GPS-INS Extended Kalman Filter
% Written by Duncan Greer
%
%

% add the gardsim stuff to the path
path(path,'../../phd/Matlab/gard/gardsim/');
 
disp('Loading data...');
load('Flight3.mat');


% get the time base for each sensor
time_imu = imudata(:,1);
time_gps = gpsdata(:,1);
time_baro = barodata(:,1);
time_mag = magdata(:,1);

% expected/ideal data rates of sensors
IMU_RATE = 200;
MAG_RATE = 20;
GPS_RATE = 10;
BARO_RATE = 120;

LOCAL_GRAVITY = 9.80;

disp('done');


% average the first 10 seconds of roll and pitch gravity vector
CoarseRoll = mean(imudata(1:10*IMU_RATE,11));
CoarsePitch = mean(imudata(1:10*IMU_RATE,12));

% find the initial dcm
C_BN = GARDSim_DCMfromEuler(CoarseRoll,CoarsePitch,0);



disp('Calculating Heading Data');

MagVariation = 11*pi/180;
if(~exist('MagHeading'))
    MagHeading = zeros(1,length(magdata));
    for(index=1:length(magdata))
       MagHeading(index) =  atan2(magdata(index,3),magdata(index,2)) + pi + MagVariation;
       if MagHeading(index) > pi
           MagHeading(index) = MagHeading(index) - (2*pi);
       end
       if MagHeading(index) < -pi
           MagHeading(index) = MagHeading(index) + (2*pi);
       end
    end
end

InitialAttitude(1) = CoarseRoll;
InitialAttitude(2) = CoarsePitch;
InitialAttitude(3) = mean(MagHeading(1:10*MAG_RATE));



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
GyroBias(1) = mean(imudata(1:30*IMU_RATE,8))*IMU_RATE;
GyroBias(2) = mean(imudata(1:30*IMU_RATE,9))*IMU_RATE;
GyroBias(3) = mean(imudata(1:30*IMU_RATE,10))*IMU_RATE;


AccelBias = zeros(3,1);

%imudata(:,8) = imudata(:,8) - GyroBias(1);
%imudata(:,9) = imudata(:,9) - GyroBias(2);
%imudata(:,10) = imudata(:,10) - GyroBias(3);



dt = 1/IMU_RATE;
gps_dt = 1 / GPS_RATE;

x_gyro_beta = 1/300;
y_gyro_beta = 1/300;
z_gyro_beta = 1/300;
x_accel_beta = 1/300;
y_accel_beta = 1/300;
z_accel_beta = 1/300;


x_gyro_Q = (IMU_RATE * 4.5e-5)^2;
y_gyro_Q = (IMU_RATE * 4.0e-5)^2;
z_gyro_Q = (IMU_RATE * 4.0e-5)^2;
x_accel_Q = (IMU_RATE * 5.0e-4)^2;
y_accel_Q = (IMU_RATE * 1.0e-4)^2;
z_accel_Q = (IMU_RATE * 1.0e-3)^2;



% rotate the pos and velocity measurements and errors to NED frame



disp('Calculating GPS Positions');
if(~exist('GPS_Pos_E'))
    for(i = 1:length(gpsdata))
        [GPS_Pos_E(i,1), GPS_Pos_E(i,2), GPS_Pos_E(i,3)] = ECEF2LLH(gpsdata(i,6:8));
        C_EN = T_ECEF2NED(GPS_Pos_E(i,2),GPS_Pos_E(i,1));

        GPSPos_EE(i,:) = abs(C_EN * gpsdata(i,9:11)');
        GPSVel_NED(i,:) = C_EN * gpsdata(i,14:16)';
        GPSVel_EE(i,:) = abs(C_EN * gpsdata(i,17:19)');
    end
end




InitialPosition = [GPS_Pos_E(1,1),GPS_Pos_E(1,2), GPS_Pos_E(1,3)];


% calculate baro height

disp('Calculating Baro Height');
QNH = barodata(1,2) + InitialPosition(3)/(29.7/3.28);
if(~exist('BaroHeight'))
    BaroHeight = zeros(1,length(barodata));
    for(index = 1:length(barodata))
        BaroHeight(index) = (QNH - barodata(index,2)) * 30/3.28;
    end
end

disp('Initialising Filter variables');

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));

R = sqrt(RM * RP);

% initialise the state vector
x_hat_plus = zeros(15,1);
% x_hat_plus(10) = GyroBias(1);
% x_hat_plus(11) = GyroBias(2);
% x_hat_plus(12) = GyroBias(3);


% initialise state covariance
P_hat_plus = zeros(15,15);
P_hat_plus(1,1) = 0.0000001;
P_hat_plus(2,2) = 0.0000001;
P_hat_plus(3,3) = 0.5;
P_hat_plus(4,4) = 0.1;
P_hat_plus(5,5) = 0.1;
P_hat_plus(6,6) = 0.1;
P_hat_plus(7,7) = 0.001;
P_hat_plus(8,8) = 0.001;
P_hat_plus(9,9) = 0.001;
P_hat_plus(10,10) = 0.00001;
P_hat_plus(11,11) = 0.00001;
P_hat_plus(12,12) = 0.00001;
P_hat_plus(13,13) = 0.00001;
P_hat_plus(14,14) = 0.00001;
P_hat_plus(15,15) = 0.00001;


saveindex = 1;
x_hat_save = zeros(length(gpsdata),15);
P_hat_save = zeros(length(gpsdata),15);
%Vel_NED_corrected = zeros(length(gpsdata),3);
%Pos_NED_corrected = zeros(length(gpsdata),3);
Acc_NED = zeros(length(imudata),3);
Vel_NED = zeros(length(imudata),3);


Vel_UVW = zeros(length(imudata),3);

disp('Calculating gravity vector');
accelerometer_filter = ones(1,2000)/2000;
gvec_roll = filter(accelerometer_filter,1,imudata(:,11));
gvec_pitch = filter(accelerometer_filter,1,imudata(:,12));

Tecef2ned= T_ECEF2NED(GPS_Pos_E(1,1),GPS_Pos_E(1,2));


OMEGA_e = 7.292115e-5; % eearth rotation rate (inertial frame)

updatetime = zeros(length(gpsdata),4);

disp('Beginning navigation...');
for(index=1:150000)%182000)
    
    % TODO perform attitude update
    
    % get sensor measurements
    omega_x = imudata(index,8)*IMU_RATE  - GyroBias(1);
    omega_y = imudata(index,9)*IMU_RATE  - GyroBias(2);
    omega_z = imudata(index,10)*IMU_RATE  - GyroBias(3);
    
%     % form the skew symmetric form
    OMEGA_b = [ 0         -omega_z     omega_y; ...
                omega_z    0          -omega_x; ...
               -omega_y    omega_x     0        ];



    if(index==1)


        % initialise the DCM
        C_BN = GARDSim_DCMfromEuler(InitialAttitude(1),InitialAttitude(2),InitialAttitude(3));
    else



        OMEGA_in = [0, -OMEGA_e * sin(-27*pi/180), 0; ...
                    OMEGA_e * sin(-27*pi/180), 0, -OMEGA_e * cos(-27*pi/180);
                    0, OMEGA_e * cos(-27*pi/180), 0];

        C_BN = C_BN + (C_BN * OMEGA_b - OMEGA_in * C_BN)*dt;


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
    Acc_NED(index,:) = (C_BN * (imudata(index,5:7)' - AccelBias/IMU_RATE))';
    %Acc_NED(index,:) = (C_BN * (imudata(index,5:7)'))';
    
    % correct vertical channel for gravity
    Acc_NED(index,3) = Acc_NED(index,3) + LOCAL_GRAVITY/IMU_RATE;
    
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
    
    % increment the velocity and position - note that the accelerations are
    % already in delta-V metres per second
    % in m/s
    
    if(index == 1)
        Vel_NED(index,:) = Acc_NED(index,:);
        Pos_LLH(index,:) = InitialPosition;
    else
        % trapezoidal velocity update
        Vel_NED(index,:) = Vel_NED(index-1,:) + Acc_NED(index,:);% + (Acc_NED(index,:)-Acc_NED(index-1,:))/2;
        % position update
        lat_dot = Vel_NED(index,1) / (RM + Pos_LLH(index-1,3));
        long_dot = Vel_NED(index,2) / (cos(Pos_LLH(index-1,1))*(RP + Pos_LLH(index-1,3)));
        Pos_LLH(index,1) = Pos_LLH(index-1,1) + lat_dot*dt;
        Pos_LLH(index,2) = Pos_LLH(index-1,2) + long_dot*dt;
        Pos_LLH(index,3) = Pos_LLH(index-1,3) + -Vel_NED(index,3)*dt;
        
    end
    
    
    if(mod(index,IMU_RATE/GPS_RATE) == 0)
    %if(0)    
        % run kalman filter update

        % find the index of gps, baro and mag data corresponding to the current
        % time
        gpsindex = find(gpsdata(:,1)>imudata(index,1),1,'first')-1;
        baroindex = find(barodata(:,1)>imudata(index,1),1,'first')-1;
        magindex = find(magdata(:,1)>imudata(index,1),1,'first')-1;
        
        if(gpsindex<1) gpsindex = 1; end
        if(baroindex<1) baroindex = 1; end
        if(magindex<1) magindex = 1; end

       
        if(gpsdata(gpsindex,4) ~= 0)
            disp(sprintf('Warning No GPS Solution at t=%f',gpsdata(gpsindex,1)));
        else
            
            %gps_dt = gpsdata(gpsindex,1) - gpslastupdate;
            
            updatetime(saveindex,1) = imudata(index,1);
            updatetime(saveindex,2) = gpsdata(gpsindex,1);
            updatetime(saveindex,3) = magdata(magindex,1);
            updatetime(saveindex,4) = barodata(baroindex,1);

            % initialise measurements
            z = zeros(7,1);
            z(1:3) = Pos_LLH(index,:) - GPS_Pos_E(gpsindex,:);
            z(4:6) = Vel_NED(index,:) - GPSVel_NED(gpsindex,:);
           
            HeadingError = psi_q(index) - MagHeading(magindex);
            if(HeadingError > pi)
                HeadingError = HeadingError - 2*pi;
            end
            if(HeadingError < -pi)
                HeadingError = HeadingError + 2*pi;
            end
            
            z(7) = -HeadingError;
            
            z(8) = Pos_LLH(index,3) - BaroHeight(baroindex);
            
            
            
            % setup state transition matrix
            
            F = zeros(9,9);
            
            % position error due to velocity error
            F(1,4) = 1/R;
            F(2,5) = 1/(R*cos(-27*pi/180.0));
            F(3,6) = -1;
            
            F(1,3) = -Vel_NED(index,1)/R^2;
            F(2,1) = Vel_NED(index,2)*sin(Pos_LLH(index,1))/(R*cos(Pos_LLH(index,1))^2);
            F(2,3) = -Vel_NED(index,2)/(R*R*cos(Pos_LLH(index,1)));
            F(2,5) = 1/(R*cos(Pos_LLH(index,1)));
            
            F(4,1) = -Vel_NED(index,2)*2*OMEGA_e*cos(Pos_LLH(index,1)) - Vel_NED(index,2)*Vel_NED(index,2)/(R*cos(Pos_LLH(index,1))^2);
            F(4,3) = (Vel_NED(index,2)^2 * tan(Pos_LLH(index,1)) - Vel_NED(index,1)*Vel_NED(index,3))/R^2;
            F(4,4) = Vel_NED(index,3)/R;
            F(4,5) = -2*(OMEGA_e*sin(Pos_LLH(index,1)) + Vel_NED(index,2)*tan(Pos_LLH(index,1))/R);
            F(4,6) = Vel_NED(index,1)/R;
            
            F(5,1) = 2*OMEGA_e*(Vel_NED(index,1)*cos(Pos_LLH(index,1)) - Vel_NED(index,3)*sin(Pos_LLH(index,1))) + (Vel_NED(index,1)*Vel_NED(index,2)/(R*cos(Pos_LLH(index,1))^2));
            F(5,3) = -Vel_NED(index,2)*(Vel_NED(index,1)*tan(Pos_LLH(index,1)) - Vel_NED(index,2)*Vel_NED(index,3))/R^2;
            F(5,4) = 2*OMEGA_e*sin(Pos_LLH(index,1)) + Vel_NED(index,2)*tan(Pos_LLH(index,1))/R;
            F(5,5) = (Vel_NED(index,1)*tan(Pos_LLH(index,1))+Vel_NED(index,3))/R;
            F(5,6) = 2*OMEGA_e * cos(Pos_LLH(index,1)) + Vel_NED(index,2)/R;
            
            
            
            % velocity error due to tilt error
             F(4,8) = -(Acc_NED(index,3)*IMU_RATE - LOCAL_GRAVITY);  %-fd
             F(4,9) = (Acc_NED(index,2)*IMU_RATE); % fe
             F(5,7) = (Acc_NED(index,3)*IMU_RATE - LOCAL_GRAVITY);   %fd
             F(5,9) = -(Acc_NED(index,1)*IMU_RATE); % -fn
             F(6,7) = -(Acc_NED(index,2)*IMU_RATE); % -fe
             F(6,8) = (Acc_NED(index,1)*IMU_RATE); % fn

             
             F(6,1) = 2*OMEGA_e*Vel_NED(index,2)*sin(Pos_LLH(index,1));
             F(6,3) = (Vel_NED(index,1)^2 + Vel_NED(index,2)^2)/R^2;  %% note: gravity correction term missing
             F(6,4) = -2*Vel_NED(index,1)/R;
             F(6,5) = -2*(OMEGA_e * cos(Pos_LLH(index,1)) + Vel_NED(index,2)/R);
             
             F(7,1) = -OMEGA_e * sin(Pos_LLH(index,1));
             F(7,3) = -Vel_NED(index,2)/R^2;
             F(8,3) = Vel_NED(index,1) / R^2;
             F(9,1) = -OMEGA_e * cos(Pos_LLH(index,1)) - Vel_NED(index,2)/(R*cos(Pos_LLH(index,1))^2);
             F(9,3) = Vel_NED(index,2)*tan(Pos_LLH(index,1))/R^2;
            
             % tilt error due to velocity error
            F(7,5) = 1/R;
            F(8,4) = -1/R;           
            F(9,5) = -tan(Pos_LLH(index,1))/R;
            

            
            % tilt error due to inertial rotation of the reference frame
            F(7,8) = -OMEGA_e*sin(Pos_LLH(index,1)) - Vel_NED(index,2)*tan(Pos_LLH(index,1))/R;
            F(7,9) = Vel_NED(index,1)/R;
            F(8,7) = OMEGA_e*sin(Pos_LLH(index,1)) + Vel_NED(index,2)*tan(Pos_LLH(index,1))/R;
            F(9,7) = -Vel_NED(index,1)/R;
            F(9,8) = -OMEGA_e*cos(Pos_LLH(index,1)) - Vel_NED(index,2)/R;
            F(8,9) = OMEGA_e*cos(Pos_LLH(index,1)) + Vel_NED(index,2)/R;
            
            
            PHI_k = eye(15,15);
            PHI_k(1:9,1:9) = eye(9,9) + F * gps_dt;

            
            % augmented state matrix
            PHI_k(7:9,10:12) = -C_BN * gps_dt;
            PHI_k(4:6,13:15) = C_BN * gps_dt;
            
            PHI_k(10,10) = 1 - x_gyro_beta * gps_dt;
            PHI_k(11,11) = 1 - y_gyro_beta * gps_dt;
            PHI_k(12,12) = 1 - z_gyro_beta * gps_dt;
            PHI_k(13,13) = 1 - x_accel_beta * gps_dt;
            PHI_k(14,14) = 1 - y_accel_beta * gps_dt;
            PHI_k(15,15) = 1 - z_accel_beta * gps_dt;

            Q = zeros(15,15);

            Q(10,10) = 2*x_gyro_beta*x_gyro_Q;
            Q(11,11) = 2*y_gyro_beta*y_gyro_Q;
            Q(12,12) = 2*z_gyro_beta*z_gyro_Q;
            Q(13,13) = 2*x_accel_beta*x_accel_Q;
            Q(14,14) = 2*y_accel_beta*y_accel_Q;
            Q(15,15) = 2*z_accel_beta*z_accel_Q;
            
            G = zeros(15,15);
            G(7:9,10:12) = -C_BN;
            G(4:6,13:15) = C_BN;
            G(10:12,10:12) = eye(3,3);
            G(13:15,13:15) = eye(3,3);
            
            % convert to discrete time form
            Q_k = PHI_k * (G * Q * G') * PHI_k' * gps_dt;
            %Q_k = (G * Q * G') * gps_dt;
            
            
            H_k = zeros(7,15);
            H_k(1,1) = 1;
            H_k(2,2) = 1;
            H_k(3,3) = 1;
            H_k(4,4) = 1;
            H_k(5,5) = 1;
            H_k(6,6) = 1;
            H_k(7,9) = 1;
            H_k(8,3) = 1;


            % the gps position error estimates are needed here
            GPS_Pos_EE = C_EN * gpsdata(gpsindex,9:11)';
            GPS_Vel_EE = C_EN * gpsdata(gpsindex,17:19)';

            R_k = zeros(7,7);
            R_k(1,1) = (10/RM)^2;
            R_k(2,2) = (10/RP*cos(Pos_LLH(index,1)))^2;
            R_k(3,3) = 10^2;
            R_k(4,4) = 1^2;
            R_k(5,5) = 1^2;
            R_k(6,6) = 2^2;

            % mag heading
            R_k(7,7) = (10 * pi/180)^2;
            
            
            % baralt
            if(imudata(index,1) < 200)
                R_k(8,8) = 1000^2;
            else
                R_k(8,8) = 2^2;
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
            
            % convert tilt error to a quaternion
            del_alpha = x_hat_plus(7);
            del_beta = x_hat_plus(8);
            del_gamma = x_hat_plus(9);

            % save bias correction
            GyroBias = x_hat_plus(10:12);
            AccelBias = x_hat_plus(13:15);
            
            % reset state vector
            x_hat_plus(1:9) = 0;


            % correct atttitude
            del_att_skew = [0         -del_gamma   del_beta; ...
                           del_gamma  0          -del_alpha; ...
                           -del_beta  del_alpha   0];
             
            C_BN = (eye(3,3) + del_att_skew) * C_BN;
            
            
             
             
        end
        
    end
    
    
    % update body velocity
    Vel_UVW(index,:) = (C_BN' * Vel_NED(index,:)')';
    
    if(mod(index,2000) == 0)
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
        
    


% plot results
plot_results;

