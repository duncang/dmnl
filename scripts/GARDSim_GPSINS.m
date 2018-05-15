% GARDSim_GPSINS.m
%
%
% Implements an integrated GPS-Inertial position solution
%
% Written by Duncan Greer
%
% $Id: GARDSim_GPSINS.m 1879 2008-07-15 05:20:21Z n2523710 $
%
%
% The GPS data is at 1 second epochs.  The INS data is calculated at 50Hz,
% or 0.02 second epochs.  The EKF will run at 1Hz, calculating the INS
% system errors
%
% The linearisation takes place about the refernce trajectory provided by
% the Inertial Navigation System.  This is a direct kalman filter
% implementation, where the system errors are not fed back to the INS.
% Seee GARDSim_GPSINS2 for the feedback implementation.
%
% x1 = latitude position error (rad)
% x2 = longitude position error (rad)
% x3 = height position error (m)
% x4 = north velocity error (m/s)
% x5 = east velocity error (m/s)
% x6 = down velicty error (m/s)

% load truth data
load 'data/long_flight/pos_truth_llh.mat';
load 'data/long_flight/vel_truth.mat';
load 'data/long_flight/att_truth.mat';

% load INS solution data
load 'data/long_flight/ins_output.mat';

% load GPS solution data
load 'data/long_flight/gps_output.mat';

% setup constants

WGS84Constants;
d2r = pi/180;
r2d = 180/pi;


% position
InitialPosition = [-27.38*pi/180 153.12*pi/180 300]';  % somewhere near brisbane in LLH, rads and meters
% velocity
InitialVelocity = [50 0 0]';  % NED velocities in m/s
% attitude
InitialAttitude = [0 0 0]';  % roll, pitch, yaw in radians

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(2));
RP = PrimeRadius(InitialPosition(2));

% gravity model - Kayton eqn 2.6
g = -9.79;  % m/s/s
%g = -0.01 * (978.049 * (1 + .00529 * sin(InitialPosition(1))^2));


%



NumberEpochs = length(gps_out_llh);
gps_dt = 1;
ins_dt = 0.02;

TimeINS = [0:ins_dt:NumberEpochs/ins_dt];
TimeGPS = [0:gps_dt:NumberEpochs/gps_dt];

dt = gps_dt;

% first run - setup initial conditions
x_hat_in = zeros(6,1);
P_in = zeros(6,6);
x_hat_save = zeros(6,NumberEpochs);
P_save = zeros(6,NumberEpochs);

disp(sprintf('[GARDSim_GPSINS] Beginning simulation with %d Epochs',NumberEpochs));
for Epoch=1:NumberEpochs

        
        % get gps velocity in NED
        gps_vel_ecef = gps_out_ecef(Epoch,5:7)';
        Long = gps_out_llh(Epoch,2) * d2r;
        Lat = gps_out_llh(Epoch,1) * d2r;
        Tecef2ned = T_ECEF2NED(Long,Lat);
        gps_vel_ned = Tecef2ned * gps_vel_ecef;
        
        % setup measurements matrix
        z(1,1) = (PKF_x_save(((Epoch-1) * 1/ins_dt) + 1,1) - gps_out_llh(Epoch,1) * d2r);
        z(2,1) = (PKF_x_save(((Epoch-1) * 1/ins_dt) + 1,4) - gps_out_llh(Epoch,2) * d2r);
        z(3,1) = (PKF_x_save(((Epoch-1) * 1/ins_dt) + 1,7) - gps_out_llh(Epoch,3));
        z(4,1) = (PKF_x_save(((Epoch-1) * 1/ins_dt) + 1,2) - gps_vel_ned(1));
        z(5,1) = (PKF_x_save(((Epoch-1) * 1/ins_dt) + 1,5) - gps_vel_ned(2));
        z(6,1) = (PKF_x_save(((Epoch-1) * 1/ins_dt) + 1,8) - gps_vel_ned(3));
        
        % setup phi matrix
        phi = eye(6,6);
        phi(1,4) = dt / (RM + InitialPosition(3));
        phi(2,5) = dt / (cos(InitialPosition(1))*(RP + InitialPosition(3)));
        phi(3,6) = dt;
        
        % setup H matrix
        H = eye(6,6);
        
        % setup Q and R matrices
        R = eye(6,6);
        R(3,3) = 10;
        R(4,4) = 10;
        R(5,5) = 5;
        
        Q = eye(6,6) * dt;
        
        % evaluate KF
        [x_hat_out, P_out] = GARD_EvaluateKF(dt, x_hat_in, P_in, phi, H, z, Q, R);
        
        % save x_hat_out
        x_hat_save(:,Epoch) = x_hat_out;
        P_save(:,Epoch) = diag(P_out);
        
        % update x_hat_in
        x_hat_in = x_hat_out;
        P_in = P_out;

        if(mod(Epoch,100) == 0)
            disp(sprintf('[GARDSim_GPSINS] Completed Epoch %d',Epoch));
        end
    
end

NumberEpochsINS = length(PKF_x_save);

% initialise storage vars to svae time
Lat_error = zeros(1,NumberEpochsINS);
Long_error = zeros(1,NumberEpochsINS);
Height_error = zeros(1,NumberEpochsINS);

disp('Applying corrections');
for Epoch = 1:NumberEpochsINS
    %if Epoch<50
    %    gpsins(Epoch,1) = PKF_x_save(Epoch,1);
    %    gpsins(Epoch,2) = PKF_x_save(Epoch,4);
    %    gpsins(Epoch,3) = PKF_x_save(Epoch,7);
    %else
        % apply corrections - positoin and velocity
        gpsins(Epoch,1) = PKF_x_save(Epoch,1) - x_hat_save(1,floor(Epoch/50)+1) - (x_hat_save(4,floor(Epoch/50)+1)*(Epoch-50*floor(Epoch/50)+1)*ins_dt / (RM+PKF_x_save(Epoch,7)));
        gpsins(Epoch,2) = PKF_x_save(Epoch,4) - x_hat_save(2,floor(Epoch/50)+1) - (x_hat_save(5,floor(Epoch/50)+1)*(Epoch-50*floor(Epoch/50)+1)*ins_dt / (cos(PKF_x_save(Epoch,1))*(RP+PKF_x_save(Epoch,7))));
        gpsins(Epoch,3) = PKF_x_save(Epoch,7) - x_hat_save(3,floor(Epoch/50)+1) + (x_hat_save(6,floor(Epoch/50)+1)*(Epoch-50*floor(Epoch/50)+1)*ins_dt);
    %end
    
    Lat_error(Epoch) = (gpsins(Epoch,1) - pos_truth_llh(2,Epoch))*RM;
    Long_error(Epoch) = (gpsins(Epoch,2) - pos_truth_llh(3,Epoch))*RP*cos(gpsins(Epoch,1));
    Height_error(Epoch) = (gpsins(Epoch,3) - pos_truth_llh(4,Epoch));
    
end

disp('[GARDSim_GPSINS] Done!');

% plot results



% latitude
figure()
hold on;
plot(TimeINS,Lat_error,'g*');
grid on;
xlabel('Simulation Time (sec)');
ylabel('North-South Latitude Error (m)');
hold off;

% longitude
figure()
hold on;
plot(TimeINS,Long_error,'g*');
grid on;
xlabel('Simulation Time (sec)');
ylabel('East-West Longitude Error (m)');
hold off;

% height
figure()
hold on;
plot(TimeINS,Height_error,'g*');
grid on;
xlabel('Simulation Time (sec)');
ylabel('Height Error (m)');
hold off;

% % latitude correction
% figure();
% hold on;
% plot(x_hat_save(1,:) * RM);
% grid on;
% xlabel('Sim Time (sec)');
% ylabel('Latitude Correction (m)');
% hold off;
% 
% % longitude correction
% figure();
% hold on;
% plot(x_hat_save(2,:) * RP * cos(InitialPosition(1)));
% grid on;
% xlabel('Sim Time (sec)');
% ylabel('Longitude Correction (m)');
% hold off;
% 
% % height correction
% figure();
% hold on;
% plot(x_hat_save(3,:));
% grid on;
% xlabel('Sim Time (sec)');
% ylabel('Height Correction (m)');
% hold off;
% 
% % north velocity correction
% figure();
% hold on;
% plot(x_hat_save(4,:));
% grid on;
% xlabel('Sim Time (sec)');
% ylabel('North Velocity Correction (m/s)');
% hold off;
% 
% % east velocity correction
% figure();
% hold on;
% plot(x_hat_save(5,:));
% grid on;
% xlabel('Sim Time (sec)');
% ylabel('East Velocity Correction (m/s)');
% hold off;
% 
% figure();
% hold on;
% plot(x_hat_save(6,:));
% grid on;
% xlabel('Sim Time (sec)');
% ylabel('Down Velocity Correction (m/s)');
% hold off;
