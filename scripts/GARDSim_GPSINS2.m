% GARDSim_GPSINS2.m
%
%
% Implements an integrated GPS-Inertial position solution
%
% Written by Duncan Greer
%
% $Id: GARDSim_GPSINS2.m 1879 2008-07-15 05:20:21Z n2523710 $
%
%
% The GPS data is at 1 second epochs.  The INS data is calculated at 50Hz,
% or 0.02 second epochs.  The EKF will run at 1Hz, calculating the INS
% system errors
%
% The linearisation takes place about the refernce trajectory provided by
% kalman filter output.  This implementation is a feedback implementation
% where the estimated errors are fed back to the INS output.  
%
% ==== STATES ====
%
% x1 = latitude position error (rad)
% x2 = longitude position error (rad)
% x3 = height position error (m)
% x4 = north velocity error (m/s)
% x5 = east velocity error (m/s)
% x6 = down velicty error (m/s)
% x7 = roll error (rad)
% x8 = pitch error (rad)
% x9 = yaw error (rad)
% x10 = x acc bias error (m/s/s)
% x11 = y acc bias error (m/s/s)
% x12 = z acc bias error (m/s/s)
% x13 = x gyro bias error (m/s/s)
% x14 = y gyro bias error (m/s/s)
% x15 = z gyro bias error (m/s/s)
%
% ==== MEASURMENTS ====
% z1 = latitude error (INS - GPS)
% z2 = longitude error (INS - GPS)
% z3 = height error (INS - GPS)
% z4 = north velocity error (INS - GPS)
% z5 = east velocity error (INS - GPS)
% z6 = down velocity error (INS - GPS)
% z7 = roll error (attitude - g_vector)
% z8 = pitch error (attitude - g_vector)
% z9 = yaw error (attitude - heading measurement)
%
% There are two processes going on - a high speed loop evaluating the INS
% output, and a low-speed loop evaluating the INS errors in a GPS-INS EKF.
% 
%

NumberStates = 15;
NumberMeasurements = 9;

% load truth data
load 'data/long_flight/pos_truth_llh.mat';
load 'data/long_flight/vel_truth.mat';
load 'data/long_flight/att_truth.mat';

% load GPS solution data
load 'data/long_flight/gps_output.mat';

% load IMU measurements data
UseNoisy = 0;

if UseNoisy == 1
    load('data/long_flight/sensors_noisy.mat');
    sensors = sensors_noisy;
    clear 'sensors_noisy';
else
    load('data/long_flight/sensors_clean.mat');
    sensors = sensors_clean;
    clear 'sensors_clean';
end

% setup constants

WGS84Constants;
d2r = pi/180;
r2d = 180/pi;


% position
InitialPosition = [-27.4*pi/180 153.2*pi/180 4500/3.28]';  % somewhere near brisbane in LLH, rads and meters
% velocity
InitialVelocity = [0 -60 0]';  % NED velocities in m/s
% attitude
InitialAttitude = [0 0 0]';  % roll, pitch, yaw in radians

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(2));
RP = PrimeRadius(InitialPosition(2));

% gravity model - Kayton eqn 2.6
g = -9.79;  % m/s/s
%g = -0.01 * (978.049 * (1 + .00529 * sin(InitialPosition(1))^2));


% setup error statistics
sigma_lat = 50 / RM;
sigma_long = 50 / RP*cos(InitialPosition(1));
sigma_h = 30;
sigma_v = 5;

% fault detection
Pfa = 1/15000;

% detection thresholds
chi = chi2inv(1-Pfa,1);

Thresh_lat = chi * sigma_lat;
Thresh_long = chi * sigma_long;
Thresh_h = chi * sigma_h;
Thresh_v = chi * sigma_v;

SD = sqrt(Thresh_lat^2 + Thresh_long ^2 + Thresh_h^2 + Thresh_v^2);


%NumberEpochsGPS = length(gps_out_llh);
NumberEpochsGPS = 1000;
gps_dt = 1;
ins_dt = 0.01;

TimeINS = [0:ins_dt:NumberEpochsGPS-1];
TimeGPS = [0:gps_dt:NumberEpochsGPS-1];

%dt = gps_dt;

NumberEpochsINS = length(TimeINS);


% first run - setup initial conditions
x_hat_in = zeros(NumberStates,1);
P_in = zeros(NumberStates,NumberStates);
x_hat_save = zeros(NumberStates,NumberEpochsGPS);
P_save = zeros(NumberStates,NumberEpochsGPS);

disp(sprintf('[GARDSim_GPSINS2] Beginning navigation loop for %d Epochs',NumberEpochsINS));

% initialise epochs
Epoch_lo = 0;
Epoch_hi = 0;

% initialise storage vars to svae time
Lat_error = zeros(1,NumberEpochsINS);
Long_error = zeros(1,NumberEpochsINS);
Height_error = zeros(1,NumberEpochsINS);
P_xn = zeros(1,NumberEpochsINS);
P_yn = zeros(1,NumberEpochsINS);
P_zn = zeros(1,NumberEpochsINS);
V_xn = zeros(1,NumberEpochsINS);
V_yn = zeros(1,NumberEpochsINS);
V_zn = zeros(1,NumberEpochsINS);

a = zeros(1,NumberEpochsINS);
b = zeros(1,NumberEpochsINS);
c = zeros(1,NumberEpochsINS);
d = zeros(1,NumberEpochsINS);

phi_q = zeros(1,NumberEpochsINS);
theta_q = zeros(1,NumberEpochsINS);
psi_q = zeros(1,NumberEpochsINS);

phi_uc = zeros(1,NumberEpochsINS);
theta_uc = zeros(1,NumberEpochsINS);
psi_uc = zeros(1,NumberEpochsINS);

% high-speed loop - INS evaluation
for Epoch_hi = 1:NumberEpochsINS
   if Epoch_hi == 1
       % first time step
       a(Epoch_hi) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
       b(Epoch_hi) = sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) - cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
       c(Epoch_hi) = cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
       d(Epoch_hi) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2);
       
       % initialise geodetic position
       Pos_LLH(Epoch_hi,:) = InitialPosition;       
        Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch_hi,2),Pos_LLH(Epoch_hi,1));
        Tned2ecef = Tecef2ned';
        dt = 0.01;
        
        P_xn(Epoch_hi) = InitialPosition(1);
        P_yn(Epoch_hi) = InitialPosition(2);
        P_zn(Epoch_hi) = InitialPosition(3);
        V_xn(Epoch_hi) = InitialVelocity(1);
        V_yn(Epoch_hi) = InitialVelocity(2);
        V_zn(Epoch_hi) = InitialVelocity(3);
        V_n = [V_xn(Epoch_hi);V_yn(Epoch_hi);V_zn(Epoch_hi)];
   else
        ins_dt = TimeINS(Epoch_hi) - TimeINS(Epoch_hi-1);
        
        % get sensor measurements
        omega_x = sensors(5,Epoch_hi);% + randn(1) * 0.0055;
        omega_y = sensors(6,Epoch_hi);% + randn(1) * 0.0055;
        omega_z = sensors(7,Epoch_hi);% + randn(1) * 0.0055;

        
        a_dot = -0.5 * (b(Epoch_hi-1) * omega_x + c(Epoch_hi-1) * omega_y + d(Epoch_hi-1) * omega_z);
        b_dot = 0.5 * (a(Epoch_hi-1) * omega_x - d(Epoch_hi-1) * omega_y + c(Epoch_hi-1) * omega_z);
        c_dot = 0.5 * (d(Epoch_hi-1) * omega_x + a(Epoch_hi-1) * omega_y - b(Epoch_hi-1) * omega_z);
        d_dot = -0.5 * (c(Epoch_hi-1) * omega_x - b(Epoch_hi-1) * omega_y - a(Epoch_hi-1) * omega_z);
        
        a(Epoch_hi) = a(Epoch_hi-1) + a_dot * dt;
        b(Epoch_hi) = b(Epoch_hi-1) + b_dot * dt;
        c(Epoch_hi) = c(Epoch_hi-1) + c_dot * dt;
        d(Epoch_hi) = d(Epoch_hi-1) + d_dot * dt;
        
        Tecef2ned= T_ECEF2NED(P_xn(Epoch_hi-1),P_yn(Epoch_hi-1));
        Tned2ecef = Tecef2ned';
        V_n = [V_xn(Epoch_hi-1);V_yn(Epoch_hi-1);V_zn(Epoch_hi-1)];
    
    
    A_xb = sensors(2,Epoch_hi);% + randn(1) * 0.04; 
    A_yb = sensors(3,Epoch_hi);% + randn(1) * 0.04;
    A_zb = sensors(4,Epoch_hi);% + randn(1) * 0.04;

    % find gravity vector for attitude estimate
    gvec_phi(Epoch_hi) =  atan2(-A_yb,sqrt(A_xb^2 + A_zb^2));% roll
    gvec_theta(Epoch_hi) = atan2(A_xb,-A_zb); % pitch
    
    c11 = (a(Epoch_hi)^2 + b(Epoch_hi)^2 - c(Epoch_hi)^2 - d(Epoch_hi)^2);
    c12 = 2 * (b(Epoch_hi)*c(Epoch_hi) - a(Epoch_hi)*d(Epoch_hi));
    c13 = 2 * (b(Epoch_hi)*d(Epoch_hi) + a(Epoch_hi)*c(Epoch_hi));
    c21 = 2 * (b(Epoch_hi)*c(Epoch_hi) + a(Epoch_hi)*d(Epoch_hi));
    c22 = (a(Epoch_hi)^2 - b(Epoch_hi)^2 + c(Epoch_hi)^2 - d(Epoch_hi)^2);
    c23 = 2 * (c(Epoch_hi)*d(Epoch_hi) - a(Epoch_hi)*b(Epoch_hi));
    c31 = 2 * (b(Epoch_hi)*d(Epoch_hi) - a(Epoch_hi)*c(Epoch_hi));
    c32 = 2 * (c(Epoch_hi)*d(Epoch_hi) + a(Epoch_hi)*b(Epoch_hi));
    c33 = (a(Epoch_hi)^2 - b(Epoch_hi)^2 - c(Epoch_hi)^2 + d(Epoch_hi)^2);

    C_bn = [c11, c12, c13; c21, c22, c23; c31, c32, c33];

    phi_q(Epoch_hi) = atan2(c32,c33);
    theta_q(Epoch_hi) = asin(-c31);
    psi_q(Epoch_hi) = atan2(c21,c11);

    % save the uncorrected attitude
    phi_uc(Epoch_hi) = phi_q(Epoch_hi);
    theta_uc(Epoch_hi) = theta_q(Epoch_hi);
    psi_uc(Epoch_hi) = psi(Epoch_hi);
    
    % rotate accelerometer measurements to nav frame using A_n = C_bn *
    % A_b

    A_n = C_bn * [A_xb;A_yb;A_zb];
    A_xn = A_n(1);
    A_yn = A_n(2);
    A_zn = A_n(3) - g;

    % perform coriolis correction


    OMEGA_e_n = Tecef2ned * [0;0;(15/3600)*d2r];
    
    Coriolis = cross(2 * OMEGA_e_n,V_n);
    A_xn = A_xn - Coriolis(1);
    A_yn = A_yn - Coriolis(2);
    A_zn = A_zn - Coriolis(3);
 
    % RECTANGULAR INTEGRATION POSITION UPDATE       
    % propogate velocity estimate
    V_xn(Epoch_hi) = V_xn(Epoch_hi-1) + A_xn * dt;
    V_yn(Epoch_hi) = V_yn(Epoch_hi-1) + A_yn * dt;
    V_zn(Epoch_hi) = V_zn(Epoch_hi-1) + A_zn * dt;

    % determine latitude and longitude rates

    lat_dot = V_xn(Epoch_hi) / (RM + P_zn(Epoch_hi-1));
    long_dot = V_yn(Epoch_hi) / (cos(P_xn(Epoch_hi-1))*(RP + P_zn(Epoch_hi-1)));
    
    % propogate position estimate
    P_xn(Epoch_hi) = P_xn(Epoch_hi-1) + lat_dot * dt;
    P_yn(Epoch_hi) = P_yn(Epoch_hi-1) + long_dot * dt;
    P_zn(Epoch_hi) = P_zn(Epoch_hi-1) + V_zn(Epoch_hi) * dt;
    
    

    
   end
    
    % low-speed loop - kalman filter
    if(mod(Epoch_hi-1,100) == 0 && Epoch_hi ~= 1)
            
            Epoch_lo = Epoch_lo + 1;

            
            if 1%(Epoch_hi < 2000*50 || Epoch_hi > 2800*50)
            
                % get gps velocity in NED
                gps_vel_ecef = gps_out_ecef(Epoch_lo,5:7)';
                Long = gps_out_llh(Epoch_lo,2) * d2r;
                Lat = gps_out_llh(Epoch_lo,1) * d2r;
                Tecef2ned = T_ECEF2NED(Long,Lat);
                gps_vel_ned(Epoch_lo,:) = (Tecef2ned * gps_vel_ecef)';

                % setup measurements matrix
                z = zeros(NumberMeasurements,1);
                z(1,1) = (P_xn(Epoch_hi) - gps_out_llh(Epoch_lo,1) * d2r);
                z(2,1) = (P_yn(Epoch_hi) - gps_out_llh(Epoch_lo,2) * d2r);
                z(3,1) = (P_zn(Epoch_hi) - gps_out_llh(Epoch_lo,3));
                z(4,1) = (V_xn(Epoch_hi) - gps_vel_ned(Epoch_lo,1));
                z(5,1) = (V_yn(Epoch_hi) - gps_vel_ned(Epoch_lo,2));
                z(6,1) = (V_zn(Epoch_hi) - gps_vel_ned(Epoch_lo,3));
                
                % x1 = latitude position error (rad)
                % x2 = longitude position error (rad)
                % x3 = height position error (m)
                % x4 = north velocity error (m/s)
                % x5 = east velocity error (m/s)
                % x6 = down velicty error (m/s)
                % x7 = roll error (rad)
                % x8 = pitch error (rad)
                % x9 = yaw error (rad)

                % setup phi matrix
                %phi = eye(NumberStates,NumberStates);
                phi = zeros(NumberStates,NumberStates);
                phi(1,4) = gps_dt / (RM + P_zn(Epoch_hi));
                phi(2,5) = gps_dt / (cos(P_xn(Epoch_hi))*(RP + P_zn(Epoch_hi)));
                phi(3,6) = -gps_dt;
                phi(4,8) = -g * dt;
                phi(5,7) = -g * dt;
                phi(7,9) = omega_y * dt;
                phi(7,5) = dt / RM;
                phi(8,9) = omega_x * dt;
                phi(8,4) = dt / RM;
                

                % setup H matrix
                H = zeros(NumberMeasurements,NumberStates);
                H(1:6,1:6) = eye(6,6);

                % setup Q and R matrices
                R = eye(NumberMeasurements,NumberMeasurements);
                
                R(1,1) = 500 * sigma_lat^2;
                R(2,2) = 500 * sigma_long^2;
                R(3,3) = sigma_h^2;
                R(4,4) = 2*sigma_v^2; 
                R(5,5) = 2*sigma_v^2;
                R(6,6) = 2*sigma_v^2;
                
                Q = eye(NumberStates,NumberStates) * gps_dt;
                dt_3 = (gps_dt ^ 3) / 3;
                dt_2 = (gps_dt ^ 2) / 2;
                Q(1,1) = dt_3/RM;
                Q(2,2) = dt_3/(RP*cos(P_xn(Epoch_hi)));
                Q(3,3) = dt_3;
                Q(4,4) = 5;
                Q(5,5) = 5;
                Q(6,6) = 2;
                Q(7,7) = dt/5;
                Q(8,8) = dt/5;
                Q(9,9) = dt/5;
                
                %Q = Q/25;




                % evaluate KF
                [x_hat_out, P_out, v_out, s2_out] = GARD_EvaluateKF(dt, x_hat_in, P_in, phi, H, z, Q, R);

                % find innovation vector
                v_save(:,Epoch_lo) = v_out;
                s2_save(Epoch_lo) = s2_out;

                % save x_hat_out
                x_hat_save(:,Epoch_lo) = x_hat_out;
                P_save(:,Epoch_lo) = diag(P_out);
                % update x_hat_in
                x_hat_in = x_hat_out;
                P_in = P_out;

                %if (Epoch_lo < 1000 || Epoch_lo > 1060)
                if(1)                        
                    P_xn(Epoch_hi) = P_xn(Epoch_hi) - x_hat_save(1,Epoch_lo);
                    P_yn(Epoch_hi) = P_yn(Epoch_hi) - x_hat_save(2,Epoch_lo);
                    P_zn(Epoch_hi) = P_zn(Epoch_hi) - x_hat_save(3,Epoch_lo);
                    V_xn(Epoch_hi) = V_xn(Epoch_hi) - x_hat_save(4,Epoch_lo);
                    V_yn(Epoch_hi) = V_yn(Epoch_hi) - x_hat_save(5,Epoch_lo);
                    V_zn(Epoch_hi) = V_zn(Epoch_hi) - x_hat_save(6,Epoch_lo);
%                     phi_q(Epoch_hi-49:Epoch_hi) = phi_q(Epoch_hi-49:Epoch_hi) - x_hat_save(7,Epoch_lo);
%                     theta_q(Epoch_hi-49:Epoch_hi) = theta_q(Epoch_hi-49:Epoch_hi) - x_hat_save(8,Epoch_lo);
%                     psi_q(Epoch_hi-49:Epoch_hi) = psi_q(Epoch_hi-49:Epoch_hi) - x_hat_save(9,Epoch_lo);
                    
%                 elseif (s2_out < 2.5e5)
%                     
%                     P_xn(Epoch_hi) = P_xn(Epoch_hi) - x_hat_save(1,Epoch_lo);
%                     P_yn(Epoch_hi) = P_yn(Epoch_hi) - x_hat_save(2,Epoch_lo);
%                     P_zn(Epoch_hi) = P_zn(Epoch_hi) - x_hat_save(3,Epoch_lo);
%                     V_xn(Epoch_hi) = V_xn(Epoch_hi) - x_hat_save(4,Epoch_lo);
%                     V_yn(Epoch_hi) = V_yn(Epoch_hi) - x_hat_save(5,Epoch_lo);
%                     V_zn(Epoch_hi) = V_zn(Epoch_hi) - x_hat_save(6,Epoch_lo);
                end
                
            end

    end% end low speed loop
    
%     % determine how long since the last KF loop
%     lapse_dt = (Epoch_hi - (Epoch_lo / ins_dt))*ins_dt;
%     
%     % update for velocity error
%     P_xn(Epoch_hi) = P_xn(Epoch_hi) - (x_hat_save(4,Epoch_lo) * lapse_dt) / (RM + P_zn(Epoch_hi));
%     P_yn(Epoch_hi) = P_yn(Epoch_hi) - (x_hat_save(5,Epoch_lo) * lapse_dt) / (cos(P_xn(Epoch_hi)) * (RP + P_zn(Epoch_hi)));
%     P_zn(Epoch_hi) = P_zn(Epoch_hi) + x_hat_save(6,Epoch_lo) * dt;
    
    
    Lat_error(Epoch_hi) = (P_xn(Epoch_hi) - pos_truth_llh(2,Epoch_hi))*RM;
    Long_error(Epoch_hi) = (P_yn(Epoch_hi) - pos_truth_llh(3,Epoch_hi))*RP*cos(P_xn(Epoch_hi));
    Height_error(Epoch_hi) = (P_zn(Epoch_hi) - pos_truth_llh(4,Epoch_hi));

    if(mod(Epoch_hi,1000) == 0)
        disp(sprintf('[GARDSim_GPSINS2] Completed Epoch %d',Epoch_hi));
    end
    
end % end high-speed loop
%NumberEpochsINS = length(PKF_x_save);



% disp('Applying corrections');
% for Epoch = 1:NumberEpochsINS
%     %if Epoch<50
%     %    gpsins(Epoch,1) = PKF_x_save(Epoch,1);
%     %    gpsins(Epoch,2) = PKF_x_save(Epoch,4);
%     %    gpsins(Epoch,3) = PKF_x_save(Epoch,7);
%     %else
%         % apply corrections - positoin and velocity
%         gpsins(Epoch,1) = PKF_x_save(Epoch,1) - x_hat_save(1,floor(Epoch/50)+1) - (x_hat_save(4,floor(Epoch/50)+1)*(Epoch-50*floor(Epoch/50)+1)*ins_dt / (RM+PKF_x_save(Epoch,7)));
%         gpsins(Epoch,2) = PKF_x_save(Epoch,4) - x_hat_save(2,floor(Epoch/50)+1) - (x_hat_save(5,floor(Epoch/50)+1)*(Epoch-50*floor(Epoch/50)+1)*ins_dt / (cos(PKF_x_save(Epoch,1))*(RP+PKF_x_save(Epoch,7))));
%         gpsins(Epoch,3) = PKF_x_save(Epoch,7) - x_hat_save(3,floor(Epoch/50)+1) + (x_hat_save(6,floor(Epoch/50)+1)*(Epoch-50*floor(Epoch/50)+1)*ins_dt);
%     %end
%     
%     Lat_error(Epoch) = (gpsins(Epoch,1) - pos_truth_llh(2,Epoch))*RM;
%     Long_error(Epoch) = (gpsins(Epoch,2) - pos_truth_llh(3,Epoch))*RP*cos(gpsins(Epoch,1));
%     Height_error(Epoch) = (gpsins(Epoch,3) - pos_truth_llh(4,Epoch));
%     
% end

disp('[GARDSim_GPSINS] Done!');


% calculate attitde error
phi_err = (att_truth(2,1:length(phi_q)) - phi_q);
theta_err = (att_truth(3,1:length(theta_q)) - theta_q);
psi_err = (att_truth(4,1:length(psi_q)) - psi_q);

for Epoch=1:NumberEpochsINS
    if(psi_err(Epoch) > 180)
        psi_err(Epoch) = psi_err(Epoch) - 2*pi;
    end
end


% plot results

% plot velocity
figure();
plot(TimeINS,V_xn,'b')
hold on;
grid on;
plot(TimeGPS(1:Epoch_lo),gps_vel_ned(:,1),'r.');
plot(TimeINS,vel_truth(2,1:NumberEpochsINS),'g');
xlabel('Time (sec)');
ylabel('North Velocity (m/s)');

figure();
plot(TimeINS,V_yn,'b')
hold on;
grid on;
plot(TimeGPS(1:Epoch_lo),gps_vel_ned(:,2),'r.');
plot(TimeINS,vel_truth(3,1:NumberEpochsINS),'g')
xlabel('Time (sec)');
ylabel('East Velocity (m/s)');

figure();
plot(TimeINS,V_zn,'b')
hold on;
grid on;
plot(TimeGPS(1:Epoch_lo),gps_vel_ned(:,3),'r.');
plot(TimeINS,vel_truth(4,1:NumberEpochsINS),'g');
xlabel('Time (sec)');
ylabel('Down Velocity (m/s)');

%latitude
figure()
hold on;
plot(TimeINS,Lat_error,'b');
grid on;
xlabel('Simulation Time (sec)');
ylabel('North-South Latitude Error (m)');
hold off;

% longitude
figure()
hold on;
plot(TimeINS,Long_error,'b');
grid on;
xlabel('Simulation Time (sec)');
ylabel('East-West Longitude Error (m)');
hold off;

% height
figure()
hold on;
plot(TimeINS,Height_error,'b');
grid on;
xlabel('Simulation Time (sec)');
ylabel('Height Error (m)');
hold off;

figure();
grid on;
hold on;
plot(P_yn*r2d,P_xn*r2d,'.');
plot(gps_out_llh(:,2),gps_out_llh(:,1),'r.');
plot(pos_truth_llh(3,:)*r2d,pos_truth_llh(2,:)*r2d,'g');
hold off;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
legend('EKF Estimate','GPS Only','True Trajectory');


% attitude errors
figure();
grid on;
hold on;
plot(att_truth(1,1:NumberEpochsINS),phi_q * r2d,'b');
plot(att_truth(1,1:NumberEpochsINS),phi_uc * r2d,'r');
plot(att_truth(1,1:NumberEpochsINS),att_truth(2,1:NumberEpochsINS) * r2d,'g');
hold off;
xlabel('Time (sec)');
ylabel('Roll Angle (deg)');

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
