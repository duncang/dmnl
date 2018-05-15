% GARDSim_INS.m
%
% GARDSim Inertial Navigation System Implementation
% Written by Duncan Greer
% 
% This module implments an inertial navigation system.  The first section
% propogates the attitude estimate using an EKF.  
%
% $Id: GARDSim_INS.m 1879 2008-07-15 05:20:21Z n2523710 $
%
% Filter States
% KF for position determination (PKF)
% states: 
%   x1 - Latitude
%   x2 - North Velocity
%   x3 - North acceleration
%   x4 - Longitude
%   x5 - East velocity
%   x6 - East acceleration
%   x7 - height (-ve Down position)
%   x8 - Down velocity
%   x9 - Down acceleration

% load sim data

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

% load truth data
load('data/long_flight/pos_truth_llh.mat');
load('data/long_flight/pos_truth_ecef.mat');
load('data/long_flight/vel_truth.mat'); % velocity truth in NED
load('data/long_flight/att_truth.mat');


r2d = 180 / pi;
d2r = pi / 180;

% time
Time = sensors(1,:);
NumberEpochs = length(Time);

% the INS mechanisation is done in a North-East-Down (NED) navigation frame

% setup initial parameters

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

PKF_x_hat_in = [InitialPosition(1);InitialVelocity(1);0;InitialPosition(2);InitialVelocity(2);0;InitialPosition(3);InitialVelocity(3);0];
PKF_P_in = zeros(9,9);

% quaternion representation
for Epoch=1:NumberEpochs
   if Epoch == 1
       % first time step
       a(Epoch) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
       b(Epoch) = sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) - cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
       c(Epoch) = cos(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2);
       d(Epoch) = cos(InitialAttitude(1)/2) * cos(InitialAttitude(2)/2) * sin(InitialAttitude(3)/2) + sin(InitialAttitude(1)/2) * sin(InitialAttitude(2)/2) * cos(InitialAttitude(3)/2);
       
       % initialise geodetic position
       Pos_LLH(Epoch,:) = InitialPosition;       
        Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch,2),Pos_LLH(Epoch,1));
        Tned2ecef = Tecef2ned';
        dt = 0.02;
   else
        dt = Time(Epoch) - Time(Epoch-1);
        
        % get sensor measurements
        omega_x = sensors(5,Epoch);
        omega_y = sensors(6,Epoch);
        omega_z = sensors(7,Epoch);

        
        a_dot = -0.5 * (b(Epoch-1) * omega_x + c(Epoch-1) * omega_y + d(Epoch-1) * omega_z);
        b_dot = 0.5 * (a(Epoch-1) * omega_x - d(Epoch-1) * omega_y + c(Epoch-1) * omega_z);
        c_dot = 0.5 * (d(Epoch-1) * omega_x + a(Epoch-1) * omega_y - b(Epoch-1) * omega_z);
        d_dot = -0.5 * (c(Epoch-1) * omega_x - b(Epoch-1) * omega_y - a(Epoch-1) * omega_z);
        
        a(Epoch) = a(Epoch-1) + a_dot * dt;
        b(Epoch) = b(Epoch-1) + b_dot * dt;
        c(Epoch) = c(Epoch-1) + c_dot * dt;
        d(Epoch) = d(Epoch-1) + d_dot * dt;
        
        Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch-1,2),Pos_LLH(Epoch-1,1));
        Tned2ecef = Tecef2ned';
        
   end      
        A_xb = sensors(2,Epoch);% + randn(1)/1000; 
        A_yb = sensors(3,Epoch);% + randn(1)/1000;
        A_zb = sensors(4,Epoch);% + randn(1)/1000;
        
        c11 = (a(Epoch)^2 + b(Epoch)^2 - c(Epoch)^2 - d(Epoch)^2);
        c12 = 2 * (b(Epoch)*c(Epoch) - a(Epoch)*d(Epoch));
        c13 = 2 * (b(Epoch)*d(Epoch) + a(Epoch)*c(Epoch));
        c21 = 2 * (b(Epoch)*c(Epoch) + a(Epoch)*d(Epoch));
        c22 = (a(Epoch)^2 - b(Epoch)^2 + c(Epoch)^2 - d(Epoch)^2);
        c23 = 2 * (c(Epoch)*d(Epoch) - a(Epoch)*b(Epoch));
        c31 = 2 * (b(Epoch)*d(Epoch) - a(Epoch)*c(Epoch));
        c32 = 2 * (c(Epoch)*d(Epoch) + a(Epoch)*b(Epoch));
        c33 = (a(Epoch)^2 - b(Epoch)^2 - c(Epoch)^2 + d(Epoch)^2);
        
        C_bn = [c11, c12, c13; c21, c22, c23; c31, c32, c33];
        
        phi_q(Epoch) = atan2(c32,c33);
        theta_q(Epoch) = asin(-c31);
        psi_q(Epoch) = atan2(c21,c11);
        
        % rotate accelerometer measurements to nav frame using A_n = C_bn *
        % A_b
        
        A_n = C_bn * [A_xb;A_yb;A_zb];
        A_xn = A_n(1);
        A_yn = A_n(2);
        A_zn = A_n(3) - g;
        
        % perform coriolis correction

        
        OMEGA_e_n = Tecef2ned * [0;0;(15/3600)*d2r];
        V_n = [PKF_x_hat_in(2);PKF_x_hat_in(5);PKF_x_hat_in(8)];
        Coriolis = cross(2 * OMEGA_e_n,V_n);
        A_xn = A_xn - Coriolis(1);
        A_yn = A_yn - Coriolis(2);
        A_zn = A_zn - Coriolis(3);
      
 % RECTANGULAR INTEGRATION POSITION UPDATE       
%         % propogate velocity estimate
%         V_xn(Epoch) = V_xn(Epoch-1) + A_xn * dt;
%         V_yn(Epoch) = V_yn(Epoch-1) + A_yn * dt;
%         V_zn(Epoch) = V_zn(Epoch-1) + A_zn * dt;
%         
%         % determine latitude and longitude rates
%         
%         % propogate position estimate
%         P_xn(Epoch) = P_xn(Epoch-1) + V_xn(Epoch) * dt;
%         P_yn(Epoch) = P_yn(Epoch-1) + V_yn(Epoch) * dt;
%         P_zn(Epoch) = P_zn(Epoch-1) + V_zn(Epoch) * dt;
        

 % KF POSITION UPDATE
       % KF for position determination (PKF)
        % states: 
        %   x1 - Latitude
        %   x2 - North Velocity
        %   x3 - North acceleration
        %   x4 - Longitude
        %   x5 - East velocity
        %   x6 - East acceleration
        %   x7 - height (-ve Down position)
        %   x8 - Down velocity
        %   x9 - Down acceleration
        
        % PHI matrix - x_hat_k_minus_minus = PHI_k * x_hat_k_minus1_plus
        PKF_PHI = eye(9,9);
        % x - channel
        PKF_PHI(1,2) = dt/(RM+PKF_x_hat_in(7));
        PKF_PHI(1,3) = 0;
        PKF_PHI(2,3) = dt;
        % y - channel
        PKF_PHI(4,5) = dt/((RP+PKF_x_hat_in(7))*cos(PKF_x_hat_in(1)));
        PKF_PHI(4,6) = 0;
        PKF_PHI(5,6) = dt;
        % z - channel
        PKF_PHI(7,8) = -dt;
        PKF_PHI(7,9) = 0;
        PKF_PHI(8,9) = dt;
        
        
        PKF_z = [A_xn;A_yn;A_zn];
        PKF_H = zeros(3,9);
        PKF_H(1,3) = 1;
        PKF_H(2,6) = 1;
        PKF_H(3,9) = 1;
        
        % Process Noise Matrix, Q
        dt_2 = (dt^2) / 2;
        dt_3 = (dt^3) / 3;
        PKF_Q = zeros(9,9);
        PKF_Q(1,1) = dt_3;
        PKF_Q(2,2) = dt_2;
        PKF_Q(3,3) = dt;
        PKF_Q(4,4) = dt_3;
        PKF_Q(5,5) = dt_2;
        PKF_Q(6,6) = dt;
        PKF_Q(7,7) = dt_3;
        PKF_Q(8,8) = dt_2;
        PKF_Q(9,9) = dt;
        
        
        % Measurement Noise Matrix, R
        PKF_R = eye(3,3) / 1000;  % 
        
        
        [PKF_x_hat_out, PKF_P_out] = GARD_EvaluateKF(dt, PKF_x_hat_in, PKF_P_in, PKF_PHI, PKF_H, PKF_z, PKF_Q, PKF_R);
        
        % update for next step
        PKF_x_hat_in = PKF_x_hat_out;
        PKF_P_in = PKF_P_out;


   % save results for analysis
        PKF_x_save(Epoch,:) = PKF_x_hat_out;
        PKF_P_save(Epoch,:) = diag(PKF_P_out);
        
        
        % find the NED to ECEF transfer matrix and rotate position vector
        % to ECEF frame, then add to ECEF position.

        
        
        %Pos_ECEF(Epoch,:) = Pos_ECEF(1,:)' + Tned2ecef * [PKF_x_hat_out(1);PKF_x_hat_out(4);PKF_x_hat_out(7)];
        %[Pos_LLH(Epoch,1),Pos_LLH(Epoch,2),Pos_LLH(Epoch,3)] = ECEF2LLH(Pos_ECEF(Epoch,:));
        Pos_LLH(Epoch,:) = [PKF_x_hat_out(1),PKF_x_hat_out(4),PKF_x_hat_out(7)];
        
        

end

% save results to file for later processing
save data/long_flight/ins_output PKF_x_save PKF_P_save phi_q theta_q psi_q

% plot attitude determination results
% 
figure();
%plot(Time,att_truth(2,:)*r2d,'b');
hold on;
plot(Time,(att_truth(2,:) - phi_q)*r2d,'r');
%legend('KF','rect');
title('Roll Angle Error');
% 
figure();
%plot(Time,theta*r2d,'b+');
hold on;
plot(Time,(att_truth(3,:) - theta_q)*r2d,'r');
%legend('KF','rect');
title('Pitch Angle Error');

figure();
%plot(Time,psi*r2d,'b+');
hold on;
plot(Time,(att_truth(4,:) - psi_q)*r2d,'r');
%legend('KF','rect');
title('Yaw Angle Error');
% 
% figure();
% hold on;
% plot(PKF_x_save(:,4),PKF_x_save(:,1),'r+');
% plot(P_yn,P_xn,'go');
% grid on;
% xlabel('Easting (m)');
% ylabel('Northing (m)');
% 
% figure();
% hold on;
% plot(PKF_x_save(:,7),'r+');
% plot(P_zn,'go');
% grid on;
% xlabel('Time (s)');
% ylabel('Height (m)')'

% % plot results in LLH
figure();
hold on;
plot(pos_truth_llh(3,:)*r2d,pos_truth_llh(2,:)*r2d);
plot(Pos_LLH(:,2)*r2d,Pos_LLH(:,1)*r2d,'r+');
grid on;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
hold off;
% 
% % plot latitude error
figure();
hold on;
plot(Time, (Pos_LLH(:,1)' - pos_truth_llh(2,:)) * RM);
% errror bounds
plot(Time, sqrt(PKF_P_save(:,1) * RM),'r');
plot(Time, -sqrt(PKF_P_save(:,1) * RM),'r');
grid on;
xlabel('Simulation Time (seconds)');
ylabel('Latitude Error (metres North-South)');
hold off;

% plot longitude
figure();
hold on;
plot(Time, (Pos_LLH(:,2) - pos_truth_llh(3,:)')* RP .* cos(Pos_LLH(:,1))); % note - correcting RP approxiimately for latitude - i think this is right??
% errror bounds
plot(Time, sqrt(PKF_P_save(:,4) * RP .* cos(Pos_LLH(:,1))),'r'); 
plot(Time, -sqrt(PKF_P_save(:,4) * RP .* cos(Pos_LLH(:,1))),'r');
grid on;
xlabel('Simulation Time (seconds)');
ylabel('Longitude Error (metres East-West)');
hold off;

% % plot height 
figure();
hold on;
plot(Time, Pos_LLH(:,3)' - pos_truth_llh(4,:));
grid on;
xlabel('Simulation Time (seconds)');
ylabel('Height Error (m)');
hold off;

% plot velocity errors
%North Velocity
figure();
hold on;
plot(PKF_x_save(:,2) - vel_truth(2,:)');
% plot error bounds
plot(2*sqrt(PKF_P_save(:,2)),'r');
plot(-2*sqrt(PKF_P_save(:,2)),'r');
hold off;
grid on;
xlabel('Sample Number');
ylabel('North Velocity Error (m/s)');
% 
figure();
hold on;
plot(PKF_x_save(:,5) - vel_truth(3,:)');
% plot error bounds
plot(2*sqrt(PKF_P_save(:,5)),'r');
plot(-2*sqrt(PKF_P_save(:,5)),'r');
hold off;
grid on;
xlabel('Sample Number');
ylabel('East Velocity Error (m/s)');
% 
figure();
hold on;
plot(PKF_x_save(:,8) - vel_truth(4,:)');
% plot error bounds
plot(2*sqrt(PKF_P_save(:,8)),'r');
plot(-2*sqrt(PKF_P_save(:,8)),'r');
hold off;
grid on;
xlabel('Sample Number');
ylabel('Down Velocity Error (m/s)');


