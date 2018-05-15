% Kalman filtered GPS position and velocity solution using an 8 state
% Kalman Filter
%
% 
% implements the 8 State GPS position and velocity kalman filter solution
% written by Peter Roberts
%
% adapted by duncan greer 27 September 2005
% $Id: GARD_KFGPSSolution8.m 1850 2008-07-14 04:52:47Z greerd $
%

% load gps constants
GPSConstants;

% this is still not quite right as it does not incorporate the range rate
% measurement properly.  see the blue book.

dt = 1;   % 1 second update rate

Sigma_PR = 10; % m - standard deviation of expected pseudorange noise

% starting position estimate  - start point of flight data
% use LLH2ECEF(-27.40*pi/180, 153.20*pi/180, 4500/3.28)                   
X_init = -5059037.69670272;
Y_init = 2555503.82695979;
Z_init = -2918268.17673033;
RxClockBias_init = 0;
          
UserPos = [X_init, Y_init, Z_init, RxClockBias_init];

% T_ECEF2NED(-27.4*pi/180,153.2*pi/180)' * [0 -60 0]'
UserVel=[  -27.61198,-53.26892, 0,0]; % 50 m/s north
   
% setup error values
PositionError = 5.0;
VelocityError = 2.0;
ClockBiasError = 1E5;
ClockDriftError = 10.0;
RangeNoiseVariance = 7.5^2; % m
RangeRateNoiseVariance = 5.0^2; % m/s

% not sure what these values are - something to do with the process noise
Sp = 1.0;
Sf = 1.0;
Sg = 1.0;

% setup state transition - time invariant if dt is constant
phi = eye(8,8);
phi(1,5) = dt;
phi(2,6) = dt;
phi(3,7) = dt;
phi(4,8) = dt;


% setup Q and R matrices (process and noise covariance)
Q = zeros(8,8);

Sp_dt3 = Sp * (dt ^ 3) / 3;
Sp_dt2 = Sp * (dt ^ 2) / 2;
Sp_dt = Sp * dt;
Sf_dt = Sf * dt;
Sg_dt3 = Sg * (dt ^ 3) / 3;
Sg_dt2 = Sg * (dt ^ 2) / 2;
Sg_dt = Sg * dt;

Q(1,1) = Sp_dt3;
Q(2,2) = Sp_dt3;
Q(3,3) = Sp_dt3;

Q(1,5) = Sp_dt2;
Q(2,6) = Sp_dt2;
Q(3,7) = Sp_dt2;
Q(4,8) = Sg_dt2;
Q(5,1) = Sp_dt2;
Q(6,2) = Sp_dt2;
Q(7,3) = Sp_dt2;
Q(8,4) = Sg_dt2;

Q(4,4) = Sf_dt + Sg_dt3;
Q(5,5) = Sp_dt;
Q(6,6) = Sp_dt;
Q(7,7) = Sp_dt;
Q(8,8) = Sg_dt;


% inialise x_hat_in, P_in
% x_hat_in = [X_init, 
%             Y_init, 
%             Z_init, 
%             RxClockBias_init,
%             50,
%             0,
%             0,
%             0]

x_hat_in = [UserPos, UserVel]';
% REMEMBER x_hat is the DELTA VALUE!!!!!!
%x_hat_in = [0 0 0 0 0 0 0 0]';

PositionVariance = PositionError ^ 2;
VelocityVariance = VelocityError ^ 2;
ClockBiasVariance = ClockBiasError ^ 2;
ClockDriftVariance = ClockDriftError ^2;
%
P_in = zeros(8,8);
P_in_dash = [PositionVariance,
        PositionVariance,
        PositionVariance,
        ClockBiasVariance,
        VelocityVariance,
        VelocityVariance,
        VelocityVariance,
        ClockDriftVariance];
for index = 1:8
    P_in(index,index) = P_in_dash(index);
end

    

load 'data/rnav_approach/PR_Simulation.mat';
NumberEpochs = length(PR_Sim);
NumberSVs = 32;

% iono model parameters for the above Nav file
ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA           
BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA  

% clear bogus pseudorange rate measurements
% for SV=1:NumberSVs
%    PRR_Sim(1,SV) = 0;
% end

% begin loop
for Epoch = 1:NumberEpochs
    GPSTime  = 259199 + Epoch;
    % clear variables which are dependant on the number of measurements
    clear R SVPos H;        
    % get the measurements vector
        
    % format the input vectors
    SVIndex = 0;
    for SV=1:NumberSVs
        if(PR_Sim(Epoch,SV) ~= 0)
            % add to PR vector
            SVIndex = SVIndex + 1;
            PR_Vec(SVIndex) = PR_Sim(Epoch,SV);
            PRR_Vec(SVIndex) = PRR_Sim(Epoch,SV);
%             if(Epoch == 1)
%                 PRR_Vec(SVIndex) = 0;
%             else
%                 if(CP_Sim(Epoch-1,SV) == 0)
%                     PRR_Vec(SVIndex) = (CP_Sim(Epoch+1,SV) - CP_Sim(Epoch,SV));
%                 else
%                     PRR_Vec(SVIndex) = (CP_Sim(Epoch,SV) - CP_Sim(Epoch-1,SV));
%                 end
%             end
            SVPos(SVIndex,:) = [SV_X_Data(Epoch,SV) SV_Y_Data(Epoch,SV) SV_Z_Data(Epoch,SV) SV_T_Data(Epoch,SV)];
            SVVel(SVIndex,:) = [SV_Xvel_Data(Epoch,SV) SV_Yvel_Data(Epoch,SV) SV_Zvel_Data(Epoch,SV) SV_Tvel_Data(Epoch,SV)];
            SVAcc(SVIndex,:) = [SV_Xacc_Data(Epoch,SV) SV_Yacc_Data(Epoch,SV) SV_Zacc_Data(Epoch,SV) SV_Tacc_Data(Epoch,SV)];
            

            % calculate hte earth rotation correction as per Kayton pg 228
            % eq 5.67

            delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) * UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
            PR_Vec(SVIndex) = PR_Vec(SVIndex) + delta_pr_omegaedot;
            
            % calculate the iono delay correction - single frequency user
            % model from ICD 200
            D_IONO = ionomodel(GPSTime, UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
            PR_Vec(SVIndex) = PR_Vec(SVIndex) - D_IONO;
            
        end
    end

    NumberMeasurements = length(PR_Vec);

    % R matrix
    R = eye(2*NumberMeasurements);
    for k = 1:NumberMeasurements
        R(k,k) = RangeNoiseVariance;
        R(k+NumberMeasurements,k+NumberMeasurements) = RangeRateNoiseVariance;
    end
    
    % determine the H matrix
    for k = 1:NumberMeasurements
        %Calculated slant ranges
    
        for m = 1:3
             ele(m) =  SVPos(k,m) - UserPos(m);
        end    
                
        r_VecCalc(k) =  norm(ele);   

        H(k,1) =  -ele(1)/r_VecCalc(k);
        H(k,2) =  -ele(2)/r_VecCalc(k);
        H(k,3) =  -ele(3)/r_VecCalc(k);
        H(k,4) = 1;   
        H(k,5) = 0;
        H(k,6) = 0;
        H(k,7) = 0;
        H(k,8) = 0;
        
        H(k+NumberMeasurements,1) = 0;
        H(k+NumberMeasurements,2) = 0;
        H(k+NumberMeasurements,3) = 0;
        H(k+NumberMeasurements,4) = 0;   
        H(k+NumberMeasurements,5) = -ele(1)/r_VecCalc(k);
        H(k+NumberMeasurements,6) = -ele(2)/r_VecCalc(k);
        H(k+NumberMeasurements,7) = -ele(3)/r_VecCalc(k);
        H(k+NumberMeasurements,8) = 1;
        
        % find apriori estimate of pseudorange
        PR_Vec_minus(k) = r_VecCalc(k) + UserPos(4) + UserVel(4) * dt;  % geometric range + c * delta_T
        
        %predicted relative velocity of sv and receiver
        r_VecCalcVel(k) = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
        Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);

        PRR_Vec_minus(k) = Relative_Velocity(k) + UserVel(4);
    end
    
    % find measurement vector - delta between predicted and measured
    % pseudorange
    z_Vec = ([PR_Vec PRR_Vec]' - [PR_Vec_minus PRR_Vec_minus]');
    
    
    % evaluate the KF
    %[x_hat_out, P_out] = GARD_EvaluateKF(dt, x_hat_in, P_in, phi, H, z_Vec, Q, R);
    % Prediction step
    x_hat_minus = phi * x_hat_in;
    
    P_minus = phi * P_in * phi' + Q;

    % correction step

    % calculate the kalman gain
    V = H * P_minus * H' + R;
    K = P_minus * H' * inv(V);

    % update equations
    %v = z - H * x_hat_minus;
    x_hat_out = x_hat_minus + K * (z_Vec);
    P_out = P_minus - K * H * P_minus;
    
    % save the filter covariance
    P_save(Epoch,:) = diag(P_out);
    
    % update user position for next iteration
      %UserPos = UserPos + x_hat_out(1:4)';
      %UserVel = UserVel + x_hat_out(5:8)';
    UserPos = x_hat_out(1:4)';
    UserVel = x_hat_out(5:8)';
    
    x_hat_in = x_hat_out;
    P_in = P_out;
    

    
    % save position history
    x_out(Epoch,:) = [UserPos UserVel];

    % clear variables
    clear PR_Vec PRR_Vec SVPos SVVel SVAcc PR_Vec_minus PRR_Vec_minus;
end

% convert positions to llh
for Epoch=1:NumberEpochs
    [Latitude(Epoch),Longitude(Epoch),Height(Epoch)] = ECEF2LLH(x_out(Epoch,1:3));
    T = T_ECEF2ENU(Longitude(Epoch),Latitude(Epoch));
    % convert to degrees
    Latitude(Epoch) = Latitude(Epoch) * 180 / pi;
    Longitude(Epoch) = Longitude(Epoch) * 180 / pi;
    
    % rotate velocity values to ENU
    Velocity(Epoch,:) = T * x_out(Epoch,5:7)';
end

gps_out_ecef = x_out;
gps_out_llh = [Latitude; Longitude; Height]';

% save data to file
save data/rnav_approach/gps_output gps_out_ecef gps_out_llh;

% output final position

load data/rnav_approach/pos_truth_ecef.mat;
load data/rnav_approach/pos_truth_llh.mat;
load data/rnav_approach/vel_truth.mat;

for Epoch=1:NumberEpochs
    X_error(Epoch) = pos_truth_ecef(2,Epoch * 100) - x_out(Epoch,1);
    Y_error(Epoch) = pos_truth_ecef(3,Epoch * 100) - x_out(Epoch,2);
    Z_error(Epoch) = pos_truth_ecef(4,Epoch * 100) - x_out(Epoch,3);
    Nvel_error(Epoch) = vel_truth(2,Epoch * 100) - Velocity(Epoch,1);
    Evel_error(Epoch) = vel_truth(3,Epoch * 100) - Velocity(Epoch,2);
    Dvel_error(Epoch) = vel_truth(4,Epoch * 100) - Velocity(Epoch,3);
end

%plot(X_error);

% calculate the mean and variance of the final errors
X_error_mean = mean(X_error);
X_error_var = var(X_error);
Y_error_mean = mean(Y_error);
Y_error_var = var(Y_error);
Z_error_mean = mean(Z_error);
Z_error_var = var(Z_error);

Nvel_error_mean = mean(Nvel_error);
Nvel_error_var = var(Nvel_error);
Evel_error_mean = mean(Evel_error);
Evel_error_var = var(Evel_error);
Dvel_error_mean = mean(Dvel_error);
Dvel_error_var = var(Evel_error);

X_noise = sqrt(X_error_mean^2 + X_error_var)
Y_noise = sqrt(Y_error_mean^2 + Y_error_var)
Z_noise = sqrt(Z_error_mean^2 + Z_error_var)
N_noise = sqrt(Nvel_error_mean^2 + Nvel_error_var)
E_noise = sqrt(Evel_error_mean^2 + Evel_error_var)
D_noise = sqrt(Dvel_error_mean^2 + Dvel_error_var)



r2d = 180/pi;

% plot position output
figure();
hold on;
plot(pos_truth_llh(3,:)*r2d,pos_truth_llh(2,:)*r2d,'g')
plot(Longitude,Latitude,'r+')
hold off;
grid on;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

% position error ecef
figure();
subplot(3,1,1),plot(X_error);
hold on;
subplot(3,1,1),plot(2*sqrt(P_save(:,1)),'r');
subplot(3,1,1),plot(-2*sqrt(P_save(:,1)),'r');
ylabel('X Position Error');
subplot(3,1,2),plot(Y_error);
hold on;
subplot(3,1,2),plot(2*sqrt(P_save(:,2)),'r');
subplot(3,1,2),plot(-2*sqrt(P_save(:,2)),'r');
ylabel('Y Position Error');
subplot(3,1,3),plot(Z_error);
hold on;
subplot(3,1,3),plot(2*sqrt(P_save(:,3)),'r');
subplot(3,1,3),plot(-2*sqrt(P_save(:,3)),'r');
xlabel('Time (sec)');
ylabel('Z Position Error');

% plot velocity error output

figure();
subplot(3,1,1),plot(Nvel_error);
axis([0 length(Velocity) -2 2])
ylabel('North Velocity Error (m)');
subplot(3,1,2),plot(Evel_error);
axis([0 length(Velocity) -2 2])
ylabel('East Velocity Error (m)');
subplot(3,1,3),plot(Dvel_error);
axis([0 length(Velocity) -2 2])
xlabel('Time (sec)');
ylabel('Down Velocity Error (m)');




% plot North Velocity
figure();
plot(vel_truth(1,:),vel_truth(2,:))
hold on
plot(Velocity(:,1),'ro')
grid on
xlabel('Epoch (s)');
ylabel('North Velocity (m/s)');
legend('Simulated Truth','KF Estimate');

% East Velocity
figure();
hold on;
plot(vel_truth(1,:),vel_truth(3,:));
plot(Velocity(:,2),'ro');
grid on;
xlabel('Epoch (s)');
ylabel('East Velocity (m/s)');
legend('Simulated Truth','KF Estimate');

% Down velocity
figure();
hold on;
plot(vel_truth(1,:),vel_truth(4,:));
plot(Velocity(:,3),'ro');
xlabel('Epoch (s)');
ylabel('Down Velocity (m/s)');
legend('Simulated Truth','KF Estimate');
