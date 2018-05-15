%% GPS Position Solution using UKF
% Written by Duncan Greer 2 Feb 2007
%
% $Id: GARDSim_GPS_UKF.m 1879 2008-07-15 05:20:21Z n2523710 $
%
%



% x1  = ECEF X
% x2  = ECEF Y
% x3  = ECEF Z
% x4  = ECEF X Velocity
% x5  = ECEF Y Velocity
% x6  = ECEF Z Velocity
% x7  = Receiver Clock Bias (Rb)
% x8  = Receiver Clock Drift Rate (Rf)
% x9 = ECEF X Acceleration
% x10 = ECEF Y Acceleration
% x11 = ECEF Z Acceleration


%% setup simulation

% set number formatting for display
format long g;

% load GPS constants
GPSConstants;

% iono model parameters for the above Nav file
ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA           
BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA  

% position
InitialPosition = [-27.38*pi/180 153.12*pi/180 300]';  % somewhere near brisbane in LLH, rads and meters
% velocity
InitialVelocity = [50 0 0]';  % NED velocities in m/s
% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(2));
RP = PrimeRadius(InitialPosition(2));

% 3 position, 3 velocity and 2 clock
NumberStates = 11;

% process noise states - 8
ProcessNoiseStates = 11;  

% use a maximum of 8 measurements (PR+PRR)
MeasurementNoiseStates=16;


Na = NumberStates+ProcessNoiseStates+MeasurementNoiseStates;

% UKF scaling parameters
alpha = 0.8; % range: 1e-3 < alpha <= 1
beta = 2;  % 2 is optimal for gaussian priors
kapa = 0;  %

LAMBDA = alpha^2 * (Na+kapa) - Na;

% load truth data
load 'data/long_flight/pos_truth_llh.mat';
load 'data/long_flight/vel_truth.mat';
load 'data/long_flight/att_truth.mat';

% load GPS solution data
load 'data/long_flight/PR_Simulation.mat';


NumberEpochsGPS = length(PR_Sim);
gps_dt = 1;

TimeGPS = [0:gps_dt:NumberEpochsGPS-1];

NumberSVs = 32;
SVDontUse = zeros(1,NumberSVs);

%% initialise covariances
Px_kminus = eye(NumberStates,NumberStates)*100;
Q = eye(ProcessNoiseStates,ProcessNoiseStates)*10*10;

R = eye(MeasurementNoiseStates,MeasurementNoiseStates);
R(1:8,1:8) = eye(8,8)*7.5*7.5;
R(9:16,9:16) = eye(8,8)*3*3;

x_hat_kminus = zeros(11,1);
x_hat_kminus(1:3,1) = llh2ecef(-27*pi/180,150*pi/180,0);
GRASOn = 0;


%% start gps solutions
for Epoch_lo = 1:NumberEpochsGPS

    % get the gps time - used for the iono correction
    GPSTime  = 259199 + Epoch_lo;


    % perform state prediction
    % generate augmented state vector
    xa_hat_kminus = [x_hat_kminus;zeros(ProcessNoiseStates,1);zeros(MeasurementNoiseStates,1)];

    Pxa_kminus = [Px_kminus          zeros(NumberStates,NumberStates) zeros(NumberStates,MeasurementNoiseStates);
              zeros(ProcessNoiseStates,NumberStates) Q           zeros(ProcessNoiseStates,MeasurementNoiseStates);
              zeros(MeasurementNoiseStates,NumberStates) zeros(MeasurementNoiseStates,ProcessNoiseStates)  R];
          
          
          
          
          
    blah = sqrt(Na+LAMBDA) * chol(Pxa_kminus);
    for i=0:Na
        if i==0
            xs_0 = xa_hat_kminus;
            W_0_m = LAMBDA / (Na + LAMBDA);
            W_0_c = W_0_m + (1 - alpha^2 + beta);
        else
            xs_i(:,i) = xa_hat_kminus + blah(:,i);
            xs_i(:,i+Na) = xa_hat_kminus - blah(:,i);
            W_i_m = 1 / (2 * (Na + LAMBDA));
            W_i_c = W_i_m;
        end
    end

    
    
    
    
    % propagate system dynamics
    for i=0:2*Na
        if i==0
            %% zeroth point
            
            % position
            xs_0_k(1) = xs_0(1) + xs_0(4)*gps_dt;
            xs_0_k(2) = xs_0(2) + xs_0(5)*gps_dt;
            xs_0_k(3) = xs_0(3) + xs_0(6)*gps_dt;
            xs_0_k(7) = xs_0(7) + xs_0(8)*gps_dt;
            
            % velocity
            xs_0_k(4) = xs_0(4) + xs_0(9)*gps_dt;
            xs_0_k(5) = xs_0(5) + xs_0(10)*gps_dt;
            xs_0_k(6) = xs_0(6) + xs_0(11)*gps_dt;
            xs_0_k(8) = xs_0(8);
            
            % acceleration
            xs_0_k(9) = xs_0(9);
            xs_0_k(10) = xs_0(10);
            xs_0_k(11) = xs_0(11);
            
            
            
            %% add process noise
            V = xs_0(NumberStates+1:NumberStates+ProcessNoiseStates)';
            xs_0_k(1:NumberStates) = xs_0_k(1:NumberStates) + V;
            
            
        else
            
            % i-points
            
            % position
            xs_i_k(1,i) = xs_i(1,i) + xs_i(4,i)*gps_dt;
            xs_i_k(2,i) = xs_i(2,i) + xs_i(5,i)*gps_dt;
            xs_i_k(3,i) = xs_i(3,i) + xs_i(6,i)*gps_dt;
            xs_i_k(7,i) = xs_i(7,i) + xs_i(8,i)*gps_dt;
            
            % velocity
            xs_i_k(4,i) = xs_i(4,i) + xs_i(9,i)*gps_dt;
            xs_i_k(5,i) = xs_i(5,i) + xs_i(10,i)*gps_dt;
            xs_i_k(6,i) = xs_i(6,i) + xs_i(11,i)*gps_dt;
            xs_i_k(8,i) = xs_i(8,i);
            
            % acceleration
            xs_i_k(9,i) = xs_i(9,i);
            xs_i_k(10,i) = xs_i(10,i);
            xs_i_k(11,i) = xs_i(11,i);
            
            
            
            %% add process noise
            V = xs_i(NumberStates+1:NumberStates+ProcessNoiseStates,i);
            xs_i_k(1:NumberStates,i) = xs_i_k(1:NumberStates,i) + V;
        end
    end
    
    
    
    
    
    
    
    
    W = eye(1,2*Na+1);
    W(1,1) = W_0_m;
    W(1,2:2*Na+1) = W_i_m;
    X = zeros(NumberStates,2*Na+1);
    X(:,1) = xs_0_k(1:NumberStates);
    X(:,2:2*Na+1) = xs_i_k(1:NumberStates,:);
    x_hat_kminus = (W*X')';

    nx = length(x_hat_kminus); nw = length(W);
    Px_kminus=((ones(nx,1)*W).*(X-x_hat_kminus*ones(1,nw)))*(X-x_hat_kminus*ones(1,nw))';
        
    % initial guess at user pos required to estimate earth rotation
    % correction
    UserPos(1:3) = [x_hat_kminus(1);x_hat_kminus(2);x_hat_kminus(3)];
 %   UserPos(4) = xs_0_k(7);
        
    % format the input measurement vectors
    SVIndex = 0;
    for SV=1:NumberSVs
        if((PR_Sim(Epoch_lo,SV) > 100) && SVDontUse(SV) == 0)  %%SV ~= ExcludeSVPRN
            % add to PR vector
            SVIndex = SVIndex + 1;

            SV_Vec(SVIndex) = SV;
            PR_Vec(SVIndex) = PR_Sim(Epoch_lo,SV);
            %PRR_Vec(SVIndex) = PRR_Sim(Epoch_lo,SV);
            if(Epoch_lo == 1)
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
            % save previous value
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

            % calculate the iono delay correction - single frequency user
            % model from ICD 200
            if(GRASOn == 0)
                ionodelay = ionomodel(GPSTime, UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
                PR_Vec(SVIndex) = PR_Vec(SVIndex) - ionodelay;
            end

        end
    end

    NumberGPSMeasurements = length(PR_Vec);



    if NumberGPSMeasurements > 8
        NumberGPSMeasurements = 8;
    end

    % formulate the measurment vector
    y_k = [PR_Vec(1:8)';PRR_Vec(1:8)'];


    for k = 1:NumberGPSMeasurements
        % find apriori estimate of pseudorange for each sigma point

        % zero-th sigma point

        
        % Get User position in ECEF
        UserPos(1:3) = [xs_0_k(1);xs_0_k(2);xs_0_k(3)];
        UserPos(4) = xs_0_k(7);
        UserVel(1:3) =  [xs_0_k(4);xs_0_k(5);xs_0_k(6)];
        UserVel(4) = xs_0_k(8);
        
        geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
        geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));

    %% calculate measurement predicition
        % pseudorange prediction
        PR_Vec_minus_0(k) = geo_range_to_sat + UserPos(4) + UserVel(4) * gps_dt;  % geometric range + c * delta_T
        %predicted relative velocity of sv and receiver
        Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
        PRR_Vec_minus_0(k) = Relative_Velocity(k) + UserVel(4);

        for i=1:2*Na
            % Get User position in ECEF
            UserPos(1:3) = [xs_i_k(1,i);xs_i_k(2,i);xs_i_k(3,i)];
            UserPos(4) = xs_i_k(7,i);
            UserVel(1:3) =  [xs_i_k(4,i);xs_i_k(5,i);xs_i_k(6,i)];
            UserVel(4) = xs_i_k(8,i);
            
            
            geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
            geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));

            PR_Vec_minus_i(k,i) = geo_range_to_sat + UserPos(4) + UserVel(4) * gps_dt + xs_i(NumberStates+ProcessNoiseStates+k-1,i);  % geometric range + c * delta_T

            %predicted relative velocity of sv and receiver
            Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
            PRR_Vec_minus_i(k,i) = Relative_Velocity(k) + UserVel(4) + xs_i(NumberStates+ProcessNoiseStates+8+k-1,i);
        end % for i=1:2*Na



    end  % for k = 1:NumberGPSMeasurements

    ys_kminus_0 = [PR_Vec_minus_0'; PRR_Vec_minus_0'];
    ys_kminus_i = [PR_Vec_minus_i(:,:); PRR_Vec_minus_i(:,:)];

    % find the sum of sigma points for the measurement prediction
    y_hat_kminus = zeros(2*NumberGPSMeasurements,1);

    w_im = 1;  % TODO: fix this value
    y_hat_kminus = y_hat_kminus + w_im * ys_kminus_0;


    for i=1:2*Na
        y_hat_kminus = y_hat_kminus + w_im * ys_kminus_i(:,i);
    end


    W = ones(1,2*Na+1);
    W(1,1) = W_0_m;
    W(1,2:2*Na+1) = W_i_m;

    Y = zeros(MeasurementNoiseStates,2*Na+1);
    Y(:,1) = ys_kminus_0(1:MeasurementNoiseStates);
    Y(:,2:2*Na+1) = ys_kminus_i(1:MeasurementNoiseStates,:);
    y_hat_kminus = (W*Y')';

    nx = length(x_hat_kminus); ny = length(y_hat_kminus); nw = length(W);
    Py_kminus=((ones(ny,1)*W).*(Y-y_hat_kminus*ones(1,nw)))*(Y-y_hat_kminus*ones(1,nw))';
    Pxy_kminus =((ones(nx,1)*W).*(X-x_hat_kminus*ones(1,nw)))*(Y-y_hat_kminus*ones(1,nw))';


    % calculate kalman gain
    K_k = Pxy_kminus * inv(Py_kminus);
    %K_k = (Pxy_kminus / Sy_kminus') / Sy_kminus;

    % apply correction
    z_k = y_k - y_hat_kminus;
    z_save(Epoch_lo,:) = z_k;
    

    
    x_hat_kplus = x_hat_kminus + K_k * (z_k);

    Px_kplus = Px_kminus - K_k * Py_kminus * K_k';


    %save results
    t_save(Epoch_lo) = Epoch_lo;
    x_hat_save(Epoch_lo,:) = x_hat_kplus;
    P_save(Epoch_lo,:) = diag(Px_kplus);


    % copy updated state to high-speed loop state for next round
    x_hat_kminus = x_hat_kplus;
    Px_kminus = Px_kplus;
    
    
        
    %% convert results to LLH
    [x_save_llh(Epoch_lo,1),x_save_llh(Epoch_lo,2),x_save_llh(Epoch_lo,3)] = ECEF2LLH(x_hat_save(Epoch_lo,1:3));
    Tecef2ned = T_ECEF2NED(x_save_llh(Epoch_lo,2),x_save_llh(Epoch_lo,1));
    
    
    P_save_llh(Epoch_lo,1:3) = abs(Tecef2ned*P_save(Epoch_lo,1:3)');
    
    
    %% calculate velocity in North-East-Down

    x_save_vel(Epoch_lo,:) = Tecef2ned*x_hat_save(Epoch_lo,4:6)';
    P_save_vel(Epoch_lo,:) = abs(Tecef2ned*P_save(Epoch_lo,4:6)');
    
    %% calculate position error
    pos_error_llh(Epoch_lo,:) = x_save_llh(Epoch_lo,:) - pos_truth_llh(2:4,Epoch_lo*100)';
    pos_error_llh(Epoch_lo,1) = pos_error_llh(Epoch_lo,1) * RM;
    pos_error_llh(Epoch_lo,2) = pos_error_llh(Epoch_lo,2) * RP * cos(pos_truth_llh(2,Epoch_lo*100));
    
    
    vel_error_ned(Epoch_lo,:) = x_save_vel(Epoch_lo,:) - vel_truth(2:4,Epoch_lo*100)';
    
    x_save_acc(Epoch_lo,:) = Tecef2ned*x_hat_save(Epoch_lo,9:11)';
    
    %% update R-matrix
%     if(Epoch_lo > 10)
%         for k = 1:NumberGPSMeasurements
%             R(k,k) = std(abs(z_save(:,k)))^2;
%             R(k+NumberGPSMeasurements,k+NumberGPSMeasurements) = std(abs(z_save(:,k+NumberGPSMeasurements)))^2;
%         end
%     end
%     
    
    
    
    
    disp(sprintf('Completed Epoch %d',Epoch_lo));

    
    
    
end



%% plot results
% figure();
% plot3(x_hat_save(:,1),x_hat_save(:,2),x_hat_save(:,3));

figure();
plot(x_save_llh(:,2)*180/pi,x_save_llh(:,1)*180/pi,'r*');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
grid on;
title('GPS UKF Position Solution');

% figure();
% subplot(3,1,1),plot(t_save,x_save_llh(:,1));
% subplot(3,1,2),plot(t_save,x_save_llh(:,2));
% subplot(3,1,3),plot(t_save,x_save_llh(:,3));

figure();

subplot(3,1,1),plot(t_save,pos_error_llh(:,1));
title('Position Error');
hold on;
subplot(3,1,1),plot(t_save,2*sqrt(P_save_llh(:,1)),'r');
subplot(3,1,1),plot(t_save,-2*sqrt(P_save_llh(:,1)),'r');

subplot(3,1,2),plot(t_save,pos_error_llh(:,2));
hold on;
subplot(3,1,2),plot(t_save,2*sqrt(P_save_llh(:,2)),'r');
subplot(3,1,2),plot(t_save,-2*sqrt(P_save_llh(:,2)),'r');

subplot(3,1,3),plot(t_save,pos_error_llh(:,3));
hold on;
subplot(3,1,3),plot(t_save,2*sqrt(P_save_llh(:,3)),'r');
subplot(3,1,3),plot(t_save,-2*sqrt(P_save_llh(:,3)),'r');


figure();
subplot(3,1,1),plot(t_save,x_save_vel(:,1));
title('Velocity');
subplot(3,1,2),plot(t_save,x_save_vel(:,2));
subplot(3,1,3),plot(t_save,x_save_vel(:,3));

figure();
subplot(3,1,1),plot(t_save,vel_error_ned(:,1));
title('Velocity Error');
subplot(3,1,2),plot(t_save,vel_error_ned(:,2));
subplot(3,1,3),plot(t_save,vel_error_ned(:,3));

figure();
subplot(3,1,1),plot(t_save,x_save_acc(:,1));
title('Acceleration');
subplot(3,1,2),plot(t_save,x_save_acc(:,2));
subplot(3,1,3),plot(t_save,x_save_acc(:,3));





