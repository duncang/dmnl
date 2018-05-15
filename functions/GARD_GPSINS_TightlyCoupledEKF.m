function results = GARD_GPSINS_TightlyCoupledEKF(imudata, gps,barodata,magdata, ...
    IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE, Sensor_Params, CompassCalibration, StartTime, StopTime, InitialPosition, InitialVelocity,InitialAttitude,GPSLeverArm)

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


ModuleName = '[GARD_GPSINS_TightlyCoupledEKF] ';


% setup constants
WGS84Constants;
d2r = pi/180;
r2d = 180/pi;
g = -9.79;  % m/s/s
GPSConstants;

NumberStates = 17;






gps_dt = 1/GPS_RATE;
ins_dt = 1/IMU_RATE;



% initialise epochs
Epoch_lo = 0;
Epoch_hi = 0;


C_BN = GARD_EulerToDCM(InitialAttitude(1),InitialAttitude(2),InitialAttitude(3));
% 

% get the initial meridian and prime radii of curvature
RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));
RMh = RM + InitialPosition(3);
RPh = RP + InitialPosition(3);


NumberSVs = 32;

SVDontUse = zeros(1,NumberSVs);

SVDontUse(29) = 1;

NumberEpochsGPS = StopTime * GPS_RATE;

TimeINS = 0:ins_dt:NumberEpochsGPS;
TimeGPS = 0:gps_dt:NumberEpochsGPS-1;


NumberEpochsINS = min(length(TimeINS),StopTime/ins_dt + 1);

% initialise storage vars to save time
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

GRASOn = 0;

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


x_hat_out_full = zeros(NumberStates,1);
x_hat_minus = x_hat_out_full;

PositionVariance = PositionError ^ 2;
VelocityVariance = VelocityError ^ 2;
ClockBiasVariance = ClockBiasError ^ 2;
ClockDriftVariance = ClockDriftError ^2;
%

sigma_pr_gnd = 5.0;
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



GyroBias = zeros(NumberEpochsGPS,3);
AccelBias = zeros(NumberEpochsGPS,3);



 
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




% two loops will run - a high speed loop calculating teh INS solution at
% 50 Hz - The data is collected at 100Hz, and a Runge-Kutta integration
% scheme is used 

Pfa = 2.22e-8;
SV_TrackTime = zeros(NumberSVs,1);

disp(strcat(ModuleName,sprintf('Beginning navigation loop for %d Epochs',NumberEpochsINS)));

%% check that start and stop times match the data


%%
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
    

    % implement the inertial navigation solution here
    if Epoch_hi == 1
%        % first time step

        q_INS(Epoch_hi,:) = EulerToQuat(InitialAttitude);
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
        F_INS(1:9,1:9) = GARD_GenerateINSFMatrix(Pos_LLH(Epoch_hi,:),Vel_NED(Epoch_hi,:),Acc_NED(Epoch_hi,:)-[0 0 9.80]);

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

        Q(10,10) = 2*x_accel_beta*x_accel_Q;
        Q(11,11) = 2*y_accel_beta*y_accel_Q;
        Q(12,12) = 2*z_accel_beta*z_accel_Q;
        Q(13,13) = 2*x_gyro_beta*x_gyro_Q;
        Q(14,14) = 2*y_gyro_beta*y_gyro_Q;
        Q(15,15) = 2*z_gyro_beta*z_gyro_Q;

%         Q(16,16) = (Sc)*ins_dt + Sf*(ins_dt^3)/3;
%         Q(17,17) = Sf*ins_dt;
%         Q(16,17) = Sf*(ins_dt^2)/2;
%         Q(17,16) = Sf*(ins_dt^2)/2;
        
        Q(16,16) = (Sc) + Sf*(ins_dt^2)/2;
        Q(17,17) = Sf;
        Q(16,17) = Sf*(ins_dt);
        Q(17,16) = Sf*(ins_dt);


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
        
        % get sensor measurements
          if(Epoch_lo < 2)
             omega_x = imudata(imudataindex,8);
             omega_y = imudata(imudataindex,9);
             omega_z = imudata(imudataindex,10);
 
             A_xb = imudata(imudataindex,5);
             A_yb = imudata(imudataindex,6);
             A_zb = imudata(imudataindex,7);
         else
            omega_x = imudata(imudataindex,8) - GyroBias(Epoch_lo,1);
            omega_y = imudata(imudataindex,9) - GyroBias(Epoch_lo,2);
            omega_z = imudata(imudataindex,10) - GyroBias(Epoch_lo,3);

            A_xb = imudata(imudataindex,5) - AccelBias(Epoch_lo,1);
            A_yb = imudata(imudataindex,6) - AccelBias(Epoch_lo,2);
            A_zb = imudata(imudataindex,7) - AccelBias(Epoch_lo,3);
          end
         
          t_INS(Epoch_hi,:) = imudata(imudataindex,1);
          
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
        C_BN = out.C_BN;
        q_INS(Epoch_hi,:) = out.q_INS;
        
        P_out_diag(:,Epoch_hi) = diag(out.P);
        %% calculate error bounds
        Pos_var(Epoch_hi,:) = P_out_diag(1:3,Epoch_hi).*[RM^2;RP^2*cos(Pos_LLH(Epoch_hi,1))^2;1];
        
        x_hat_minus = out.x;
        P_in_full = out.P;
        


    end % 
    
    
    
        % low-speed loop - kalman filter
        if(mod(Epoch_hi,IMU_RATE/GPS_RATE) == 0)
        %if(0)  %% uncomment this line to turn off the KF update
            Epoch_lo = Epoch_lo + 1;
        
        
            
            % find the gps data index of this epoch
            gpsdataindex = find(gpstime(:,3) == gpstime(60,3) + StartTime + (Epoch_lo-1)*GPS_RATE,1);
            
            % update the current imu data index
            imudataindex_i = find(imudata(:,1) < gpstime(gpsdataindex,1), 1,'last');
            
           
            
           % if(Epoch_lo == StartTime + 60)  % enable fault detection
                FDEnabled = 1;
            %end

            if Epoch_lo == 1
                continue;
            end
            
            clear PR_Vec PRR_Vec H R K V z_Vec NumberGPSMeasurements PR_Vec_minus PRR_Vec_minus Tecef2ned Tned2ecef SV_Vec SVPos SVVel;

            
            % calculate PR error
            %if Epoch_lo > 150
            %    PR_error = (Epoch_lo - 150)*3.0; % add 0.5 m/s ramp fault from time 150
            %else
                PR_error = 0.0;
            %end
            
            % generate current rotation matrix from earth to navigation frame
            % coordinates
            Tecef2ned= T_ECEF2NED(Pos_LLH(Epoch_hi,1),Pos_LLH(Epoch_hi,2));
            Tned2ecef = Tecef2ned';
        
            % get the gps time - used for the iono correction
            %GPSTime  = gps.RangeData(Epoch_lo).GPSSec;
            

            % UserPos and UserVel are the current state vector
            % which the EKF will linearise on
            UserPos(1:3) =  LLH2ECEF(Pos_LLH(Epoch_hi,1),Pos_LLH(Epoch_hi,2),Pos_LLH(Epoch_hi,3));
            
            % adjust for lever arm
            UserPos(1:3) = UserPos(1:3) + (Tned2ecef * C_BN * GPSLeverArm)';
            
            UserPos(4) = UserClock(Epoch_hi,1);
            
            UserVel(1:3) = Tned2ecef * Vel_NED(Epoch_hi,:)';
            UserVel(1:3) = UserVel(1:3) + (Tned2ecef * C_BN * cross(Omega_in',GPSLeverArm))'; 
            
            UserVel(4) = UserClock(Epoch_hi,2);
            
            
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


                        if SVIndex == 1
                           PR_Vec(SVIndex) = PR_Vec(SVIndex) + PR_error; 
                        end
                        

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
                            TropoDelay(Epoch_lo,SV) =  GARD_TropoDelay(SV_Elevation(Epoch_lo,SV),Pos_LLH(Epoch_hi,3));
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
                z_Vec = ([PR_Vec_raw(1:NumberGPSMeasurements) PRR_Vec_raw(1:NumberGPSMeasurements)]' - [PR_Vec_minus PRR_Vec_minus]');                
                
                
                H(1:NumberGPSMeasurements,1:3) = H(1:NumberGPSMeasurements,1:3)*Tecef2ned';
                H(NumberGPSMeasurements+1:2*NumberGPSMeasurements,4:6) = H(NumberGPSMeasurements+1:2*NumberGPSMeasurements,4:6)*Tecef2ned';
                
%                 %% add compass measurement
%                 Heading_error = CompassSensor(Epoch_hi) - psi_q(Epoch_hi);
%                 if(Heading_error > pi)
%                     Heading_error = Heading_error - 2*pi;
%                 end
%                 if(Heading_error < -pi)
%                     Heading_error = Heading_error + 2*pi;
%                 end
                
                %HeadingIndex = size(H,1)+1;
                %H(HeadingIndex,9) = 0;%-1;
                %z_Vec(HeadingIndex) = 0;% Heading_error;
                %R(HeadingIndex,HeadingIndex) = 100;%CompassVariance*10;
                
                


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

                
                C_BN = (eye(3,3) - del_att_skew) * C_BN; % DCM correction
                
                
                C_BN = GARD_OrthogonaliseDCM(C_BN);


                q_INS(Epoch_hi,:) = GARD_DCMToQuat(C_BN);
                
                % calculate euler angles for reference
                phi_q(Epoch_hi) = atan2(C_BN(3,2),C_BN(3,3));
                theta_q(Epoch_hi) = asin(-C_BN(3,1));
                psi_q(Epoch_hi) = atan2(C_BN(2,1),C_BN(1,1));
                
                
                
                x_out = [Pos_LLH(Epoch_hi,:)';Vel_NED(Epoch_hi,:)';q_INS(Epoch_hi,:)'];
                
                
                %%% END FULL SOLUTION
                
                if FDEnabled == 1
                
                %% start sub solutions
                NumberSubMeasurements = NumberGPSMeasurements-1;
                NumberSubFilters = nchoosek(NumberGPSMeasurements,NumberSubMeasurements);

                R_Sub = eye(2*NumberSubMeasurements);
                for k = 1:NumberSubMeasurements
                    R_Sub(k,k) = RangeNoiseVariance;
                    R_Sub(k+NumberGPSMeasurements,k+NumberGPSMeasurements) = RangeRateNoiseVariance;
                end
                
                % initialise sub filters
                if ~exist('P_minus_Sub','var')
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


                    % zero the row of the H matrix corresponding to the
                    % excluded satellite
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

                    EKF_lambda_ss_H(Epoch_lo,j) = (EKF_Beta_ss_H(:,j)' * pinv(EKF_B_ss_H(:,:,j)) * EKF_Beta_ss_H(:,j));

                    EKF_B_lambda_H = abs(eigs(EKF_B_ss_H(:,:,j)));
                    %EKF_TD_H(Epoch_lo,j) = sqrt(max(EKF_B_lambda_H)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));
                    EKF_TD_H(Epoch_lo,j) = chi2inv(1-PFalseAlarm/NumberGPSMeasurements,NumberGPSMeasurements-4);
                    
                    if(sqrt(EKF_lambda_ss_H(Epoch_lo,j)) > EKF_TD_H(Epoch_lo,j))
                        disp(sprintf('%s H-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
                    end

                    EKF_Beta_ss_V(:,j) = NED_ss(3);
                    EKF_B_ss_V(:,:,j) = abs(P_ned_Sub(3,3) - P_ned(3,3));

                    EKF_lambda_ss_V(Epoch_lo,j) = (EKF_Beta_ss_V(:,j)' * pinv(EKF_B_ss_V(:,:,j)) * EKF_Beta_ss_V(:,j));

                    EKF_B_lambda_V = eigs(EKF_B_ss_V(:,:,j));
                    %EKF_TD_V(Epoch_lo,j) = sqrt(max(EKF_B_lambda_V)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));
                    EKF_TD_V(Epoch_lo,j) = chi2inv(1-PFalseAlarm/NumberGPSMeasurements,NumberGPSMeasurements-4);
                    
                    if(sqrt(EKF_lambda_ss_V(Epoch_lo,j)) > EKF_TD_V(Epoch_lo,j))
                        disp(sprintf('%s V-Fault detected at Epoch %d on Sub-filter %d',ModuleName,Epoch_lo,j));
                    end

                    %% this HPL represents the fault-free (H0) hypothesis. 
%                     Sigma_MD_H = abs(norminv(P_MD_H,0,1));
%                     Sigma_MD_V = abs(norminv(P_MD_V,0,1));
%                     EKF_HPL_sub(Epoch_lo,j) = Sigma_MD_H * sqrt(P_ned(1,1) + P_ned(2,2)) + EKF_TD_H(Epoch_lo,j);
%                     EKF_VPL_sub(Epoch_lo,j) = Sigma_MD_V * sqrt(P_ned(3,3)) + EKF_TD_V(Epoch_lo,j);
%                     


                        EKF_P_lambda_H = abs(eigs(P_ned(1:2,1:2)));
                        EKF_P_lambda_V = abs(eigs(P_ned(3,3)));

                        EKF_K_H0_H = CEP_TableLookup_PmdH(min(EKF_P_lambda_H) / max(EKF_P_lambda_H));
                        EKF_K_H0_V = CEP_TableLookup_Pmd_V(1); % only have one eigen value for vert... 

                        EKF_HPL_H0_sub(Epoch_lo,j) = EKF_K_H0_H * sqrt(max(EKF_P_lambda_H));
                        EKF_VPL_H0_sub(Epoch_lo,j) = EKF_K_H0_V * sqrt(max(EKF_P_lambda_V));

                        %% todo - calculate the Fault-in progress (h1) hypothesis based
                        %% HPL_H1.

                        % HPL_H1 = max (HPE_B, HPE_NP, HPE_NB)
                        Pbias_H = sqrt(RAIM_lambda_H(NumberGPSMeasurements-4));
                        Pbias_V = sqrt(RAIM_lambda_V(NumberGPSMeasurements-4));

                        HPE_B = Pbias_H * sqrt(max(EKF_B_lambda_H));
                        eigs_sub = abs(eigs(P_ned_Sub(1:2,1:2)));
                        K_NP_H = CEP_TableLookup_PmdH(min(eigs_sub)/max(eigs_sub));
                        HPE_NP = K_NP_H * sqrt(max(eigs_sub));

                        K_NB_H = CEP_TableLookup_PmdH(min(EKF_B_lambda_H) / max(EKF_B_lambda_H));
                        HPE_NB = K_NB_H * sqrt(max(EKF_B_lambda_H));


                        EKF_HPE(j) = HPE_B + HPE_NP - HPE_NB;

                        VPE_B = Pbias_V * sqrt(max(EKF_B_lambda_V));

                        eigs_sub = abs(eigs(P_ned_Sub(3,3)));
                        K_NP_V = CEP_TableLookup_Pmd_V(1);
                        VPE_NP = K_NP_V * sqrt(max(eigs_sub));

                        K_NB_V = CEP_TableLookup_Pmd_V(1);
                        VPE_NB = K_NB_V * sqrt(max(EKF_B_lambda_V));

                        EKF_VPE(j) = VPE_B + VPE_NP - VPE_NB;


                    
                    end

            
                    clear z_Vec_Sub EKF_H_Sub R_Sub K_Sub V_Sub ys_kminus_Sub EKF_y_k_Sub PR_Vec_minus_Sub PRR_Vec_minus_Sub PR_SubVec PRR_SubVec; 
                
%             EKF_HPL(Epoch_lo) = max(EKF_HPL_sub(Epoch_lo,:));
%             EKF_VPL(Epoch_lo) = max(EKF_VPL_sub(Epoch_lo,:));
                    EKF_HPL_H1(Epoch_lo) = max(EKF_HPE);
                    EKF_VPL_H1(Epoch_lo) = max(EKF_VPE);

                    EKF_HPL_H0(Epoch_lo) = max(EKF_HPL_H0_sub(Epoch_lo,:));
                    EKF_VPL_H0(Epoch_lo) = max(EKF_VPL_H0_sub(Epoch_lo,:));

                    EKF_HPL(Epoch_lo) = max(EKF_HPL_H0(Epoch_lo),EKF_HPL_H1(Epoch_lo));
                    EKF_VPL(Epoch_lo) = max(EKF_VPL_H0(Epoch_lo),EKF_VPL_H1(Epoch_lo));
                
                end % if FDEnabled
        
        %end % this is the end of the if statement which decides if this is the first epoch or not
        
        

        if mod(Epoch_lo,10) == 0
           disp(sprintf('Completed Epoch %d',Epoch_lo)); 
        end
        
    end % end low-speed loop
    
    
    
end % end high-speed loop




        
disp(strcat(ModuleName,'Done!'));



results.Pos_LLH = Pos_LLH;
results.Vel_NED = Vel_NED;
results.Acc_NED = Acc_NED;

results.omega_b = omega_b;
results.A_b = A_b;

results.gvec_phi = gvec_phi;
results.gvec_theta = gvec_theta;

if exist('EKF_HPL','var')
    results.EKF_HPL = EKF_HPL;
end
if exist('EKF_VPL','var')
    results.EKF_VPL = EKF_VPL;
end

results.phi_q = phi_q;
results.theta_q = theta_q;
results.psi_q = psi_q;
results.UserClock = UserClock;
results.GyroBias = GyroBias;
results.AccelBias = AccelBias;
results.SatellitesUsed = SatellitesUsed;


results.imudatatime = imudatatime;

results.P_out = P_out_diag;
results.Pos_var = Pos_var;

results.EKF_lambda_ss_H = EKF_lambda_ss_H;
results.EKF_lambda_ss_V = EKF_lambda_ss_V;

results.EKF_TD_H = EKF_TD_H;
results.EKF_TD_V = EKF_TD_V;

results.SV_Vec = SV_Vec_save;



