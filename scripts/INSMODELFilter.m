clear all
close all
warning off


%purpose of this m file is to show if and how much the ins accel and pqr
%data can be filtered by the aeromodel, since the ins is corrupted by
%engine noise etc, temp, vibration, but aeromodel predicts pure
%aerodynamics

%This script runs INS and MODEL and GPS
%but runs a KF first to filter the INS accels and gyros with aeromodel measurements before going into
%the filtered INS and GPS filter.

%
% % 28.2.07
% % Q values for the IMU come from gyro specs, and are the same for full and
% % sub filter solution.
% %thefeore Q values for the model will come from the pqr and related to
% %uncertcomainty in parameter estimates.
% % R values come from predicted measurement noise eg 7.5^2 and should be the same for full and sub filter solutions% if test statistics are large, there is likely a problem with sub filter
% % processing somewhere
% %
% %
% % %$Id: INSMODELFilter.m 1884 2008-07-15 05:54:33Z n2523710 $
% % %Troy Bruggemann  11 August 2005
GPSConstants; %All generic constants put in this script, as global variables;
SensorNoiseParameters;
% %
% % %This approach uses loosely coupled approach, is a much simpler
% % %implemetnation than tightly coupled.
% %
feature accel on
% %
% % % load 'AeromodelWorkspace\P_saveModel.mat';
% % % load 'AeromodelWorkspace\x_state_total_plusModel.mat';
% %
% % %load the data
LoadData;
%
%
% %With GPS need to assume that the GPS antenna is at the centre of gravity of the airplane. Incorporate GPS antenna location from cofg in body axes later.
% %load the truth data
%
%
% %startepoch = 10*Rate;  %start at 10 seconds
% %101 should be 2 seconds into the simulation, ie 1 +2*50 = 101.
startepochHighRate = 201;  %start epoch is 10 for 1 Hz Rate and 19 for 2 Hz
% %data is at 50 Hz so take every 50th for 1 second
%
% %startepochHighRate = 951;
%endepochHighRate = 12100;
endepochHighRate = 6050;
%
%
% %startepochHighRate = 6001;  %start epoch is 10 for 1 Hz Rate and 19 for 2 Hz
% %data is at 50 Hz so take every 50th for 1 second
%
% %endepochHighRate = 10001;
%
%
% %endepochHighRate = 20190;
% %endepochHighRate = 40000;
%
% %endepochHighRate = 25960;
%
% %this is the epochs at 1 Hz
startepoch = 3;  %need to start this at 2 seconds, which is 101 ie. i*50+1 in 50 Hz
%endepoch = 120;
endepoch = 60;
% %endepoch = 780;
%
%
% %endepoch = 400;
% %startepoch = 123; %want it to start at 120 seconds but need to start it 3 second later   %need to start this at 2 seconds, which is 101 ie. i*50+1 in 50 Hz
% %endepoch = 200;
%
% %endepoch = 350;
% %endepoch = 750;
%
% %reads in data, formatting
%
ReadData;
GenerateNoise;

%get data from air data sensors
AirDataSensors;


%initialise with zeros for GPSsensorSub


Xpos_GPSSub = zeros(32,endepoch);
Ypos_GPSSub = zeros(32,endepoch);
Zpos_GPSSub = zeros(32,endepoch);
dt_GPSSub = zeros(32,endepoch);


%get data from GPS
GPSSensorSub;


INSSensor;
%
% %===========================================
% %WIND ESTIMATION
% %===========================================
%
% WindEstimateFromSensors;



%
%   startepochHighRate = 951;
%   endepochHighRate = 6050;
%
%    %endepochHighRate = 20190;
%
%    startepoch = 20;
%    endepoch = 120;
%    endepoch = 400;
%
%
%
% endepochHighRate = 6052;
% %
% %
% % %startepochHighRate = 6001;  %start epoch is 10 for 1 Hz Rate and 19 for 2 Hz
% % %data is at 50 Hz so take every 50th for 1 second
% %
% endepochHighRate = 7010;
% %
% %
% % %endepochHighRate = 20000;
% % %endepochHighRate = 40000;
% %
% % endepochHighRate = 25960;
%
% %this is the epochs at 1 Hz
% startepoch = 3;  %need to start this at 2 seconds, which is 101 ie. i*50+1 in 50 Hz
% endepoch = 120;
%startepoch = 20;

%endepoch = 140;
%startepoch = 123; %want it to start at 120 seconds but need to start it 3 second later   %need to start this at 2 seconds, which is 101 ie. i*50+1 in 50 Hz
%endepoch = 200;

%endepoch = 350;
%endepoch = 750;


%=======================================

%===============================================================
%estimate WIND North, East Wind from Truth data
%==============================================================

%
% for i = 1:120
%
%     ss(i) = norm(Quaternions_truth1Hz(:,i));
%
% end
%
% plot(ss(2:120));

%===============================================================
%END WIND ESTIMATE
%==============================================================


%============================================================
%AERODYNAMIC MODEL
%============================================================

%=====================================================================
%Atmospheric Parameters
%=====================================================================

%Standard Atmosphere in SI units

%Altitude % ASsume GA aircraft operating betweeen 0 and 15,000 feet.

%GravityTruth(i)



%get the bias and noise from the data.

%   Signal1 = INSSensorNoise1Hz(1,:);
%   [CurveFit1, BestFit] = LeastSquaresBestFit(Signal1, 10, 1);
%
%     Signal2 = INSSensorNoise1Hz(2,:);
%   [CurveFit2, BestFit] = LeastSquaresBestFit(Signal2, 10, 1);
%
%     Signal3 = INSSensorNoise1Hz(3,:);
%   [CurveFit3, BestFit] = LeastSquaresBestFit(Signal3, 10, 1);
%
%     Signal4 = INSSensorNoise1Hz(4,:);
%   [CurveFit4, BestFit] = LeastSquaresBestFit(Signal4, 10, 1);
%
%     Signal5 = INSSensorNoise1Hz(5,:);
%   [CurveFit5, BestFit] = LeastSquaresBestFit(Signal5, 10, 1);
%
%     Signal6 = INSSensorNoise1Hz(6,:);
%   [CurveFit6, BestFit] = LeastSquaresBestFit(Signal6, 10, 1);
%
%
%
%   TruthINSBias(1,:) = CurveFit1;
%    TruthINSBias(2,:) = CurveFit2;
%     TruthINSBias(3,:) = CurveFit3;
%      TruthINSBias(4,:) = CurveFit4;
%       TruthINSBias(5,:) = CurveFit5;
%        TruthINSBias(6,:) = CurveFit6;
%
% %
% TruthINSNoise(1,:) = Signal1-CurveFit1;
% TruthINSNoise(2,:) = Signal2-CurveFit2;
% TruthINSNoise(3,:) = Signal3-CurveFit3;
% TruthINSNoise(4,:) = Signal4-CurveFit4;
% TruthINSNoise(5,:) = Signal5-CurveFit5;
% TruthINSNoise(6,:) = Signal6-CurveFit6;

% TruthINSBias(1,:) = 0;
% TruthINSBias(2,:) = 0;
% TruthINSBias(3,:) = 0;
% TruthINSBias(4,:) = 0;
% TruthINSBias(5,:) = 0;
% TruthINSBias(6,:) = 0;
%
%
% TruthINSNoise(1,:) = 0;
% TruthINSNoise(2,:) = 0;
% TruthINSNoise(3,:) = 0;
% TruthINSNoise(4,:) = 0;
% TruthINSNoise(5,:) = 0;
% TruthINSNoise(6,:) = 0;

%
% jjj = xcorr(TruthINSNoise(1,:),TruthINSNoise(1,:));
%
%
% jjj2 = xcorr(TruthINSBias(1,:),TruthINSBias(1,:));
% jjj3 = xcorr(Signal1,Signal1);





%start point

%run the integration processor



gps_dt = 1;
model_dt = 0.01;
KF_dt = 1;  %update rate for KF

%estimate GPS velocities

GRASOn = 0;

% setup error values
if(GRASOn == 1)
    PositionErrorGPS = 1.0;
    PositionErrorModel = 20.0;
    VelocityErrorGPS = 2.0;
    VelocityErrorModel = 4.0;
    ClockBiasError = 1E9;   %this is in metres
    ClockDriftError = 100.0;  %this is in metres
    RangeNoiseVariance = 3.0^2; % m
    RangeRateNoiseVariance = 2.0^2; % m/s
else
    PositionErrorGPS = 5.0;
    PositionErrorModel = 20.0;
    VelocityErrorGPS = 2.0;
    VelocityErrorModel = 40.0;
    ClockBiasError = 3;
    ClockDriftError = 1.0;
    RangeNoiseVariance = 10^2; % m
    RangeRateNoiseVariance = 5.0^2; % m/s

end




% %add error to INS
% for i = startepoch:endepoch
%
% ax_b_INS1Hz(i) = ax_b_INS1Hz(i) ;
% ay_b_INS1Hz(i) = ay_b_INS1Hz(i);
% az_b_INS1Hz(i) =  az_b_INS1Hz(i);
%
%
% p_INS1Hz(i) = p_INS1Hz(i);
% q_INS1Hz(i) = q_INS1Hz(i);
% r_INS1Hz(i) = r_INS1Hz(i);
%
%
% end


iprev = startepoch - 5;

%initialize
P_inSubTotal = zeros(32,17,17,endepoch);

for i = startepoch:endepoch



    N = N_save(i);

    %clear variables which are re-used each epoch and have varying vector
    %or matrix sizes depending on the number of measurements

    clear SatPosEKF
    clear SatVelEKF
    clear PseudorangesEKF
    clear PseudorangeRatesEKF
    clear H
    clear R
    clear r_VecCalc
    clear PR_Vec_minus
    clear r_VecCalcVel
    clear Relative_velocity
    clear PRR_Vec_minus
    clear z
    clear V
    clear K
    clear delta_pr_omegaedot
    clear z_INS
    clear z_M
    clear z_INSandMODEL

    clear PRCalcINS
    clear PR_rate_CalcINS


    clear PRCalcMODEL
    clear PR_rate_CalcMODEL

    clear PRCalcMODELandINS
    clear PR_rate_CalcMODELandINS




    %     NumberMeasurements = 6 + (6 + 6 + 3) ; %This is 6 INS measurements, 6 model measurements %appended to 6 more for model accel and gyro and 3 wind in NED
    %     NumberStates = (9 + 6)+ (9 + 6 +3); %total 33 states. the extra 3 are 3 wind components in north east down velocity

    %NumberMeasurements = 6 + (6 + 6 + 3) ; %This is 6 INS measurements, 6 model measurements %appended to 6 more for model accel and gyro and 3 wind in NED

    NumberMeasurements = 2*N_save(i);
    %NumberStates = (9 + 6+ 9 +6 +  9 +6 + 2+2 + 2 ); %2+ 2 is the clock bias and drift from the INS, and model %total 51 states. the extra 3 are 3 wind components in north east down velocity

    NumberStates = (9 +2+6);

    %just use a constant value, this is used to convert from metres to
    %radians and vice versa
    ReConst = 6348545.90022133;

    %run at first time - initialisation

    if i == startepoch


        %
        %         P_in_INS = zeros(9,9);
        %  P_in_M = zeros(9,9);




        P_in_dash = [0.052^2,   %INS att
            0.052^2,
            0.052^2,
            20^2,  %INS vel NED
            20^2,
            20^2,
            (50/ReConst)^2,  %lat error rad
            (50/ReConst)^2,  %lon error rad
            (50)^2,  %down error m
            8^2,  %GPS clock bias error from ins model and gps m/s
            0.4^2, %GPS clock drift error m/s
            0.001,  %x accel bias INS
            0.001,   %y accel bias INS
            0.001,   %z accel bias  INS
            0.001,  %x gyro bias INS
            0.001,   %y gyro bias INS
            0.001 ]*100;   %z gyro bias INS and model ];



        for index = 1:NumberStates
            P_in(index,index) = P_in_dash(index);
        end

        P_in = P_in*10;

        %P_in_M = P_in_INS;

        %initial starting

        %assume the INS is calibrated against the GPS

        %I can replace Roll pitch yaw euler angles here with quaternions for the
        %attitude representation.

        %estimate velocities from truth

        TMatrix_ECEF2NED = T_ECEF2NED( Lon_truth1Hz(i), Lat_truth1Hz(i));

        VelocityNED = TMatrix_ECEF2NED*[Xvel_truth1Hz(i),Yvel_truth1Hz(i),Zvel_truth1Hz(i)]';

        %use the truth as starting point.

        %         omega_xMODEL1Hz(i) = p_INS1Hz(i);
        %         omega_yMODEL1Hz(i) = q_INS1Hz(i);
        %         omega_zMODEL1Hz(i) = r_INS1Hz(i);
        %
        %         ax_b_MODEL1Hz(i) = ax_b_INS1Hz(i);
        %         ay_b_MODEL1Hz(i) = ay_b_INS1Hz(i);
        %         az_b_MODEL1Hz(i) = az_b_INS1Hz(i);


        omega_xMODEL1Hz(i) = 0;
        omega_yMODEL1Hz(i) = 0;
        omega_zMODEL1Hz(i) = 0;

        ax_b_MODEL1Hz(i) = 0;
        ay_b_MODEL1Hz(i) = 0;
        az_b_MODEL1Hz(i) = 0;


        omega_xMODEL1HzandINS(i) = 0;
        omega_yMODEL1HzandINS(i) = 0;
        omega_zMODEL1HzandINS(i) = 0;

        ax_b_MODEL1HzandINS(i) = 0;
        ay_b_MODEL1HzandINS(i) = 0;
        az_b_MODEL1HzandINS(i) = 0;


        %minus indicates apriori

        x_state_total_minusINS(:,i) = [Quaternions_truth1Hz(1,i),   %q0
            Quaternions_truth1Hz(2,i),    %q1
            Quaternions_truth1Hz(3,i),     %q2
            Quaternions_truth1Hz(4,i),     %q3
            V_n_truth1Hz(i),    %Vn
            V_e_truth1Hz(i),    %Ve
            V_d_truth1Hz(i),   %Vd
            Lat_truth1Hz(i),  %lat
            Lon_truth1Hz(i), %lon
            Hgt_truth1Hz(i)];


        x_state_total_GPSINS(:,i) = [0,0]; %clock bias, clock drift


        % calibrate initial gyro bias
        GyroBias(1) = mean(p_INS1Hz(1:30));
        GyroBias(2) = mean(q_INS1Hz(1:30));
        GyroBias(3) = mean(r_INS1Hz(1:30))


        x_state_total_INSBIAS(:,i) = [ GyroBias(1), GyroBias(2), GyroBias(3),0,0,0];
        x_state_total_plusINS(:,i) = x_state_total_minusINS(:,i);


        %for filtering INS with model


        %initial P_A
        P_in_A(1,1) = 10;
        P_in_A(2,2) = 10;
        P_in_A(3,3) = 10;
        P_in_A(4,4) = 10;
        P_in_A(5,5) = 10;
        P_in_A(6,6) = 10;


    else

        %use truth for clock state
        x_state_total_GPSINS(1,i) = PosTruthGPS(i,4);
        x_state_total_GPSINS(2,i) = VelTruthGPS(i,4);

        %x_state_total_INSBIAS(:,i) = [  TruthINSBias(4,i) , TruthINSBias(5,i) , TruthINSBias(6,i) , TruthINSBias(1,i) , TruthINSBias(2,i) , TruthINSBias(3,i) ];


        x_state_total_INSBIAS(:,i) = [ 0, 0, 0,0,0,0];



        %=======================================================
        %INS PROPAGATION
        %=======================================================
        %this is run at 50 Hz


        X_state_inINS(1:10,i) = x_state_total_plusINS(1:10,i-1);


        X_state_outTempINS_F(:,1) = X_state_inINS(1:10,i);

        X_state_outTempINS(:,1) = X_state_inINS(1:10,i);


        P_minusOld = P_in;

        for pp = 1:100
            
            
            m = (i-1)*100+pp;

            omega_x_INS = p_INS_50Hz((i-1)*100+pp) - x_state_total_INSBIAS(1,i);
            omega_y_INS = q_INS_50Hz((i-1)*100+pp)- x_state_total_INSBIAS(2,i);
            omega_z_INS = r_INS_50Hz((i-1)*100+pp)- x_state_total_INSBIAS(3,i);

            A_xb_INS = ax_b_INS_50Hz((i-1)*100+pp)- x_state_total_INSBIAS(4,i);
            A_yb_INS = ay_b_INS_50Hz((i-1)*100+pp)- x_state_total_INSBIAS(5,i);
            A_zb_INS = az_b_INS_50Hz((i-1)*100+pp)- x_state_total_INSBIAS(6,i);

            A_b_INS = [A_xb_INS, A_yb_INS, A_zb_INS];
            omega_b_INS = [omega_x_INS, omega_y_INS, omega_z_INS];

            dt_INS = 0.01;

            llh_dot = [0,0,0];% (this is unused) in the function
            %gravity = -9.80;

            %Get Gravity estimate (use truth at the moment)
            g_INS = GravityTruth((i-1)*100+pp); %this is at 50 Hz
            %g_INS = 0;
            %g_INS = Earth_Gravity(X_state_outTempINS(8:10,pp)',g_para1,g_para2,g_para3); %uses lat lon and hgt to estimate g

            [X_state_outTempINS_F(:,pp+1)] = INS_Mechanization2(X_state_outTempINS_F(:,pp), A_b_INS, omega_b_INS, g_INS, dt_INS,llh_dot);


            %convert quaternion to euler
            quat_M = [ X_state_outTempINS_F(1,pp), X_state_outTempINS_F(2,pp),X_state_outTempINS_F(3,pp),X_state_outTempINS_F(4,pp)];
            euler_M = QuatToEuler(quat_M);


            TMatrixBtoNED_M = T_Body2NED(euler_M(1),euler_M(2), euler_M(3));

            %convert NED velocities to body
            Vel_Body_M = TMatrixBtoNED_M'*[X_state_outTempINS_F(5,pp),X_state_outTempINS_F(6,pp),X_state_outTempINS_F(7,pp)]';


            u_b_M = Vel_Body_M(1);
            v_b_M = Vel_Body_M(2);
            w_b_M = Vel_Body_M(3);


            %         V_T0  = sqrt(u_b_M^2 + v_b_M^2 + w_b_M^2);  %this should be airspeed not body, is approximate only, source of error?
            %         beta0  = asin(v_b_M/V_T0);
            %         alpha0 = atan2(w_b_M,u_b_M);



            V_T0  = Airspeed_truth1Hz(i);
            beta0  = Beta_truth1Hz(i);
            alpha0 = Alpha_truth1Hz(i);


            %         Airspeed_truth1Hz(i) =  Airspeed_truth(i*50+1);
            %     Beta_truth1Hz(i) =  Beta_truth(i*50+1);
            %     Alpha_truth1Hz(i) =  Alpha_truth(i*50+1);
            %


            phi0 = euler_M(1);
            theta0 = euler_M(2);
            psi0 = euler_M(3);

            p0 = p_INS_50Hz(m-1);
            q0 = q_INS_50Hz(m-1);
            r0 = r_INS_50Hz(m-1);



            pn0 = 0;
            pe0 = 0;
            h0 = 0;  %note this is up not down
            q00 =  X_state_outTempINS_F(1,pp);
            q10 =  X_state_outTempINS_F(2,pp);
            q20 =  X_state_outTempINS_F(3,pp);
            q30 =  X_state_outTempINS_F(4,pp);
            Lat0 = X_state_outTempINS_F(8,pp);
            Lon0 =  X_state_outTempINS_F(9,pp);
            Hgt0 =  X_state_outTempINS_F(10,pp);

            %Wind_NED = X_state_inMODEL(17:19,i)';

            Wind_NED = [0 0 0];

            %Parameter_Noise = [CL_noise(i), CD_noise(i),CY_noise(i),Cl_noise(i),Cm_noise(i),Cn_noise(i)];
            
            
             Parameter_Noise = [0, 0,0,0,0,0];

            %this is the start state
            x0IN(1,:) = [V_T0, beta0 ,alpha0, phi0, theta0, psi0, p0,q0,r0,pn0, pe0, h0, q00, q10, q20, q30, Lat0, Lon0, Hgt0,u_b_M,v_b_M,w_b_M];

            g_M = Earth_Gravity([Lat0,Lon0,Hgt0]',g_para1,g_para2,g_para3); %uses lat lon and hgt to estimate g


            Latprev = Lat0;
            Lonprev = Lon0;
            Hgtprev = Hgt0;


            %=========================================================
            %Now run the aeromodel prediction and EKF, to filter p_INS_50Hz
            %etc and then recalculate the X_state_outTempINS for use with
            %the GPS


           % m = (i-1)*100+pp;  %note this is an index, not mass

            %[X_state_outTemp(:,pp+1)] = INS_Mechanization2(X_state_outTemp(:,pp), A_b, omega_b, g, ins_dt,llh_dot);

            %delta_x_control(m,:) = [Elevator_truth(m), Aileron_truth(m), Rudder_truth(m), Throttle_truth(m)]; %this throttle needs to be scaled, not suree what it means at the moment, but probably needs to represent a change in thrust, since it is multiplied by the Xdelta_t term in the model    %control inputs will be the current control inputs at i

            %delta_x_control(m,:) = [0, 0, 0,0]; %this throttle needs to be scaled, not suree what it means at the moment, but probably needs to represent a change in thrust, since it is multiplied by the Xdelta_t term in the model    %control inputs will be the current control inputs at i

            %delta_x_control(m,:) = [Elevator_noise(i), Aileron_noise(i), Rudder_noise(i),Throttle_noise(i)]; %this throttle needs to be scaled, not suree what it means at the moment, but probably needs to represent a change in thrust, since it is multiplied by the Xdelta_t term in the model    %control inputs will be the current control inputs at i

            delta_x_control(m,:) = [Elevator_truth(m) + Elevator_noise(i), Aileron_truth(m)+ Aileron_noise(i), Rudder_truth(m) + Rudder_noise(i), Throttle_truth(m)+ Throttle_noise(i)]; %this throttle needs to be scaled, not suree what it means at the moment, but probably needs to represent a change in thrust, since it is multiplied by the Xdelta_t term in the model    %control inputs will be the current control inputs at i


            AeroCoeffsTruth = [CD_truth(m), CY_truth(m), CL_truth(m), Cl_truth(m), Cm_truth(m), Cn_truth(m)];
            EngCoeffsTruth = [MAP_truthEng(m), m_dotAir_truthEng(m),m_dotFuel_truthEng(m),BSFC_truthEng(m),P_truthEng(m), OMEGA_truthEng(m)] ;
            PropCoeffsTruth = [J_truth(m), CT_truth(m) ,CP_truth(m)];

            ForcesTruth = [Faero_truth(1,m),Faero_truth(2,m),Faero_truth(3,m),Fprop_truth(1,m),Fprop_truth(2,m),Fprop_truth(3,m)];

            MomentsTruth = [  Mcg_truth(1,m),Mcg_truth(2,m),Mcg_truth(3,m),Mprop_truth(1,m),Mprop_truth(2,m),Mprop_truth(3,m)];

            MassTruth = Mass_truth(m);

            InertiaTruth = [Jx_truth(m), Jy_truth(m), Jz_truth(m), Jxz_truth(m)];
            CGposTruth = [CGxpos_truth(m), CGypos_truth(m), CGzpos_truth(m)];
            MachTruth =  Mach_truth(m);

            %calculate MACH number for next time around
            %M0 = x0IN(pp,1)/AtmosTruth(4,m);

            M0 = MachTruth;

            Wind_NED = [u_wind_truth(i), v_wind_truth(i), w_wind_truth(i)];

            AeroCoeffCorrections = [0 0 0 0 0 0];

            %AeroCoeffCorrections = -x_state_total_MODELBIAS(:,i);


            %correct p q r before it goes in:

            % x0IN(pp,7:9) = x0IN(pp,7:9) - x_state_total_MODELBIAS(1:3,1);

            [dx, AeroCoeffs,Forces, Moments] = AerodynamicModelStevens(M0,g_M,AtmosTruth(:,m),x0IN(pp,:),delta_x_control(m,:), AeroCoeffsTruth, EngCoeffsTruth,PropCoeffsTruth,ForcesTruth, MomentsTruth, MassTruth,InertiaTruth,CGposTruth,MachTruth,Parameter_Noise, Wind_NED,AeroCoeffCorrections);

            dt_M = 0.01;

            %integrate to get u v w p q r

            x0IN(pp+1,:) = x0IN(pp,:) + dt_M*dx(1:22);


            p = x0IN(pp+1,7);
            q = x0IN(pp+1,8);
            r = x0IN(pp+1,9);


            omega_xMODEL100Hz(m) = p;
            omega_yMODEL100Hz(m) = q;
            omega_zMODEL100Hz(m) = r;

            ax_b_MODEL100Hz(m)= dx(20) ;
            ay_b_MODEL100Hz(m)= dx(21) ;
            az_b_MODEL100Hz(m)= dx(22) ;


            phi = x0IN(pp,4);
            theta = x0IN(pp,5);
            psi = x0IN(pp,6);

            phi_dot = p +tan(theta)*(q*sin(phi) + r*cos(phi));
            theta_dot = q*cos(phi) - r*sin(phi);
            psi_dot = (q*sin(phi) + r*cos(phi))/cos(theta);


            %integrate to get phi theta psi
            x0IN(pp+1,4)  =  x0IN(pp,4) + dt_M*phi_dot;
            x0IN(pp+1,5)  = x0IN(pp,5)+ dt_M*theta_dot;
            x0IN(pp+1,6)  =  x0IN(pp,6)+ dt_M*psi_dot;


            u = x0IN(pp+1,20);
            v = x0IN(pp+1,21);
            w = x0IN(pp+1,22);


            phi =  x0IN(pp+1,4);
            theta =  x0IN(pp+1,5);
            psi =  x0IN(pp+1,6);



            Cphi = cos(phi);
            Ctheta = cos(theta);
            Cpsi = cos(psi);

            Sphi = sin(phi);
            Stheta = sin(theta);
            Spsi = sin(psi);


            Wind_NED(1) = 0;
            Wind_NED(2) = 0;
            Wind_NED(3) = 0;


            pn_dot = u*Ctheta*Cpsi + v*(-Cphi*Spsi + Sphi*Stheta*Cpsi) + w*(Sphi*Spsi+Cphi*Stheta*Cpsi);% - Wind_NED(1) + u_wind; %note u_wind are in NED at the moment until i get the stevens textbook

            pe_dot = u*Ctheta*Spsi + v*(Cphi*Cpsi + Sphi*Stheta*Spsi) + w*(-Sphi*Cpsi + Cphi*Stheta*Spsi);% - Wind_NED(2) + v_wind;

            h_dot = u*Stheta - v*Sphi*Ctheta - w*Cphi*Ctheta;% + Wind_NED(3) + w_wind;  %note that this appears to be 'up' not 'down' as I first thought.


            %integrate to get p n ed
            Lat = Lat0;
            Hgt = Hgt0;

            a = 6378137.0;   % semi-major axis (metres)
            f = 1/298.2572; % flattening
            e2 = f * (2-f); % eccentricity squared
            e = sqrt(e2);   % first eccentricity

            Rmeridian = a * (1 - e2) / (sqrt(1 - e2 * sin(Lat)^2))^3;

            %find the normal radius of curvature
            Rnormal = a / sqrt(1 - e2 * sin(Lat)^2);

            Lat_dot = pn_dot/(Rmeridian + Hgt);
            Lon_dot = pe_dot/((Rnormal + Hgt)*cos(Lat));
            Hgt_dot = h_dot;

            %g_M = Earth_Gravity([x0IN(pp,17),x0IN(pp,18),x0IN(pp,19)]',g_para1,g_para2,g_para3); %uses lat lon and hgt to estimate g

            Lat = Latprev + dt_M*Lat_dot;
            Lon = Lonprev + dt_M*Lon_dot;
            Hgt = Hgtprev + dt_M*Hgt_dot;


            %update
            Latprev = Lat;
            Lonprev = Lon;
            Hgtprev = Hgt;



            dt_A = 1/100;
            %Run the KF for INS and model integration


            %states are  p q r ax ay azfrom INS
            
            
            

            F_A = zeros(6,6);

            
            
            PHI_A = expm(F_A*dt_A);

%             Q_A(1,1) = (1/300)*0.00552^2;
%             Q_A(2,2) = (1/300)*0.00552^2;
%             Q_A(3,3) = (1/300)*0.00552^2;
%             Q_A(4,4) = (1/300)*0.0124^2;
%             Q_A(5,5) = (1/300)*0.0124^2;
%             Q_A(6,6) = (1/300)*0.0124^2;
% 
% 
%             %this is the noise of the aeromodel measurements
%             R_A(1,1) = 10*(1/300)*0.00552^2;
%             R_A(2,2) = 10*(1/300)*0.00552^2;
%             R_A(3,3) = 10*(1/300)*0.00552^2;
%             R_A(4,4) = 10*(1/300)*0.0124^2;
%             R_A(5,5) = 10*(1/300)*0.0124^2;
%             R_A(6,6) = 10*(1/300)*0.0124^2;
%             
            



            
            W_A(1,1) = 0.00245^2;
            W_A(2,2) = 0.00245^2;
            W_A(3,3) = 0.00245^2;
            W_A(4,4) = 0.00043^2;
            W_A(5,5) = 0.00043^2;
            W_A(6,6) = 0.00043^2;
            
            
            
            Q_A = PHI_A*W_A*PHI_A'*dt_A;


            %this is the noise of the aeromodel measurements
            R_A(1,1) = 0.1^2;
            R_A(2,2) = 0.1^2;
            R_A(3,3) = 0.1^2;
            R_A(4,4) = 0.1^2;
            R_A(5,5) = 0.1^2;
            R_A(6,6) = 0.1^2;



            state_Aminus = [p_INS_50Hz(m), q_INS_50Hz(m),r_INS_50Hz(m),ax_b_INS_50Hz(m),ay_b_INS_50Hz(m),az_b_INS_50Hz(m)];


            %state_Aminus = state_A + state_A*dt_A;


            H_A = eye(6); %6 states, 6 measurements


            %measurements z are difference between apriori states and aeromodel

            zModel_A = [omega_xMODEL100Hz(m), omega_yMODEL100Hz(m),omega_zMODEL100Hz(m),ax_b_MODEL100Hz(m),ay_b_MODEL100Hz(m),az_b_MODEL100Hz(m)];
            zINS_A = [p_INS_50Hz(m), q_INS_50Hz(m),r_INS_50Hz(m),ax_b_INS_50Hz(m),ay_b_INS_50Hz(m),az_b_INS_50Hz(m)];

            %z_A =  zINS_A - zModel_A;
            

            %which way should this be? 
             z_A =  zModel_A- zINS_A;

            %for the INS, the acceleration dots (jerks) and p q r dots are
            %functions of the noise since the noise is driving as inputs. For the aerodynamic model they are
            %functions of the aircraft dynamics


            P_minus_A = PHI_A * P_in_A * PHI_A' + Q_A;


            V_A = H_A * P_minus_A * H_A' + R_A;
            K_A = P_minus_A * H_A' * inv(V_A);

            x_hat_out_A = K_A * (z_A');   %this is the error estimate in acceleration


            P_out_A = (eye(6) - K_A*H_A)*P_minus_A*(eye(6)-K_A*H_A)' + K_A*R_A*K_A';

            P_in_A = P_out_A; %update for next time around

            %     K_save_A(:,1:NumberMeasurements,i) = K;
            %     H_save(1:NumberMeasurements,:,i) = H;

            % save x_hat_out
            x_hat_save_A(:,m) = x_hat_out_A;
            %P_save(:,:,i) = P_out;
            P_minus_save_A(:,:,m) = P_minus_A;

            %===============================================
            %Apply INS corrections to INS total state
            %===============================================


            %correct state
            x_state_total_A = state_Aminus + x_hat_out_A';

            %now Run the INS propagation again for the INS GPS integration


            %these are the filtered INS measurements (hopefully)
            p_INS_50Hz_F(m) = x_state_total_A(1);
            q_INS_50Hz_F(m) = x_state_total_A(2);
            r_INS_50Hz_F(m) = x_state_total_A(3);
            ax_b_INS_50Hz_F(m) = x_state_total_A(4);
            ay_b_INS_50Hz_F(m) = x_state_total_A(5);
            az_b_INS_50Hz_F(m) = x_state_total_A(6);


            %Now repropagate the INS with the filtered measurements

            omega_x_INS = p_INS_50Hz_F((i-1)*100+pp) ;
            omega_y_INS = q_INS_50Hz_F((i-1)*100+pp);
            omega_z_INS = r_INS_50Hz_F((i-1)*100+pp);

            A_xb_INS = ax_b_INS_50Hz_F((i-1)*100+pp);
            A_yb_INS = ay_b_INS_50Hz_F((i-1)*100+pp);
            A_zb_INS = az_b_INS_50Hz_F((i-1)*100+pp);

            A_b_INS = [A_xb_INS, A_yb_INS, A_zb_INS];
            omega_b_INS = [omega_x_INS, omega_y_INS, omega_z_INS];

            dt_INS = 0.01;

            llh_dot = [0,0,0];% (this is unused) in the function
            %gravity = -9.80;

            %Get Gravity estimate (use truth at the moment)
            g_INS = GravityTruth((i-1)*100+pp); %this is at 50 Hz
            %g_INS = 0;
            %g_INS = Earth_Gravity(X_state_outTempINS(8:10,pp)',g_para1,g_para2,g_para3); %uses lat lon and hgt to estimate g

            [X_state_outTempINS(:,pp+1)] = INS_Mechanization2(X_state_outTempINS(:,pp), A_b_INS, omega_b_INS, g_INS, dt_INS,llh_dot);



        end
        
        
         x_state_total_plusINS(1:10,i) = X_state_outTempINS(:,101);



    end
end



%plot original INS accelerations

%plot(ax_b_INS_50Hz_F);

%plot INS accelerations


figure;
 plot(p_INS_50Hz_F(300:1800));
 hold;

 plot(p_INS_50Hz(300:1800),'r');
 plot(GyroTruth100Hz(1,300:1800),'g');
 
pause;


figure;
 plot(p_INS_50Hz_F(300:1800)-GyroTruth100Hz(1,300:1800));
 hold;

 plot(p_INS_50Hz(300:1800)-GyroTruth100Hz(1,300:1800),'r');

 pause;

figure;
 plot(p_INS_50Hz_F(300:1800));
 hold;

 plot(p_INS_50Hz(300:1800),'r');
 plot(GyroTruth100Hz(1,300:1800),'g');
 
pause;


figure;
 plot(p_INS_50Hz_F(300:1800)-GyroTruth100Hz(1,300:1800));
 hold;

 plot(p_INS_50Hz(300:1800)-GyroTruth100Hz(1,300:1800),'r');

 
 
 Errorfiltered = p_INS_50Hz_F(300:1800)-GyroTruth100Hz(1,300:1800);
 
 Errornonfiltered = p_INS_50Hz(300:1800)-GyroTruth100Hz(1,300:1800);
 

 
 %now see if coasting performance of model and ins is better than ins alone




