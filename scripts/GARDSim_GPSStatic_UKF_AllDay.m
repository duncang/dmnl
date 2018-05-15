%% GPS Position Solution using UKF
% Written by Duncan Greer 2 Feb 2007
%
% $Id: GARDSim_GPSStatic_UKF_AllDay.m 1879 2008-07-15 05:20:21Z n2523710 $
%
%
%
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
%

%% setup simulation

% set number formatting for display
format long g;

% load GPS constants
GPSConstants;

%% Load data

load('data\rinex\Ashtech_CGRS_EESE_21Feb06.mat');

% find the number of rinex sets 
NumberSets = size(Feb21Data.ObsData,2);

% position
[InitialPosition(1),InitialPosition(2),InitialPosition(3)] = ECEF2LLH(Feb21Data.ObsData(1).ApproxPos);
% velocity
InitialVelocity = [0 0 0]';  % NED velocities in m/s
    
% roof truth
meanECEF = [-5046773.060  2568452.907 -2925282.788];

ApproxPos = Feb21Data.ObsData(1).ApproxPos;
[ApproxPos_LLH(1) ApproxPos_LLH(2) ApproxPos_LLH(3)] = ECEF2LLH(ApproxPos);
    
Tecef2ned = T_ECEF2NED(ApproxPos_LLH(2),ApproxPos_LLH(1));

% 3 position, 3 velocity and 2 clock
NumberStates = 8;

% process noise states - 8
ProcessNoiseStates = 8;  

% use a maximum of 6 measurements (PR+PRR)
MeasurementNoiseStates=12;
    
%% initialise covariances
Px_kminus = eye(NumberStates,NumberStates);
Px_kminus(1:3,1:3) = eye(3,3)*100^2;
Px_kminus(4:6,4:6) = eye(3,3)*10^2;
Px_kminus(7,7) = 100^2;
Px_kminus(8,8) = 10^2;

x_hat_kminus = zeros(NumberStates,1);
x_hat_kminus(1:3,1) = ApproxPos;
x_hat_kminus(7,1) = 0;
GRASOn = 0;

EKF_x_hat_kplus = x_hat_kminus;
EKF_Px_kplus = eye(NumberStates,NumberStates);
EKF_Px_kplus(1:3,1:3) = eye(3,3)*100^2;
EKF_Px_kplus(5:7,5:7) = eye(3,3)*10^2;
EKF_Px_kplus(4,4) = 100^2;
EKF_Px_kplus(8,8) = 10^2;


%% values for RAIM FDE
% 
% PFalseAlarm = 2.22e-8;
% a =  [31.2919346425215; ...
%           35.2463471007126; ...
%           38.4953686364599; ...
%           41.4010185618922; ...
%            44.090770631714; ...
%           46.6276102566246; ...
%           49.0483514687134; ...
%           51.3768201044926; ...
%           53.6295337202138; ...
%           55.8185260702011; ...
%           57.9529000257331; ...
%           60.0397472340809; ...
%           62.0847253044164; ...
%            64.092436838505; ...
%            66.066688033301; ...
%            68.010670670063; ...
%           69.9270936151405; ...
%           71.8182800246577; ...
%           73.6862406187531; ...
%           75.5327299572836];
% 
% lambdatrue = [75.4140000000001; ...
%           80.2599999999999; ...
%           83.9499999999999; ...
%           87.0779999999999; ...
%           89.8469999999999; ...
%           92.3639999999999; ...
%           94.6879999999999; ...
%           96.8589999999999; ...
%                     98.906; ...
%                    100.848; ...
%                    102.698; ...
%                    104.468; ...
%                     106.17; ...
%                     107.81; ...
%                    109.396; ...
%                     110.93; ...
%                     112.42; ...
%                    113.865; ...
%                    115.275; ...
%                    116.647];
% 

PFalseAlarm = 2e-5;
PMissedDetection = 0.047;

a =   [         18.189293484087;...
          21.6395565688235;...
          24.4623581112611;...
          26.9869539367368;...
          29.3272057089061;...
          31.5385558129051];

lambdatrue = [35.3000000000002;...
          38.9000000000003;...
          41.6000000000003;...
          43.8000000000004;...
          45.7000000000004;...
          47.5000000000004];


                   
for SetNumber = 1:NumberSets
    disp(sprintf('Processing set #%d',SetNumber));
    
    NumberEpochsGPS = size(Feb21Data.ObsData(SetNumber).Observations.C1(:,:),2);
    gps_dt = 1;

    NumberSVs = size(Feb21Data.ObsData(SetNumber).Observations.C1(:,:),1);
    SVDontUse = zeros(1,NumberSVs);

    ALPHA = Feb21Data.NavigationData(SetNumber).Iono_ALPHA;
    BETA = Feb21Data.NavigationData(SetNumber).Iono_ALPHA;

    %% arrange pseudorange measurements

    PR_Sim = Feb21Data.ObsData(SetNumber).Observations.C1(:,:)';
    PRR_Sim = -Feb21Data.ObsData(SetNumber).Observations.D1(:,:)' * L1_Wavelength;
    CP_Sim = Feb21Data.ObsData(SetNumber).Observations.L1(:,:)';

    SV_Ephemeris = Feb21Data.NavigationData(SetNumber).Ephemeris;
    GPSTime_Week = Feb21Data.ObsData(SetNumber).GPSTime_Week;
    GPSTime_Sec = Feb21Data.ObsData(SetNumber).GPSTime_Sec;
    ApproxPos = Feb21Data.ObsData(SetNumber).ApproxPos;
    


    % get the initial meridian and prime radii of curvature
    RM = MeridianRadius(InitialPosition(2));
    RP = PrimeRadius(InitialPosition(2));



    Na = NumberStates+ProcessNoiseStates+MeasurementNoiseStates;

    % UKF scaling parameters
    alpha = 0.1; % range: 1e-3 < alpha <= 1
    beta = 2;  % 2 is optimal for gaussian priors
    kapa = 0;  %

    LAMBDA = alpha^2 * (Na+kapa) - Na;



    TimeGPS = [0:gps_dt:NumberEpochsGPS-1];






    GPS_PR_UERE = 7.5;
    GPS_PRR_UERE = 0.2;

    R = eye(MeasurementNoiseStates,MeasurementNoiseStates);
    R(1:6,1:6) = eye(6,6)*GPS_PR_UERE^2;
    R(7:12,7:12) = eye(6,6)*GPS_PRR_UERE^2;





    SV_Azimuth = zeros(NumberEpochsGPS,NumberSVs);
    SV_Elevation =  zeros(NumberEpochsGPS,NumberSVs);

    if(SetNumber == 8 | SetNumber==11)
      SVDontUse(4) = 1;
      SVDontUse(10) = 1;
    end
    
    URALimit = 1;

    %% loop through and check for SV unhealthy (Ephemeris param 27)
    for i=1:size(SV_Ephemeris,1)

       % check for unhealthy flag 
       if(SV_Ephemeris(i,27) ~= 0)
           disp(sprintf('SV %d Ephemeris Unhealthy',SV_Ephemeris(i,1)));
          SVDontUse(SV_Ephemeris(i,1)) = 1;
       end

       % check for large URA

       if(SV_Ephemeris(i,26) > URALimit)
           disp(sprintf('SV %d: URA=%d Exceeds Limit (%d)',SV_Ephemeris(i,1),SV_Ephemeris(i,26),URALimit));
          SVDontUse(SV_Ephemeris(i,1)) = 1;
       end

    end




    %% parameters for EKF

    % power spectral densitites of the process noise - 
    % Sp - position error
    % Sf - clock bias error
    % Sg - clock drift error

    Sp = 10^2;  %% Static receiver
    Sf = 2*2e-16*Speedoflight^2;
    Sg = 8*pi^2*2e-17*Speedoflight^2;
    % Sf = 10;
    % Sg = 1;

    % setup state transition - time invariant if dt is constant
    EKF_phi = eye(NumberStates,NumberStates);
    EKF_phi(1,5) = gps_dt;
    EKF_phi(2,6) = gps_dt;
    EKF_phi(3,7) = gps_dt;
    EKF_phi(4,8) = gps_dt;


    % setup Q and R matrices (process and noise covariance)
    EKF_Q = zeros(NumberStates,NumberStates);

    Sp_dt3 = Sp * (gps_dt ^ 3) / 3;
    Sp_dt2 = Sp * (gps_dt ^ 2) / 2;
    Sp_dt = Sp * gps_dt;
    Sf_dt = Sf * gps_dt;
    Sg_dt3 = Sg * (gps_dt ^ 3) / 3;
    Sg_dt2 = Sg * (gps_dt ^ 2) / 2;
    Sg_dt = Sg * gps_dt;

    EKF_Q(1,1) = Sp_dt3;
    EKF_Q(2,2) = Sp_dt3;
    EKF_Q(3,3) = Sp_dt3;

    EKF_Q(1,5) = Sp_dt2;
    EKF_Q(2,6) = Sp_dt2;
    EKF_Q(3,7) = Sp_dt2;
    EKF_Q(4,8) = Sg_dt2;
    EKF_Q(5,1) = Sp_dt2;
    EKF_Q(6,2) = Sp_dt2;
    EKF_Q(7,3) = Sp_dt2;
    EKF_Q(8,4) = Sg_dt2;

    EKF_Q(4,4) = Sf_dt + Sg_dt3;
    EKF_Q(5,5) = Sp_dt;
    EKF_Q(6,6) = Sp_dt;
    EKF_Q(7,7) = Sp_dt;
    EKF_Q(8,8) = Sg_dt;


    EKF_R = R;
    
    EKF_R_Sub = eye(MeasurementNoiseStates-2,MeasurementNoiseStates-2);
    EKF_R_Sub(1:5,1:5) = eye(5,5)*GPS_PR_UERE^2;
    EKF_R_Sub(6:10,6:10) = eye(5,5)*GPS_PRR_UERE^2;






    UKF_Q = zeros(ProcessNoiseStates,ProcessNoiseStates);
    UKF_Q(1,1) = Sp_dt3;
    UKF_Q(2,2) = Sp_dt3;
    UKF_Q(3,3) = Sp_dt3;
    UKF_Q(4,4) = Sp_dt;
    UKF_Q(5,5) = Sp_dt;
    UKF_Q(6,6) = Sp_dt;

    UKF_Q(1,4) = Sp_dt2;
    UKF_Q(2,5) = Sp_dt2;
    UKF_Q(3,6) = Sp_dt2;

    UKF_Q(4,1) = Sp_dt2;
    UKF_Q(5,2) = Sp_dt2;
    UKF_Q(6,3) = Sp_dt2;

    UKF_Q(8,7) = Sg_dt2;
    UKF_Q(7,8) = Sg_dt2;

    UKF_Q(7,7) = Sf_dt + Sg_dt3;
    UKF_Q(8,8) = Sg_dt;

    %% start gps solutions
    for Epoch_lo = 1:NumberEpochsGPS

        % get the gps time - used for the iono correction
        %GPSTime  = 259199 + Epoch_lo;
        GPSTime = GPSTime_Sec(Epoch_lo);



        % perform state prediction
        % generate augmented state vector
        xa_hat_kminus = [x_hat_kminus;zeros(ProcessNoiseStates,1);zeros(MeasurementNoiseStates,1)];

        Pxa_kminus = [Px_kminus          zeros(NumberStates,NumberStates) zeros(NumberStates,MeasurementNoiseStates);
                  zeros(ProcessNoiseStates,NumberStates) UKF_Q           zeros(ProcessNoiseStates,MeasurementNoiseStates);
                  zeros(MeasurementNoiseStates,NumberStates) zeros(MeasurementNoiseStates,ProcessNoiseStates)  R];





        blah = sqrt(Na+LAMBDA) * chol(Pxa_kminus);
        for i=0:Na
            if i==0
                xs_0 = xa_hat_kminus;
                W_0_m = LAMBDA / (Na + LAMBDA);
                W_0_c = W_0_m + (1 - alpha^2 + beta);
            else
                xs_i(:,i) = xa_hat_kminus + blah(i,:)';
                xs_i(:,i+Na) = xa_hat_kminus - blah(i,:)';
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
                xs_0_k(4) = xs_0(4);
                xs_0_k(5) = xs_0(5);
                xs_0_k(6) = xs_0(6);
                xs_0_k(8) = xs_0(8);
    %             xs_0_k(4) = xs_0(4) + xs_0(9)*gps_dt;
    %             xs_0_k(5) = xs_0(5) + xs_0(10)*gps_dt;
    %             xs_0_k(6) = xs_0(6) + xs_0(11)*gps_dt;
    %             xs_0_k(8) = xs_0(8);
    %             
    %             % acceleration
    %             xs_0_k(9) = xs_0(9);
    %             xs_0_k(10) = xs_0(10);
    %             xs_0_k(11) = xs_0(11);
    %             


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
    %             xs_i_k(4,i) = xs_i(4,i) + xs_i(9,i)*gps_dt;
    %             xs_i_k(5,i) = xs_i(5,i) + xs_i(10,i)*gps_dt;
    %             xs_i_k(6,i) = xs_i(6,i) + xs_i(11,i)*gps_dt;
    %             xs_i_k(8,i) = xs_i(8,i);
                xs_i_k(4,i) = xs_i(4,i);
                xs_i_k(5,i) = xs_i(5,i);
                xs_i_k(6,i) = xs_i(6,i);
                xs_i_k(8,i) = xs_i(8,i);
                % acceleration
    %             xs_i_k(9,i) = xs_i(9,i);
    %             xs_i_k(10,i) = xs_i(10,i);
    %             xs_i_k(11,i) = xs_i(11,i);



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
                PRR_Vec(SVIndex) = PRR_Sim(Epoch_lo,SV);
                if(Epoch_lo == 1)
                    %PRR_Vec(SVIndex) = 0;
                    %PR_csc(SVIndex) = PR_Vec(SVIndex);
                else
                    if(CP_Sim(Epoch_lo-1,SV) == 0)
                        %PRR_Vec(SVIndex) = (CP_Sim(Epoch_lo+1,SV) - CP_Sim(Epoch_lo,SV));
                        %PR_csc(SVIndex) = PR_Vec(SVIndex);
                    else
                        %PRR_Vec(SVIndex) = (CP_Sim(Epoch_lo,SV) - CP_Sim(Epoch_lo-1,SV));

                     %%% PERFORM CARRIER PHASE SMOOTHING of PR
                     %%% HERE


                       % alpha = gps_dt/100;
                        %P_proj = PR_prev + (L1_Wavelength/(2*pi))*(CP_Sim(Epoch_lo) - CP_Sim(Epoch_lo-1));
                        %PR_csc(SVIndex) = alpha * PR_Vec(SVIndex) + (1-alpha)*P_proj;                            

                    end
                end
                % save previous value
                %PR_prev = PR_csc(SVIndex);

                % dont use the carrier smoothed code
                %PR_Vec = PR_csc;                  



    %             SVPos(SVIndex,:) = [SV_X_Data(Epoch_lo,SV) SV_Y_Data(Epoch_lo,SV) SV_Z_Data(Epoch_lo,SV) SV_T_Data(Epoch_lo,SV)*c];
    %             SVVel(SVIndex,:) = [SV_Xvel_Data(Epoch_lo,SV) SV_Yvel_Data(Epoch_lo,SV) SV_Zvel_Data(Epoch_lo,SV) SV_Tvel_Data(Epoch_lo,SV)*c];
    %             SVAcc(SVIndex,:) = [SV_Xacc_Data(Epoch_lo,SV) SV_Yacc_Data(Epoch_lo,SV) SV_Zacc_Data(Epoch_lo,SV) SV_Tacc_Data(Epoch_lo,SV)*c];

                [SVPos(SVIndex,1), SVPos(SVIndex,2), SVPos(SVIndex,3), SVPos(SVIndex,4), ValidPosData(Epoch_lo,SV)] = ...
                   GPSOrbitPropagator(GPSTime_Week(Epoch_lo), GPSTime_Sec(Epoch_lo) - PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris, 7500);

                 [SVVel(SVIndex,1), SVVel(SVIndex,2), SVVel(SVIndex,3), SVVel(SVIndex,4), ...
                  SVAcc(SVIndex,1), SVAcc(SVIndex,2), SVAcc(SVIndex,3), SVAcc(SVIndex,4), ValidVelData(Epoch_lo,SV)] = ...
                  GPSOrbitPropagatorVelocities(GPSTime_Week(Epoch_lo),GPSTime_Sec(Epoch_lo)-PR_Vec(SVIndex)/Speedoflight, SV, SV_Ephemeris);

                SVPos(SVIndex,4) = SVPos(SVIndex,4) * Speedoflight;
                SVVel(SVIndex,4) = SVVel(SVIndex,4) * Speedoflight;
                SVAcc(SVIndex,4) = SVAcc(SVIndex,4) * Speedoflight;

                [SV_Azimuth(Epoch_lo,SV), SV_Elevation(Epoch_lo,SV)] = AzEl(UserPos(1:3), SVPos(SVIndex,1:3));

                % calculate the iono delay correction - single frequency user
                % model from ICD 200
                if(GRASOn == 0)
                    ionodelay = ionomodel(GPSTime, UserPos(1:3), SVPos(SVIndex,1:3), ALPHA, BETA);
                    PR_Vec(SVIndex) = PR_Vec(SVIndex) - ionodelay;
                end

                % calculate hte earth rotation correction as per Kayton pg 228
                % eq 5.67


                PR_Vec_raw(SVIndex) = PR_Vec(SVIndex);  % save a raw (uncorrected copy) of the PR vector for use in the LSQ algorithm later.
                delta_pr_omegaedot = -(OMEGAedot / Speedoflight) * (SVPos(SVIndex,1) *UserPos(2) - SVPos(SVIndex,2) * UserPos(1));
                PR_Vec(SVIndex) = PR_Vec(SVIndex) + delta_pr_omegaedot + SVPos(SVIndex,4);
                PRR_Vec(SVIndex) = PRR_Vec(SVIndex) + SVVel(SVIndex,4);


            end
        end

        NumberGPSMeasurements = length(PR_Vec);



        if NumberGPSMeasurements > 6
            NumberGPSMeasurements = 6;
        end

        % calculate the number of sub-filters required for single satellite
        % failture detection
        NumberSubMeasurements = NumberGPSMeasurements-1;
        NumberSubFilters = nchoosek(NumberGPSMeasurements,NumberSubMeasurements);
        
        % formulate the measurment vector
        y_k = [PR_Vec(1:NumberGPSMeasurements)';PRR_Vec(1:NumberGPSMeasurements)'];


        %% do the full-filter solution
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
            PR_Vec_minus_0(k) = geo_range_to_sat + UserPos(4) + xs_0(NumberStates+ProcessNoiseStates+k);% + UserVel(4) * gps_dt;  % geometric range + c * delta_T
            %predicted relative velocity of sv and receiver
            Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
            PRR_Vec_minus_0(k) = Relative_Velocity(k) + UserVel(4) + xs_0(NumberStates+ProcessNoiseStates+MeasurementNoiseStates/2+k);

            for i=1:2*Na
                % Get User position in ECEF
                UserPos(1:3) = [xs_i_k(1,i);xs_i_k(2,i);xs_i_k(3,i)];
                UserPos(4) = xs_i_k(7,i);
                UserVel(1:3) =  [xs_i_k(4,i);xs_i_k(5,i);xs_i_k(6,i)];
                UserVel(4) = xs_i_k(8,i);


                geo_range_to_sat = sqrt((SVPos(k,1) - UserPos(1))^2 + (SVPos(k,2) - UserPos(2))^2 + (SVPos(k,3) - UserPos(3))^2);
                geo_vel_to_sat = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));

                PR_Vec_minus_i(k,i) = geo_range_to_sat + UserPos(4) + xs_i(NumberStates+ProcessNoiseStates+k,i);% + UserVel(4) * gps_dt + xs_i(NumberStates+ProcessNoiseStates+k-1,i);  % geometric range + c * delta_T

                %predicted relative velocity of sv and receiver
                Relative_Velocity(k) = geo_vel_to_sat/geo_range_to_sat;
                PRR_Vec_minus_i(k,i) = Relative_Velocity(k) + UserVel(4) + xs_i(NumberStates+ProcessNoiseStates+MeasurementNoiseStates/2+k,i);
            end % for i=1:2*Na

            % calculate EKF H-matrix
            UserPos(1:4) = EKF_x_hat_kplus(1:4);
            UserVel(1:4) = EKF_x_hat_kplus(5:8);
            

            
            for m = 1:3
                 ele(m) =  SVPos(k,m) - UserPos(m);
            end    

            r_VecCalc(k) =  norm(ele);   

            EKF_H(k,1) =  -ele(1)/r_VecCalc(k);
            EKF_H(k,2) =  -ele(2)/r_VecCalc(k);
            EKF_H(k,3) =  -ele(3)/r_VecCalc(k);
            EKF_H(k,4) = 1;   
            EKF_H(k,5) = 0;
            EKF_H(k,6) = 0;
            EKF_H(k,7) = 0;
            EKF_H(k,8) = 0;

            EKF_H(k+NumberGPSMeasurements,1) = 0;
            EKF_H(k+NumberGPSMeasurements,2) = 0;
            EKF_H(k+NumberGPSMeasurements,3) = 0;
            EKF_H(k+NumberGPSMeasurements,4) = 0;   
            EKF_H(k+NumberGPSMeasurements,5) = -ele(1)/r_VecCalc(k);
            EKF_H(k+NumberGPSMeasurements,6) = -ele(2)/r_VecCalc(k);
            EKF_H(k+NumberGPSMeasurements,7) = -ele(3)/r_VecCalc(k);
            EKF_H(k+NumberGPSMeasurements,8) = 1;

            % find apriori estimate of pseudorange
            PR_Vec_minus(k) = r_VecCalc(k) + UserPos(4) + UserVel(4) * gps_dt;  % geometric range + c * delta_T
            %predicted relative velocity of sv and receiver
            r_VecCalcVel(k) = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));
            Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);
            PRR_Vec_minus(k) = Relative_Velocity(k) + UserVel(4);


        end  % for k = 1:NumberGPSMeasurements
%         if(Epoch_lo==1 | Epoch_lo==2 | Epoch_lo==3)
%             disp(sprintf('Epoch_lo=%d: UserVel=%f %f %f %f',Epoch_lo,UserVel(1),UserVel(2),UserVel(3),UserVel(4)));
%         end

        SV_Used(Epoch_lo,:) = SV_Vec(1:NumberGPSMeasurements);
        SV_Visible(Epoch_lo,1:length(PR_Vec)) = SV_Vec;
        
        ys_kminus_0 = [PR_Vec_minus_0'; PRR_Vec_minus_0'];
        ys_kminus_i = [PR_Vec_minus_i(:,:); PRR_Vec_minus_i(:,:)];

        % find the sum of sigma points for the measurement prediction


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

        %% Py_kminus is equivalent to HPxH'+R 
        %% Pxy_kminus PxH' 
        %% therefor Hequiv = (inv(Px_kminus)*Pxy_kminus)'
        % calculate 'H' by normalizing Pxy
        Hequiv = (inv(Px_kminus) * Pxy_kminus)';
        Hequiv2(1:NumberGPSMeasurements,1:3) = Hequiv(1:NumberGPSMeasurements,1:3);
        Hequiv2(1:NumberGPSMeasurements,4) = Hequiv(1:NumberGPSMeasurements,7);



        [lat,lon] = ECEF2LLH(ApproxPos);
        Tecef2ned2 = T_ECEF2NED(lon,lat);
        Tecef2ned2(4,1:3) = [0 0 0];
        Tecef2ned2(1:4,4) = [0 0 0 1];
        H_ltp = Hequiv2 * Tecef2ned2';

        % calculate covariance matrix (AA) and transform to local coords
        AA = inv(H_ltp' * H_ltp);

        % calculate DOPs
        var_x = AA(1,1);
        var_y = AA(2,2);
        var_z = AA(3,3);
        var_dt = AA(4,4);
        DOP_UKF(Epoch_lo,1) = sqrt(var_x + var_y + var_z + var_dt); % GDOP
        DOP_UKF(Epoch_lo,2) = sqrt(var_x + var_y + var_z);          % PDOP
        DOP_UKF(Epoch_lo,3) = sqrt(var_x + var_y);                  % HDOP
        DOP_UKF(Epoch_lo,4) = sqrt(var_z);                          % VDOP
        DOP_UKF(Epoch_lo,5) = sqrt(var_dt);                         % TDOP


        % calculate kalman gain
        K_k = Pxy_kminus * inv(Py_kminus);
        %K_k = (Pxy_kminus / Sy_kminus') / Sy_kminus;

        % apply correction
        z_k = y_k - y_hat_kminus;
        z_save(Epoch_lo,:) = z_k;



        x_hat_kplus = x_hat_kminus + K_k * (z_k);

        Px_kplus = Px_kminus - K_k * Py_kminus * K_k';

        %% calculate EKF solution for comparisson

        EKF_x_hat_kminus = EKF_x_hat_kplus;
        EKF_Px_kminus = EKF_Px_kplus;

        %EKF_z = [PR_Vec(1:NumberGPSMeasurements)';PRR_Vec(1:NumberGPSMeasurements)'];
        EKF_z = ([PR_Vec(1:NumberGPSMeasurements) PRR_Vec(1:NumberGPSMeasurements)]' - [PR_Vec_minus(1:NumberGPSMeasurements) PRR_Vec_minus(1:NumberGPSMeasurements)]');

        PR_meas(Epoch_lo,:) = [PR_Vec(1:NumberGPSMeasurements) PRR_Vec(1:NumberGPSMeasurements)]';
        PR_pred(Epoch_lo,:) = [PR_Vec_minus(1:NumberGPSMeasurements) PRR_Vec_minus(1:NumberGPSMeasurements)]';
        
        EKF_z_save(Epoch_lo,:) = EKF_z;

        [EKF_x_hat_kplus, EKF_Px_kplus, EKF_v_out, EKF_s2_out] = ...
            GARD_EvaluateKF(gps_dt, EKF_x_hat_kminus, EKF_Px_kminus, ...
            EKF_phi, EKF_H, EKF_z, EKF_Q, EKF_R);

        
     
        
        %% do the sub filter solutions
        for SubSolution = 1:NumberSubFilters
            
            
            % form the PR and CP vectors with one measurement left
            % out

            PR_SubVec(1:SubSolution-1) = PR_Vec(1:SubSolution-1);
            PR_SubVec(SubSolution:NumberSubMeasurements) = PR_Vec(SubSolution+1:NumberGPSMeasurements);

            PRR_SubVec(1:SubSolution-1) = PRR_Vec(1:SubSolution-1);
            PRR_SubVec(SubSolution:NumberSubMeasurements) = PRR_Vec(SubSolution+1:NumberGPSMeasurements);

            SVPos_Sub(1:SubSolution-1,:) = SVPos(1:SubSolution-1,:);
            SVPos_Sub(SubSolution:NumberSubMeasurements,:) = SVPos(SubSolution+1:NumberGPSMeasurements,:);

            SVVel_Sub(1:SubSolution-1,:) = SVVel(1:SubSolution-1,:);
            SVVel_Sub(SubSolution:NumberSubMeasurements,:) = SVVel(SubSolution+1:NumberGPSMeasurements,:);

                    
            % formulate the measurment vector
            y_k_Sub = [PR_SubVec(1:NumberSubMeasurements)';PRR_SubVec(1:NumberSubMeasurements)'];


            
            for k = 1:NumberSubMeasurements
                % find apriori estimate of pseudorange for each sigma point

                % zero-th sigma point


                % Get User position in ECEF
                UserPos(1:3) = [xs_0_k(1);xs_0_k(2);xs_0_k(3)];
                UserPos(4) = xs_0_k(7);
                UserVel(1:3) =  [xs_0_k(4);xs_0_k(5);xs_0_k(6)];
                UserVel(4) = xs_0_k(8);

                geo_range_to_sat = sqrt((SVPos_Sub(k,1) - UserPos(1))^2 + (SVPos_Sub(k,2) - UserPos(2))^2 + (SVPos_Sub(k,3) - UserPos(3))^2);
                geo_vel_to_sat = (SVVel_Sub(k,1) - UserVel(1))*(SVPos_Sub(k,1)-UserPos(1)) + ...
                                 (SVVel_Sub(k,2) - UserVel(2))*(SVPos_Sub(k,2)-UserPos(2)) + ...
                                 (SVVel_Sub(k,3) - UserVel(3))*(SVPos_Sub(k,3)-UserPos(3));

            %% calculate measurement predicition
                % pseudorange prediction
                PR_Vec_minus_0_Sub(k) = geo_range_to_sat + UserPos(4) + xs_0(NumberStates+ProcessNoiseStates+k);% + UserVel(4) * gps_dt;  % geometric range + c * delta_T
                %predicted relative velocity of sv and receiver
                Relative_Velocity_Sub(k) = geo_vel_to_sat/geo_range_to_sat;
                PRR_Vec_minus_0_Sub(k) = Relative_Velocity_Sub(k) + UserVel(4) + xs_0(NumberStates+ProcessNoiseStates+MeasurementNoiseStates/2+k);

                for i=1:2*Na
                    % Get User position in ECEF
                    UserPos(1:3) = [xs_i_k(1,i);xs_i_k(2,i);xs_i_k(3,i)];
                    UserPos(4) = xs_i_k(7,i);
                    UserVel(1:3) =  [xs_i_k(4,i);xs_i_k(5,i);xs_i_k(6,i)];
                    UserVel(4) = xs_i_k(8,i);


                    geo_range_to_sat = sqrt((SVPos_Sub(k,1) - UserPos(1))^2 + (SVPos_Sub(k,2) - UserPos(2))^2 + (SVPos_Sub(k,3) - UserPos(3))^2);
                    geo_vel_to_sat = (SVVel_Sub(k,1) - UserVel(1))*(SVPos_Sub(k,1)-UserPos(1)) + ...
                                     (SVVel_Sub(k,2) - UserVel(2))*(SVPos_Sub(k,2)-UserPos(2)) + ...
                                     (SVVel_Sub(k,3) - UserVel(3))*(SVPos_Sub(k,3)-UserPos(3));

                    PR_Vec_minus_i_Sub(k,i) = geo_range_to_sat + UserPos(4) + xs_i(NumberStates+ProcessNoiseStates+k,i);% + UserVel(4) * gps_dt + xs_i(NumberStates+ProcessNoiseStates+k-1,i);  % geometric range + c * delta_T

                    %predicted relative velocity of sv and receiver
                    Relative_Velocity_Sub(k) = geo_vel_to_sat/geo_range_to_sat;
                    PRR_Vec_minus_i_Sub(k,i) = Relative_Velocity_Sub(k) + UserVel(4) + xs_i(NumberStates+ProcessNoiseStates+MeasurementNoiseStates/2+k,i);
                end % for i=1:2*Na

                % calculate EKF H-matrix
                UserPos(1:4) = EKF_x_hat_kplus(1:4);
                UserVel(1:4) = EKF_x_hat_kplus(5:8);



                for m = 1:3
                     ele(m) =  SVPos(k,m) - UserPos(m);
                end    

                r_VecCalc_Sub(k) =  norm(ele);   

                EKF_H_Sub(k,1) =  -ele(1)/r_VecCalc(k);
                EKF_H_Sub(k,2) =  -ele(2)/r_VecCalc(k);
                EKF_H_Sub(k,3) =  -ele(3)/r_VecCalc(k);
                EKF_H_Sub(k,4) = 1;   
                EKF_H_Sub(k,5) = 0;
                EKF_H_Sub(k,6) = 0;
                EKF_H_Sub(k,7) = 0;
                EKF_H_Sub(k,8) = 0;

                EKF_H_Sub(k+NumberSubMeasurements,1) = 0;
                EKF_H_Sub(k+NumberSubMeasurements,2) = 0;
                EKF_H_Sub(k+NumberSubMeasurements,3) = 0;
                EKF_H_Sub(k+NumberSubMeasurements,4) = 0;   
                EKF_H_Sub(k+NumberSubMeasurements,5) = -ele(1)/r_VecCalc(k);
                EKF_H_Sub(k+NumberSubMeasurements,6) = -ele(2)/r_VecCalc(k);
                EKF_H_Sub(k+NumberSubMeasurements,7) = -ele(3)/r_VecCalc(k);
                EKF_H_Sub(k+NumberSubMeasurements,8) = 1;

                % find apriori estimate of pseudorange
                PR_Vec_minus_Sub(k) = r_VecCalc(k) + UserPos(4) + UserVel(4) * gps_dt;  % geometric range + c * delta_T
                %predicted relative velocity of sv and receiver
                r_VecCalcVel_Sub(k) = (SVVel_Sub(k,1) - UserVel(1))*(SVPos_Sub(k,1)-UserPos(1)) + ...
                                      (SVVel_Sub(k,2) - UserVel(2))*(SVPos_Sub(k,2)-UserPos(2)) + ...
                                      (SVVel_Sub(k,3) - UserVel(3))*(SVPos_Sub(k,3)-UserPos(3));
                Relative_Velocity_Sub(k) = r_VecCalcVel_Sub(k)/r_VecCalc_Sub(k);
                PRR_Vec_minus_Sub(k) = Relative_Velocity_Sub(k) + UserVel(4);


            end  % for k = 1:NumberGPSMeasurements


            SV_Used(Epoch_lo,:) = SV_Vec(1:NumberGPSMeasurements);
            SV_Visible(Epoch_lo,1:length(PR_Vec)) = SV_Vec;

            ys_kminus_0_Sub = [PR_Vec_minus_0_Sub'; PRR_Vec_minus_0_Sub'];
            ys_kminus_i_Sub = [PR_Vec_minus_i_Sub(:,:); PRR_Vec_minus_i_Sub(:,:)];

            % find the sum of sigma points for the measurement prediction


            W = ones(1,2*Na+1);
            W(1,1) = W_0_m;
            W(1,2:2*Na+1) = W_i_m;

            Y = zeros(MeasurementNoiseStates-2,2*Na+1);
            Y(:,1) = ys_kminus_0_Sub(1:MeasurementNoiseStates-2);
            Y(:,2:2*Na+1) = ys_kminus_i_Sub(1:MeasurementNoiseStates-2,:);
            y_hat_kminus_Sub = (W*Y')';

            nx = length(x_hat_kminus); ny = length(y_hat_kminus_Sub); nw = length(W);
            Py_kminus_Sub=((ones(ny,1)*W).*(Y-y_hat_kminus_Sub*ones(1,nw)))*(Y-y_hat_kminus_Sub*ones(1,nw))';
            Pxy_kminus_Sub =((ones(nx,1)*W).*(X-x_hat_kminus*ones(1,nw)))*(Y-y_hat_kminus_Sub*ones(1,nw))';

            % calculate kalman gain
            K_k_Sub = Pxy_kminus_Sub * inv(Py_kminus_Sub);
            %K_k = (Pxy_kminus / Sy_kminus') / Sy_kminus;

            % apply correction
            z_k_Sub = y_k_Sub - y_hat_kminus_Sub;
            %z_save_Sub(Epoch_lo,:) = z_k;



            x_hat_kplus_Sub(:,SubSolution) = x_hat_kminus + K_k_Sub * (z_k_Sub);

            Px_kplus_Sub(:,:,SubSolution) = Px_kminus - K_k_Sub * Py_kminus_Sub * K_k_Sub';

            
 

            EKF_z_Sub = ([PR_SubVec(1:NumberGPSMeasurements-1) PRR_SubVec(1:NumberGPSMeasurements-1)]' - ...
                [PR_Vec_minus_Sub(1:NumberGPSMeasurements-1) PRR_Vec_minus_Sub(1:NumberGPSMeasurements-1)]');

            [EKF_x_hat_kplus_Sub(:,SubSolution), EKF_Px_kplus_Sub(:,:,SubSolution), EKF_v_out_Sub, EKF_s2_out_Sub] = ...
            GARD_EvaluateKF(gps_dt, EKF_x_hat_kminus, EKF_Px_kminus, ...
            EKF_phi, EKF_H_Sub, EKF_z_Sub, EKF_Q, EKF_R_Sub);
            
            

        end %% sub solutions

        % form the solution separation vectors
        for j = 1:NumberSubFilters
            % note that only the position estimates are used - not
            % clock or velocity
            Beta_ss(:,j) = x_hat_kplus(1:3) - x_hat_kplus_Sub(1:3,j);
            B_ss(:,:,j) = Px_kplus_Sub(1:3,1:3,j) - Px_kplus(1:3,1:3);

            lambda_ss(Epoch_lo,j) = Beta_ss(:,j)' * pinv(B_ss(:,:,j)) * Beta_ss(:,j);

            B_lambda = eigs(B_ss(:,:,j));
            TD(Epoch_lo,j) = sqrt(max(B_lambda)) * abs(norminv(PFalseAlarm/NumberGPSMeasurements,0,1));

            if(lambda_ss(Epoch_lo,j) > TD(Epoch_lo,j))
                disp(sprintf('Fault detected at set %d/%d on Sub-filter %d',SetNumber,Epoch_lo,j));
            end
            
            UKF_PPL_sub(Epoch_lo,j) = sqrt(Px_kplus_Sub(1,1,j)^2 + Px_kplus_Sub(2,2,j)^2 + Px_kplus_Sub(3,3,j)^2) + TD(Epoch_lo,j);


            P_ned = Tecef2ned * Px_kplus_Sub(1:3,1:3,j);

            
            %% this HPL represents the fault-free (H0) hypothesis.  
            UKF_HPL_sub(Epoch_lo,j) = sqrt(P_ned(1,1)^2 + P_ned(2,2)^2) + TD(Epoch_lo,j);
            UKF_VPL_sub(Epoch_lo,j) = sqrt(P_ned(3,3)^2) + TD(Epoch_lo,j);

            %% todo - calculate the Fault-in progress (h1) hypothesis based
            %% HPL_H1.
            
            
            EKF_Beta_ss(:,j) = EKF_x_hat_kplus(1:3) - EKF_x_hat_kplus_Sub(1:3,j);
            EKF_B_ss(:,:,j) = EKF_Px_kplus_Sub(1:3,1:3,j) - EKF_Px_kplus(1:3,1:3);

            EKF_lamda_ss(j) = EKF_Beta_ss(:,j)' * pinv(EKF_B_ss(:,:,j)) * EKF_Beta_ss(:,j);

            EKF_B_lambda = eigs(EKF_B_ss(:,:,j));
            EKF_TD(Epoch_lo,j) = sqrt(max(EKF_B_lambda)) * abs(norminv(PFalseAlarm/NumberSubMeasurements,0,1));

            EKF_PPL_sub(Epoch_lo,j) = sqrt(EKF_Px_kplus_Sub(1,1,j)^2 + EKF_Px_kplus_Sub(2,2,j)^2 + EKF_Px_kplus_Sub(3,3,j)^2) + EKF_TD(Epoch_lo,j);


            EKF_P_ned = Tecef2ned * EKF_Px_kplus_Sub(1:3,1:3,j);

            EKF_HPL_sub(Epoch_lo,j) = sqrt(EKF_P_ned(1,1)^2 + EKF_P_ned(2,2)^2) + EKF_TD(Epoch_lo,j);
            EKF_VPL_sub(Epoch_lo,j) = sqrt(EKF_P_ned(3,3)^2) + EKF_TD(Epoch_lo,j);
            

        end
        UKF_PPL(Epoch_lo) = max(UKF_PPL_sub(Epoch_lo,:));
        UKF_HPL(Epoch_lo) = max(UKF_HPL_sub(Epoch_lo,:));
        UKF_VPL(Epoch_lo) = max(UKF_VPL_sub(Epoch_lo,:));
        
        EKF_PPL(Epoch_lo) = max(EKF_PPL_sub(Epoch_lo,:));
        EKF_HPL(Epoch_lo) = max(EKF_HPL_sub(Epoch_lo,:));
        EKF_VPL(Epoch_lo) = max(EKF_VPL_sub(Epoch_lo,:));
                
        %% calculate solution using least squares for comparisson
        [LSQ_Solution(Epoch_lo,:), LSQ_Variance(Epoch_lo,:), LSQ_NumIterations(Epoch_lo), ...
            LSQ_ResidualVector(Epoch_lo,:), LSQ_M, LSQ_Fail(Epoch_lo), ...
            LSQ_limit(Epoch_lo), LSQ_DOP(Epoch_lo,:)] = GARD_LSQ([ApproxPos 0],NumberGPSMeasurements, PR_Vec_raw(1:NumberGPSMeasurements),SVPos(1:NumberGPSMeasurements,:));

       %% do RAIM Parity on LSQ solution - Note, you must run
       %% GARDSim_CalculateThresholdPbias before using this.
       LSQ_M = LSQ_M*Tecef2ned2';
       [LSQ_RAIM_BadGeometry(Epoch_lo), LSQ_RAIM_ALERT(Epoch_lo), LSQ_RAIM_SLOPE_Max(Epoch_lo), ...
           LSQ_RAIM_r(Epoch_lo), LSQ_RAIM_Td(Epoch_lo), LSQ_RAIM_HPL(Epoch_lo),LSQ_RAIM_VPL(Epoch_lo), LSQ_RAIM_FaultySatFDI(Epoch_lo)] = ...
           GARDSim_RAIMParity(a, lambdatrue, NumberGPSMeasurements,PFalseAlarm,GPS_PR_UERE,556,LSQ_ResidualVector(Epoch_lo,:)',LSQ_M);

        SSE(Epoch_lo) = LSQ_ResidualVector(Epoch_lo,:)*LSQ_ResidualVector(Epoch_lo,:)';

        if(LSQ_RAIM_ALERT(Epoch_lo))
           disp(sprintf('RAIM Error Detected at epoch %d',Epoch_lo)); 
           % pause for 60 seconds
           %pause(10);
        end


        
        

        EKF_x_hat_save(Epoch_lo,:) = EKF_x_hat_kplus;
        EKF_P_save(Epoch_lo,:) = diag(EKF_Px_kplus);



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


        P_save_NED(Epoch_lo,1:3) = abs(Tecef2ned*P_save(Epoch_lo,1:3)');
        S_save_NED(Epoch_lo,1:3) = sqrt(P_save_NED(Epoch_lo,1:3));

        EKF_P_save_NED(Epoch_lo,1:3) = abs(Tecef2ned*EKF_P_save(Epoch_lo,1:3)');
        EKF_S_save_NED(Epoch_lo,1:3) = sqrt(EKF_P_save_NED(Epoch_lo,1:3));

        %% calculate velocity in North-East-Down

        x_save_vel(Epoch_lo,:) = Tecef2ned*x_hat_save(Epoch_lo,4:6)';
        P_save_vel(Epoch_lo,:) = abs(Tecef2ned*P_save(Epoch_lo,4:6)');
        S_save_vel(Epoch_lo,:) = sqrt(P_save_vel(Epoch_lo,1:3));

        %% calculate position error
        ECEFError_UKF(Epoch_lo,:) = x_hat_save(Epoch_lo,1:3) - meanECEF;
        ECEFError_LSQ(Epoch_lo,:) = LSQ_Solution(Epoch_lo,1:3) - meanECEF;
        ECEFError_EKF(Epoch_lo,:) = EKF_x_hat_save(Epoch_lo,1:3) - meanECEF;

        %% calculate velocity in North-East-Down
        [lat,lon] = ECEF2LLH(meanECEF);
        Tecef2ned = T_ECEF2NED(lon,lat);

        NEDError_UKF(Epoch_lo,:) = Tecef2ned * ECEFError_UKF(Epoch_lo,:)';
        NEDError_LSQ(Epoch_lo,:) = Tecef2ned * ECEFError_LSQ(Epoch_lo,:)';
        NEDError_EKF(Epoch_lo,:) = Tecef2ned * ECEFError_EKF(Epoch_lo,:)';

        TwoDError_UKF(Epoch_lo) = norm(NEDError_UKF(Epoch_lo,1:2)); 
        ThreeDError_UKF(Epoch_lo) = norm(NEDError_UKF(Epoch_lo,:));

        TwoDError_EKF(Epoch_lo) = norm(NEDError_EKF(Epoch_lo,1:2)); 
        ThreeDError_EKF(Epoch_lo) = norm(NEDError_EKF(Epoch_lo,:));
        
        TwoDError_LSQ(Epoch_lo) = norm(NEDError_LSQ(Epoch_lo,1:2)); 
        ThreeDError_LSQ(Epoch_lo) = norm(NEDError_LSQ(Epoch_lo,:));
                
        
    %     pos_error_llh(Epoch_lo,:) = x_save_llh(Epoch_lo,:) - pos_truth_llh(2:4,Epoch_lo*100)';
    %     pos_error_llh(Epoch_lo,1) = pos_error_llh(Epoch_lo,1) * RM;
    %     pos_error_llh(Epoch_lo,2) = pos_error_llh(Epoch_lo,2) * RP * cos(pos_truth_llh(2,Epoch_lo*100));


    %     vel_error_ned(Epoch_lo,:) = x_save_vel(Epoch_lo,:) - vel_truth(2:4,Epoch_lo*100)';

    %% uncomment the below line if hte 11-state filter is used.
    %    x_save_acc(Epoch_lo,:) = Tecef2ned*x_hat_save(Epoch_lo,9:11)';

        %% update R-matrix
    %     if(Epoch_lo > 10)
    %         for k = 1:NumberGPSMeasurements
    %             R(k,k) = std(abs(z_save(:,k)))^2;
    %             R(k+NumberGPSMeasurements,k+NumberGPSMeasurements) = std(abs(z_save(:,k+NumberGPSMeasurements)))^2;
    %         end
    %     end
    %     



        if(mod(Epoch_lo,100)==0)
            disp(sprintf('Completed Set #%d, Epoch %d',SetNumber,Epoch_lo));
        end



    end

    
    %% save results
    SaveIndexStart = ((SetNumber-1)*NumberEpochsGPS)+1;
    SaveIndexEnd = ((SetNumber)*NumberEpochsGPS);
    All_t_save(SaveIndexStart:SaveIndexEnd,:) = t_save'+(SetNumber-1)*NumberEpochsGPS*gps_dt;
    
    All_GPSTime_Sec(SaveIndexStart:SaveIndexEnd,:) = GPSTime_Sec';
    All_GPSTime_Week(SaveIndexStart:SaveIndexEnd,:) = GPSTime_Week';
    
    All_NEDError_UKF(SaveIndexStart:SaveIndexEnd,:) = NEDError_UKF;
    All_NEDError_EKF(SaveIndexStart:SaveIndexEnd,:) = NEDError_EKF;
    All_NEDError_LSQ(SaveIndexStart:SaveIndexEnd,:) = NEDError_LSQ;
    
    All_UKF_S_save_NED(SaveIndexStart:SaveIndexEnd,:) = S_save_NED;
    All_EKF_S_save_NED(SaveIndexStart:SaveIndexEnd,:) = EKF_S_save_NED;
    All_LSQ_DOP(SaveIndexStart:SaveIndexEnd,:) = LSQ_DOP;
    
    
    All_x_hat_save(SaveIndexStart:SaveIndexEnd,:) = x_hat_save;
    All_LSQ_Solution(SaveIndexStart:SaveIndexEnd,:) = LSQ_Solution;
    All_EKF_x_hat_save(SaveIndexStart:SaveIndexEnd,:) = EKF_x_hat_save;
    
    
    All_z_save(SaveIndexStart:SaveIndexEnd,:) = z_save;
    All_LSQ_ResidualVector(SaveIndexStart:SaveIndexEnd,:) = LSQ_ResidualVector;
    All_EKF_z_save(SaveIndexStart:SaveIndexEnd,:) = EKF_z_save;
    
    All_PR_meas(SaveIndexStart:SaveIndexEnd,:) = PR_meas;
    All_PR_pred(SaveIndexStart:SaveIndexEnd,:) = PR_pred;
    
    All_LSQ_RAIM_HPL(SaveIndexStart:SaveIndexEnd,:) = LSQ_RAIM_HPL';
    All_LSQ_RAIM_VPL(SaveIndexStart:SaveIndexEnd,:) = LSQ_RAIM_VPL';
    
    
    All_UKF_HPL(SaveIndexStart:SaveIndexEnd,:) = UKF_HPL';
    All_UKF_VPL(SaveIndexStart:SaveIndexEnd,:) = UKF_VPL';
    All_UKF_PPL(SaveIndexStart:SaveIndexEnd,:) = UKF_PPL';
    
    All_EKF_HPL(SaveIndexStart:SaveIndexEnd,:) = EKF_HPL';
    All_EKF_VPL(SaveIndexStart:SaveIndexEnd,:) = EKF_VPL';
    All_EKF_PPL(SaveIndexStart:SaveIndexEnd,:) = EKF_PPL';
    
    All_lambda_ss(SaveIndexStart:SaveIndexEnd,:) = lambda_ss;
    All_TD(SaveIndexStart:SaveIndexEnd,:) = TD;
    
    
    All_TwoDError_UKF(SaveIndexStart:SaveIndexEnd,:) = TwoDError_UKF';
    All_ThreeDError_UKF(SaveIndexStart:SaveIndexEnd,:) = ThreeDError_UKF';
    All_TwoDError_EKF(SaveIndexStart:SaveIndexEnd,:) = TwoDError_EKF';
    All_ThreeDError_EKF(SaveIndexStart:SaveIndexEnd,:) = ThreeDError_EKF';
    All_TwoDError_LSQ(SaveIndexStart:SaveIndexEnd,:) = TwoDError_LSQ';
    All_ThreeDError_LSQ(SaveIndexStart:SaveIndexEnd,:) = ThreeDError_LSQ';
    
    %% plot static positioning errors - UKF vs LSQ
%     figure();
%     subplot(3,1,1),plot(t_save/60,abs(NEDError_UKF(:,1)),'b'); hold on;
%     subplot(3,1,1),plot(t_save/60,2*S_save_NED(:,1),'b--');
%     %subplot(3,1,1),plot(t_save/60,-2*S_save_NED(:,1),'b--');
%     subplot(3,1,1),plot(t_save/60,abs(NEDError_LSQ(:,1)),'r');
%     subplot(3,1,1),plot(t_save/60,LSQ_DOP(:,3)*GPS_PR_UERE*2,'r--');
%     %subplot(3,1,1),plot(t_save/60,-LSQ_DOP(:,3)*GPS_PR_UERE*2,'r--');
% 
%     subplot(3,1,1),plot(t_save/60,abs(NEDError_EKF(:,1)),'g');
%     subplot(3,1,1),plot(t_save/60,2*EKF_S_save_NED(:,1),'g--');
% 
%     %axis([t_save(1)/60,t_save(Epoch_lo)/60,00,40]);
%     grid on;
%     title('Static Positioning UKF and LSQ Errors');
%     legend('UKF Error','UKF 2-\sigma bound','LSQ Error','LSQ 2-\sigma bound','EKF Error','EKF 2-\sigma bound');
%     ylabel('North Error (m)');
%     xlabel('Test Time (mins)');
% 
%     subplot(3,1,2),plot(t_save/60,abs(NEDError_UKF(:,2)),'b'); hold on;
%     subplot(3,1,2),plot(t_save/60,2*S_save_NED(:,2),'b--');
%     %subplot(3,1,2),plot(t_save/60,-2*S_save_NED(:,2),'b--');
%     subplot(3,1,2),plot(t_save/60,abs(NEDError_LSQ(:,2)),'r'); 
%     subplot(3,1,2),plot(t_save/60,LSQ_DOP(:,3)*GPS_PR_UERE*2,'r--');
%     %subplot(3,1,2),plot(t_save/60,-LSQ_DOP(:,3)*GPS_PR_UERE*2,'r--');
% 
%     subplot(3,1,2),plot(t_save/60,abs(NEDError_EKF(:,2)),'g');
%     subplot(3,1,2),plot(t_save/60,2*EKF_S_save_NED(:,2),'g--');
% 
%     %axis([t_save(1)/60,t_save(Epoch_lo)/60,0,40]);
%     grid on;
%     ylabel('East Error (m)');
%     xlabel('Test Time (mins)');
%     subplot(3,1,3),plot(t_save/60,abs(NEDError_UKF(:,3)),'b'); hold on;
%     subplot(3,1,3),plot(t_save/60,2*S_save_NED(:,3),'b--');
%     %subplot(3,1,3),plot(t_save/60,-2*S_save_NED(:,3),'b--');
%     subplot(3,1,3),plot(t_save/60,abs(NEDError_LSQ(:,3)),'r');
%     subplot(3,1,3),plot(t_save/60,LSQ_DOP(:,4)*GPS_PR_UERE*2,'r--');
%     %subplot(3,1,3),plot(t_save/60,-LSQ_DOP(:,4)*GPS_PR_UERE*2,'r--');
% 
%     subplot(3,1,3),plot(t_save/60,abs(NEDError_EKF(:,3)),'g');
%     subplot(3,1,3),plot(t_save/60,2*EKF_S_save_NED(:,3),'g--');
% 
%     %axis([t_save(1)/60,t_save(Epoch_lo)/60,0,40]);
%     grid on;
%     ylabel('Down Error (m)');
%     xlabel('Test Time (mins)');



    % 
    % figure();
    % subplot(3,1,1),plot(t_save/60,NEDError_UKF(:,1),'b'); hold on;
    % subplot(3,1,1),plot(t_save/60,2*S_save_NED(:,1),'r');
    % subplot(3,1,1),plot(t_save/60,-2*S_save_NED(:,1),'r'); 
    % 
    % grid on;
    % title('UKF Error Bounds');
    % ylabel('North Error (m)');
    % subplot(3,1,2),plot(t_save/60,NEDError_UKF(:,2),'b'); hold on;
    % subplot(3,1,2),plot(t_save/60,2*S_save_NED(:,2),'r'); 
    % subplot(3,1,2),plot(t_save/60,-2*S_save_NED(:,2),'r'); 
    % 
    % grid on;
    % ylabel('East Error (m)');
    % subplot(3,1,3),plot(t_save/60,NEDError_UKF(:,3),'b'); hold on;
    % subplot(3,1,3),plot(t_save/60,2*S_save_NED(:,3),'r');
    % subplot(3,1,3),plot(t_save/60,-2*S_save_NED(:,3),'r');
    % 
    % grid on;
    % ylabel('Down Error (m)');
    % xlabel('Test Time (mins)');

    % 
    % figure();
    % plot(t_save/60,x_hat_save(:,7))
    % grid on;
    % xlabel('Test Time (mins)');
    % ylabel('Clock Estimates (m)');
    % title('Receiver Clock Bias');

    %% plot residuals
    % figure();
    % subplot(2,1,1),plot(t_save/60,z_save(:,1:NumberGPSMeasurements));
    % axis([0 Epoch_lo/60 -3 3]);
    % ylabel('UKF Residuals');
    % xlabel('Test Time (mins)');
    % legend(num2str((SV_Vec(1:NumberGPSMeasurements)')));
    % grid on;
    % subplot(2,1,2),plot(t_save/60,LSQ_ResidualVector(:,1:NumberGPSMeasurements));
    % %axis([0 Epoch_lo/60 -3 3]);
    % ylabel('LSQ Residuals');
    % xlabel('Test Time (mins)');
    % legend(num2str((SV_Vec(1:NumberGPSMeasurements)')));
    % grid on;


%     for i=1:NumberGPSMeasurements
%         figure();
%         plot(t_save/60,z_save(:,i),'b'); hold on;
%         plot(t_save/60,LSQ_ResidualVector(:,i),'r');
%         plot(t_save/60,EKF_z_save(:,i),'g');
%         xlabel('Test Time (mins)');
%         ylabel('Residual (m)');
%         legend('UKF','LSQ','EKF');
%         title(sprintf('Measurement Residual for SV%d',SV_Vec(i)))
%         grid on;
%         %axis([t_save(1)/60 t_save(NumberEpochsGPS)/60 -2 2]);
%     end

    % figure();
    % plot(t_save/60,z_save(:,1:NumberGPSMeasurements));
    % hold on;
    % plot(t_save/60,LSQ_ResidualVector(:,1:NumberGPSMeasurements));
    % legend(num2str((SV_Vec(1:NumberGPSMeasurements)')));
    % grid on;
    % xlabel('Test Time (min)');
    % ylabel('Residuals');

%     %% plot clock estimate
%     figure();
%     plot(t_save/60,x_hat_save(:,7),'b');
%     hold on;
%     plot(t_save/60,LSQ_Solution(:,4),'r');
%     grid on;
%     plot(t_save/60,EKF_x_hat_save(:,4),'g');
%     xlabel('Test Time (mins)');
%     ylabel('Clock Estimates (m)');
%     %axis([0 Epoch_lo/60 6 16]);
%     title('Receiver Clock Bias');
%     legend('UKF','LSQ','EKF');

    % %% plot DOPs
    % figure();
    % subplot(3,1,1),plot(t_save/60,DOP_UKF(:,3),'b'); hold on;
    % subplot(3,1,1),plot(t_save/60,LSQ_DOP(:,3),'r'); grid on;
    % ylabel('HDOP');
    % title('DOP Values Calculated by UKF and Least-Squares');
    % legend('UKF','LSQ');
    % subplot(3,1,2),plot(t_save/60,DOP_UKF(:,4),'b'); hold on;
    % subplot(3,1,2),plot(t_save/60,LSQ_DOP(:,4),'r'); grid on;
    % ylabel('VDOP');
    % 
    % subplot(3,1,3),plot(t_save/60,DOP_UKF(:,5),'b'); hold on;
    % subplot(3,1,3),plot(t_save/60,LSQ_DOP(:,5),'r'); grid on;
    % ylabel('TDOP');
    % xlabel('Test Time (minutes)');



    %% calculate the RMS error
    %Height_RMS = sqrt(sum(Height_error.^2)/NumberEpochsINS)'
    North_RMS_LSQ(SetNumber) = sqrt(sum(NEDError_LSQ(:,1).^2)/NumberEpochsGPS);
    North_RMS_EKF(SetNumber) = sqrt(sum(NEDError_EKF(:,1).^2)/NumberEpochsGPS);
    North_RMS_UKF(SetNumber) = sqrt(sum(NEDError_UKF(:,1).^2)/NumberEpochsGPS);

    East_RMS_LSQ(SetNumber) = sqrt(sum(NEDError_LSQ(:,2).^2)/NumberEpochsGPS);
    East_RMS_EKF(SetNumber) = sqrt(sum(NEDError_EKF(:,2).^2)/NumberEpochsGPS);
    East_RMS_UKF(SetNumber) = sqrt(sum(NEDError_UKF(:,2).^2)/NumberEpochsGPS);

    Down_RMS_LSQ(SetNumber) = sqrt(sum(NEDError_LSQ(:,3).^2)/NumberEpochsGPS);
    Down_RMS_EKF(SetNumber) = sqrt(sum(NEDError_EKF(:,3).^2)/NumberEpochsGPS);
    Down_RMS_UKF(SetNumber) = sqrt(sum(NEDError_UKF(:,3).^2)/NumberEpochsGPS);


    North_RMS_LSQ_sub(SetNumber) = sqrt(sum(NEDError_LSQ(30:length(NEDError_LSQ),1).^2)/(NumberEpochsGPS-30));
    North_RMS_EKF_sub(SetNumber) = sqrt(sum(NEDError_EKF(30:length(NEDError_LSQ),1).^2)/(NumberEpochsGPS-30));
    North_RMS_UKF_sub(SetNumber) = sqrt(sum(NEDError_UKF(30:length(NEDError_LSQ),1).^2)/(NumberEpochsGPS-30));
    
    East_RMS_LSQ_sub(SetNumber) = sqrt(sum(NEDError_LSQ(30:length(NEDError_LSQ),2).^2)/(NumberEpochsGPS-30));
    East_RMS_EKF_sub(SetNumber) = sqrt(sum(NEDError_EKF(30:length(NEDError_LSQ),2).^2)/(NumberEpochsGPS-30));
    East_RMS_UKF_sub(SetNumber) = sqrt(sum(NEDError_UKF(30:length(NEDError_LSQ),2).^2)/(NumberEpochsGPS-30));

    Down_RMS_LSQ_sub(SetNumber) = sqrt(sum(NEDError_LSQ(30:length(NEDError_LSQ),3).^2)/(NumberEpochsGPS-30));
    Down_RMS_EKF_sub(SetNumber) = sqrt(sum(NEDError_EKF(30:length(NEDError_LSQ),3).^2)/(NumberEpochsGPS-30));
    Down_RMS_UKF_sub(SetNumber) = sqrt(sum(NEDError_UKF(30:length(NEDError_LSQ),3).^2)/(NumberEpochsGPS-30));


    
end % numbersets


% plot(All_NEDError_LSQ(:,1))
% hold on;
% plot(All_NEDError_UKF(:,1),'r');
% plot(All_NEDError_EKF(:,1),'g');

figure();
subplot(3,1,1),plot(All_t_save/60,abs(All_NEDError_UKF(:,1)),'b'); hold on;
subplot(3,1,1),plot(All_t_save/60,2*All_UKF_S_save_NED(:,1),'b--');
subplot(3,1,1),plot(All_t_save/60,abs(All_NEDError_LSQ(:,1)),'r');
subplot(3,1,1),plot(All_t_save/60,All_LSQ_DOP(:,3)*GPS_PR_UERE*2,'r--');
subplot(3,1,1),plot(All_t_save/60,abs(All_NEDError_EKF(:,1)),'g');
subplot(3,1,1),plot(All_t_save/60,2*All_EKF_S_save_NED(:,1),'g--');
axis([All_t_save(1)/60 All_t_save(length(All_t_save))/60 0 50]);
grid on;
title('Static Positioning UKF and LSQ Errors');
legend('UKF Error','UKF 2-\sigma bound','LSQ Error','LSQ 2-\sigma bound','EKF Error','EKF 2-\sigma bound');
ylabel('North Error (m)');
xlabel('Test Time (mins)');

subplot(3,1,2),plot(All_t_save/60,abs(All_NEDError_UKF(:,2)),'b'); hold on;
subplot(3,1,2),plot(All_t_save/60,2*All_UKF_S_save_NED(:,2),'b--');
subplot(3,1,2),plot(All_t_save/60,abs(All_NEDError_LSQ(:,2)),'r'); 
subplot(3,1,2),plot(All_t_save/60,All_LSQ_DOP(:,3)*GPS_PR_UERE*2,'r--');
subplot(3,1,2),plot(All_t_save/60,abs(All_NEDError_EKF(:,2)),'g');
subplot(3,1,2),plot(All_t_save/60,2*All_EKF_S_save_NED(:,2),'g--');
axis([All_t_save(1)/60 All_t_save(length(All_t_save))/60 0 50]);
grid on;
ylabel('East Error (m)');
xlabel('Test Time (mins)');
subplot(3,1,3),plot(All_t_save/60,abs(All_NEDError_UKF(:,3)),'b'); hold on;
subplot(3,1,3),plot(All_t_save/60,2*All_UKF_S_save_NED(:,3),'b--');
subplot(3,1,3),plot(All_t_save/60,abs(All_NEDError_LSQ(:,3)),'r');
subplot(3,1,3),plot(All_t_save/60,All_LSQ_DOP(:,4)*GPS_PR_UERE*2,'r--');

subplot(3,1,3),plot(All_t_save/60,abs(All_NEDError_EKF(:,3)),'g');
subplot(3,1,3),plot(All_t_save/60,2*All_EKF_S_save_NED(:,3),'g--');
axis([All_t_save(1)/60 All_t_save(length(All_t_save))/60 0 50]);
grid on;
ylabel('Down Error (m)');
xlabel('Test Time (mins)');

%% calculate the percentage of samples with error > 1-sigma, 2-sigma,
%% 3-sigma
UKF_Bound_Performance_1s(1) = 1-sum(abs(All_NEDError_UKF(:,1)) > All_UKF_S_save_NED(:,1))/length(All_NEDError_UKF(:,1));
UKF_Bound_Performance_1s(2) = 1-sum(abs(All_NEDError_UKF(:,2)) > All_UKF_S_save_NED(:,2))/length(All_NEDError_UKF(:,2));
UKF_Bound_Performance_1s(3) = 1-sum(abs(All_NEDError_UKF(:,3)) > All_UKF_S_save_NED(:,3))/length(All_NEDError_UKF(:,3))

EKF_Bound_Performance_1s(1) = 1-sum(abs(All_NEDError_EKF(:,1)) > All_EKF_S_save_NED(:,1))/length(All_NEDError_EKF(:,1));
EKF_Bound_Performance_1s(2) = 1-sum(abs(All_NEDError_EKF(:,2)) > All_EKF_S_save_NED(:,2))/length(All_NEDError_EKF(:,2));
EKF_Bound_Performance_1s(3) = 1-sum(abs(All_NEDError_EKF(:,3)) > All_EKF_S_save_NED(:,3))/length(All_NEDError_EKF(:,3))

LSQ_Bound_Performance_1s(1) = 1-sum(abs(All_NEDError_LSQ(:,1)) > All_LSQ_DOP(:,3)*GPS_PR_UERE)/length(All_NEDError_LSQ(:,1));
LSQ_Bound_Performance_1s(2) = 1-sum(abs(All_NEDError_LSQ(:,2)) > All_LSQ_DOP(:,3)*GPS_PR_UERE)/length(All_NEDError_LSQ(:,2));
LSQ_Bound_Performance_1s(3) = 1-sum(abs(All_NEDError_LSQ(:,3)) > All_LSQ_DOP(:,4)*GPS_PR_UERE)/length(All_NEDError_LSQ(:,3))


UKF_Bound_Performance_2s(1) = 1-sum(abs(All_NEDError_UKF(:,1)) > 2*All_UKF_S_save_NED(:,1))/length(All_NEDError_UKF(:,1));
UKF_Bound_Performance_2s(2) = 1-sum(abs(All_NEDError_UKF(:,2)) > 2*All_UKF_S_save_NED(:,2))/length(All_NEDError_UKF(:,2));
UKF_Bound_Performance_2s(3) = 1-sum(abs(All_NEDError_UKF(:,3)) > 2*All_UKF_S_save_NED(:,3))/length(All_NEDError_UKF(:,3))

EKF_Bound_Performance_2s(1) = 1-sum(abs(All_NEDError_EKF(:,1)) > 2*All_EKF_S_save_NED(:,1))/length(All_NEDError_EKF(:,1));
EKF_Bound_Performance_2s(2) = 1-sum(abs(All_NEDError_EKF(:,2)) > 2*All_EKF_S_save_NED(:,2))/length(All_NEDError_EKF(:,2));
EKF_Bound_Performance_2s(3) = 1-sum(abs(All_NEDError_EKF(:,3)) > 2*All_EKF_S_save_NED(:,3))/length(All_NEDError_EKF(:,3))

LSQ_Bound_Performance_2s(1) = 1-sum(abs(All_NEDError_LSQ(:,1)) > All_LSQ_DOP(:,3)*GPS_PR_UERE*2)/length(All_NEDError_LSQ(:,1));
LSQ_Bound_Performance_2s(2) = 1-sum(abs(All_NEDError_LSQ(:,2)) > All_LSQ_DOP(:,3)*GPS_PR_UERE*2)/length(All_NEDError_LSQ(:,2));
LSQ_Bound_Performance_2s(3) = 1-sum(abs(All_NEDError_LSQ(:,3)) > All_LSQ_DOP(:,4)*GPS_PR_UERE*2)/length(All_NEDError_LSQ(:,3))

UKF_Bound_Performance_3s(1) = 1-sum(abs(All_NEDError_UKF(:,1)) > 3*All_UKF_S_save_NED(:,1))/length(All_NEDError_UKF(:,1));
UKF_Bound_Performance_3s(2) = 1-sum(abs(All_NEDError_UKF(:,2)) > 3*All_UKF_S_save_NED(:,2))/length(All_NEDError_UKF(:,2));
UKF_Bound_Performance_3s(3) = 1-sum(abs(All_NEDError_UKF(:,3)) > 3*All_UKF_S_save_NED(:,3))/length(All_NEDError_UKF(:,3))

EKF_Bound_Performance_3s(1) = 1-sum(abs(All_NEDError_EKF(:,1)) > 3*All_EKF_S_save_NED(:,1))/length(All_NEDError_EKF(:,1));
EKF_Bound_Performance_3s(2) = 1-sum(abs(All_NEDError_EKF(:,2)) > 3*All_EKF_S_save_NED(:,2))/length(All_NEDError_EKF(:,2));
EKF_Bound_Performance_3s(3) = 1-sum(abs(All_NEDError_EKF(:,3)) > 3*All_EKF_S_save_NED(:,3))/length(All_NEDError_EKF(:,3))

LSQ_Bound_Performance_3s(1) = 1-sum(abs(All_NEDError_LSQ(:,1)) > All_LSQ_DOP(:,3)*GPS_PR_UERE*3)/length(All_NEDError_LSQ(:,1));
LSQ_Bound_Performance_3s(2) = 1-sum(abs(All_NEDError_LSQ(:,2)) > All_LSQ_DOP(:,3)*GPS_PR_UERE*3)/length(All_NEDError_LSQ(:,2));
LSQ_Bound_Performance_3s(3) = 1-sum(abs(All_NEDError_LSQ(:,3)) > All_LSQ_DOP(:,4)*GPS_PR_UERE*3)/length(All_NEDError_LSQ(:,3))

%% calculate full day RMS errors
All_North_RMS_LSQ = sqrt(sum(All_NEDError_LSQ(:,1).^2)/length(All_NEDError_LSQ))
All_North_RMS_EKF = sqrt(sum(All_NEDError_EKF(:,1).^2)/length(All_NEDError_EKF))
All_North_RMS_UKF = sqrt(sum(All_NEDError_UKF(:,1).^2)/length(All_NEDError_UKF))

All_East_RMS_LSQ = sqrt(sum(All_NEDError_LSQ(:,2).^2)/length(All_NEDError_LSQ))
All_East_RMS_EKF = sqrt(sum(All_NEDError_EKF(:,2).^2)/length(All_NEDError_EKF))
All_East_RMS_UKF = sqrt(sum(All_NEDError_UKF(:,2).^2)/length(All_NEDError_UKF))

All_Down_RMS_LSQ = sqrt(sum(All_NEDError_LSQ(:,3).^2)/length(All_NEDError_LSQ))
All_Down_RMS_EKF = sqrt(sum(All_NEDError_EKF(:,3).^2)/length(All_NEDError_EKF))
All_Down_RMS_UKF = sqrt(sum(All_NEDError_UKF(:,3).^2)/length(All_NEDError_UKF))


All_North_RMS_LSQ_sub = sqrt(sum(All_NEDError_LSQ(100:86400,1).^2)/(length(All_NEDError_LSQ)-100))
All_North_RMS_EKF_sub = sqrt(sum(All_NEDError_EKF(100:86400,1).^2)/(length(All_NEDError_EKF)-100))
All_North_RMS_UKF_sub = sqrt(sum(All_NEDError_UKF(100:86400,1).^2)/(length(All_NEDError_UKF)-100))

All_East_RMS_LSQ_sub = sqrt(sum(All_NEDError_LSQ(100:86400,2).^2)/(length(All_NEDError_LSQ)-100))
All_East_RMS_EKF_sub = sqrt(sum(All_NEDError_EKF(100:86400,2).^2)/(length(All_NEDError_EKF)-100))
All_East_RMS_UKF_sub = sqrt(sum(All_NEDError_UKF(100:86400,2).^2)/(length(All_NEDError_UKF)-100))

All_Down_RMS_LSQ_sub = sqrt(sum(All_NEDError_LSQ(100:86400,3).^2)/(length(All_NEDError_LSQ)-100))
All_Down_RMS_EKF_sub = sqrt(sum(All_NEDError_EKF(100:86400,3).^2)/(length(All_NEDError_EKF)-100))
All_Down_RMS_UKF_sub = sqrt(sum(All_NEDError_UKF(100:86400,3).^2)/(length(All_NEDError_UKF)-100))

All_TwoDRMS_LSQ = sqrt(sum(All_TwoDError_LSQ.^2)/length(All_TwoDError_LSQ));
All_TwoDRMS_EKF = sqrt(sum(All_TwoDError_EKF.^2)/length(All_TwoDError_LSQ));
All_TwoDRMS_UKF = sqrt(sum(All_TwoDError_UKF.^2)/length(All_TwoDError_LSQ));

All_TwoDRMS_LSQ_sub = sqrt(sum(All_TwoDError_LSQ(100:86400).^2)/length(All_TwoDError_LSQ));
All_TwoDRMS_EKF_sub = sqrt(sum(All_TwoDError_EKF(100:86400).^2)/length(All_TwoDError_LSQ));
All_TwoDRMS_UKF_sub = sqrt(sum(All_TwoDError_UKF(100:86400).^2)/length(All_TwoDError_LSQ));


All_ThreeDRMS_LSQ_sub = sqrt(sum(All_ThreeDError_LSQ(100:86400).^2)/length(All_ThreeDError_LSQ));
All_ThreeDRMS_EKF_sub = sqrt(sum(All_ThreeDError_EKF(100:86400).^2)/length(All_ThreeDError_LSQ));
All_ThreeDRMS_UKF_sub = sqrt(sum(All_ThreeDError_UKF(100:86400).^2)/length(All_ThreeDError_LSQ));


%% plot variation in rms performance
figure();
plot(North_RMS_UKF,'b','LineWidth',1.5)
hold on;
plot(North_RMS_EKF,'g','LineWidth',1.5)
plot(North_RMS_LSQ,'r','LineWidth',1.5)
grid on;
legend('UKF','EKF','LSQ');
xlabel('Hour');
ylabel('RMS Error (m)');
title('North Channel RMS Variation over a day');
axis([1 24 0 35]);


figure();
plot(East_RMS_UKF,'b','LineWidth',1.5)
hold on;
plot(East_RMS_EKF,'g','LineWidth',1.5)
plot(East_RMS_LSQ,'r','LineWidth',1.5)
grid on;
legend('UKF','EKF','LSQ');
xlabel('Hour');
ylabel('RMS Error');
title('East Channel RMS Variation over a day');
axis([1 24 0 35]);


figure();
plot(Down_RMS_UKF,'b','LineWidth',1.5)
hold on;
plot(Down_RMS_EKF,'g','LineWidth',1.5)
plot(Down_RMS_LSQ,'r','LineWidth',1.5)
grid on;
legend('UKF','EKF','LSQ');
xlabel('Hour');
ylabel('RMS Error');
title('Vertical Channel RMS Variation over a day');
axis([1 24 0 35]);

%% plot a chart of covariance versus position error - similar to a stanford
%% plot
% 
% figure();
% plot(abs(All_NEDError_UKF(:,1)),2*All_UKF_S_save_NED(:,1),'b.');
% hold on;
% plot(abs(All_NEDError_LSQ(:,1)),All_LSQ_DOP(:,3)*GPS_PR_UERE*2,'r.');
% plot(abs(All_NEDError_EKF(:,1)),2*All_EKF_S_save_NED(:,1),'g.');
% plot([0:1:30],[0:1:30],'k-');
% %axis([0 30 0 30]);
% grid on;
% xlabel('Position Error (m)');
% ylabel('2-\sigma Error Bound (m)');
% title('North Error Protection');
% legend('UKF','LSQ','EKF');
% 
% figure();
% plot(abs(All_NEDError_UKF(:,2)),2*All_UKF_S_save_NED(:,2),'b.');
% hold on;
% plot(abs(All_NEDError_LSQ(:,2)),All_LSQ_DOP(:,3)*GPS_PR_UERE*2,'r.');
% plot(abs(All_NEDError_EKF(:,2)),2*All_EKF_S_save_NED(:,2),'g.');
% plot([0:1:30],[0:1:30],'k-');
% %axis([0 30 0 30]);
% grid on;
% xlabel('Position Error (m)');
% ylabel('2-\sigma Error Bound (m)');
% title('East Error Protection');
% legend('UKF','LSQ','EKF');
% 
% figure();
% plot(abs(All_NEDError_UKF(:,3)),2*All_UKF_S_save_NED(:,3),'b.');
% hold on;
% plot(abs(All_NEDError_LSQ(:,3)),All_LSQ_DOP(:,4)*GPS_PR_UERE*2,'r.');
% plot(abs(All_NEDError_EKF(:,3)),2*All_EKF_S_save_NED(:,3),'g.');
% plot([0:1:30],[0:1:30],'k-');
% %axis([0 30 0 30]);
% grid on;
% xlabel('Position Error (m)');
% ylabel('2-\sigma Error Bound (m)');
% title('Down Error Protection');
% legend('UKF','LSQ','EKF');

% %% plot clock estimate
% figure();
% plot(All_t_save/60,All_x_hat_save(:,7),'b');
% hold on;
% plot(All_t_save/60,All_LSQ_Solution(:,4),'r');
% grid on;
% plot(All_t_save/60,All_EKF_x_hat_save(:,4),'g');
% xlabel('Test Time (mins)');
% ylabel('Clock Estimates (m)');
% title('Receiver Clock Bias');
% legend('UKF','LSQ','EKF');


%% 
figure;
plot(All_t_save/60,All_UKF_HPL);
hold on;
plot(All_t_save/60,All_EKF_HPL,'r');
plot(All_t_save/60,All_LSQ_RAIM_HPL,'g');
grid on;
legend('UKF','EKF','LSQ');

%% plot stanford diagrams
%% figure
HAL = 20;
figure();hold on;
plot(sqrt(All_NEDError_UKF(100:length(All_NEDError_UKF),1).^2+All_NEDError_UKF(100:length(All_NEDError_UKF),2).^2),All_UKF_HPL(100:length(All_NEDError_UKF)),'b.');
plot(sqrt(All_NEDError_EKF(100:length(All_NEDError_EKF),1).^2+All_NEDError_EKF(100:length(All_NEDError_EKF),2).^2),All_EKF_HPL(100:length(All_NEDError_EKF)),'r.');
plot(sqrt(All_NEDError_LSQ(:,1).^2+All_NEDError_LSQ(:,2).^2),All_LSQ_RAIM_HPL,'g.');
%area([HAL 100],[HAL HAL],'FaceColor','r');
plot([0 100],[HAL HAL],'k','LineWidth',2);
plot([HAL HAL],[0 HAL],'k','LineWidth',2);
axis([0 100 0 100]);
text(12,5,'MI');
text(52,5,'HMI');
text(52,35,'MI');
plot([0:3000],[0:3000],'k-','LineWidth',2);
xlabel('Position Error (m)');
ylabel('Protection Level (m)');
legend('UKF','EKF','LSQ');
grid on;
title('Horizontal Stanford Plot');

VAL = 50;
figure();
plot(abs(All_NEDError_UKF(100:length(All_NEDError_UKF),3)),All_UKF_VPL(100:length(All_NEDError_UKF)),'b.');
hold on;
plot(abs(All_NEDError_EKF(100:length(All_NEDError_EKF),3)),All_EKF_VPL(100:length(All_NEDError_EKF)),'r.');
plot(sqrt(All_NEDError_LSQ(:,1).^2+All_NEDError_LSQ(:,2).^2),All_LSQ_RAIM_HPL,'g.');
plot([0 100],[VAL VAL],'k','LineWidth',2);
plot([VAL VAL],[0 VAL],'k','LineWidth',2);
axis([0 100 0 100]);
plot([0:3000],[0:3000],'k-','LineWidth',2);
text(32,15,'MI');
text(72,25,'HMI');
text(72,65,'MI');
xlabel('Position Error (m)');
ylabel('Protection Level (m)');
legend('UKF','EKF','LSQ');
grid on;
title('Vertical Stanford Plot');


%% generate a 2-d frequency plot of error vs pl
PLRange = [0:0.1:100];
NumPLBins = length(PLRange)-1;
for PLBin = 1:NumPLBins
    % get a histogram for this PL range
    i = 1;
    blah(i) = 0;
    for Epoch_lo = 1:length(All_UKF_HPL)
       if (All_UKF_HPL(Epoch_lo) > PLRange(PLBin) && All_UKF_HPL(Epoch_lo) < PLRange(PLBin+1))
         blah(i) = sqrt(All_NEDError_UKF(Epoch_lo,1).^2+All_NEDError_UKF(Epoch_lo,2).^2);
         i = i+1;
       end
    end
    
    histPL(:,PLBin) = hist(blah,PLRange);
    clear blah;
end



figure();
plot(All_t_save/60,abs(All_NEDError_UKF(:,1)),'b'); hold on;
plot(All_t_save/60,2*All_UKF_S_save_NED(:,1),'b--');
plot(All_t_save/60,All_UKF_HPL,'b:');

plot(All_t_save/60,abs(All_NEDError_LSQ(:,1)),'r');
plot(All_t_save/60,All_LSQ_DOP(:,3)*GPS_PR_UERE*2,'r--');
plot(All_t_save/60,All_LSQ_RAIM_HPL,'r:');

plot(All_t_save/60,abs(All_NEDError_EKF(:,1)),'g');
plot(All_t_save/60,2*All_EKF_S_save_NED(:,1),'g--');
plot(All_t_save/60,All_EKF_HPL,'g:');

%axis([All_t_save(1)/60 All_t_save(length(All_t_save))/60 0 556]);
grid on;
title('Static Positioning UKF and LSQ Errors');
legend('UKF Error','UKF 2-\sigma bound','UKF HPL','LSQ Error','LSQ 2-\sigma bound','LSQ HPL','EKF Error','EKF 2-\sigma bound','EKF HPL');
ylabel('North Error (m)');
xlabel('Test Time (mins)');

figure();
plot(All_t_save/60,All_UKF_HPL,'b');
hold on;
plot(All_t_save/60,All_EKF_HPL,'g');
plot(All_t_save/60,All_LSQ_RAIM_HPL,'r');
grid on;
axis([0 1440 0 100]);
xlabel('Time (mins)');
ylabel('Protection Level (m)');
title('Horizontal Protection Level');

figure;
plot(All_t_save/60,All_UKF_VPL,'b');
hold on;
plot(All_t_save/60,All_EKF_VPL,'g');
plot(All_t_save/60,All_LSQ_RAIM_VPL,'r');
grid on;
axis([0 1440 0 100]);
xlabel('Time (mins)');
ylabel('Protection Level (m)');
title('Vertical Protection Level');





