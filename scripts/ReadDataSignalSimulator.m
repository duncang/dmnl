


TruthOffset = 39;
% convert simtruth time (col 1) to gps seconds
for i=1:length(sim_truth)
    hrs_in = floor(sim_truth(i,1)/3600/1000);
    min_in = floor(sim_truth(i,1)/60/1000);
    sec_in = sim_truth(i,1)/1000 - hrs_in*3600 - min_in * 60;
    [gpsweek_out(i), gpssec_out(i)] = ftime([07 03 01 hrs_in min_in sec_in]);
    sim_truth(i,1) = gpssec_out(i);
end

bestpos_2 = GARD_ReadBestPos('data/Simulator_Data/Mar_07/bestpos_1416_2.txt');

% a hack for this file = approx pos is wrong
% ApproxPos = [-5046719.69001578,2568403.35951166,-2925318.76002602];
% meanECEF =  [-5046773.35774802  2568446.08440315 -2925289.01760974];
% 
% ALPHA = [9.3132E-09  0.0000E+00 -5.9605E-08  0.0000E+00];
% BETA = [9.0112D+04  0.0000E+00 -1.9661E+05  0.0000E+00];

ALPHA(1) = 4.6566128999999998e-009 ;
ALPHA(2) = 1.4901161000000001e-008 ;
ALPHA(3) = -5.9604600000000002e-008;
ALPHA(4) = -5.9604600000000002e-008;
BETA(1) = 7.9872000000000000e+004;
BETA(2) = 6.5536000000000000e+004;
BETA(3) = -6.5536000000000000e+004;
BETA(4) = -3.9321600000000000e+005;
            
% load observation and nav data
% [GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ValidDataRinexObs,ApproxPos, Novatel_C1, Novatel_L1, Novatel_D1, Novatel_S1, Novatel_P2, Novatel_L2, Novatel_D2, Novatel_S2] = ReadRinexNovatel('data/Ground_Test_Data/2Feb2007/02020330.07O');
% [SV_Ephemeris, ALPHA, BETA] = freadnav('data/Ground_Test_Data/2Feb2007/02020330.07N');
% 
NumberEpochsGPS = size(Novatel_C1,2);
gps_dt = 1;

NumberSVs = size(Novatel_C1,1);
SVDontUse = zeros(1,NumberSVs);
% 
% % propogate satellite orbits and velocities
% disp('Calculating Satellite Orbits...');
% 
% % preallocate storage vectors for speed
% SV_X_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Y_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Z_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_T_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Xvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Yvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Zvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Tvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Xacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Yacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Zacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Tacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% 
% for i=1:NumberEpochsGPS
%     for j=1:NumberSVs
%         [SV_X_Data(i,j), SV_Y_Data(i,j), SV_Z_Data(i,j), SV_T_Data(i,j), ValidPosData(i,j)] = ...
%                GPSOrbitPropagator(GPSTime_Week(i), GPSTime_Sec(i), j, SV_Ephemeris, 7500);
%         [SV_Xvel_Data(i,j), SV_Yvel_Data(i,j), SV_Zvel_Data(i,j), SV_Tvel_Data(i,j), ...
%          SV_Xacc_Data(i,j), SV_Yacc_Data(i,j), SV_Zacc_Data(i,j), SV_Tacc_Data(i,j), ValidVelData(i,j)] = ...
%          GPSOrbitPropagatorVelocities(GPSTime_Week(i),GPSTime_Sec(i), j, SV_Ephemeris);
%     end
%     
%     if(mod(i,1000) == 0)
%         disp(sprintf('Completed Epoch %d',i));
%     end
% end


%% arrange pseudorange measurements

PR_Sim = Novatel_C1';
PRR_Sim = -Novatel_D1' * L1_Wavelength;
CP_Sim = Novatel_L1';


%Load in Roll Pitch Yaw data and IMU acceleration data

%Interpolate data to get it at 50 Hz. 


for i = 1:endepochReadIn   
    
    
      %this is aerosim output from its airplane block. Note this is the same as
    %XposCalc.
    Xpos_truth1Hz(i) =  sim_truth(i,2);
    Ypos_truth1Hz(i) =  sim_truth(i,3);
    Zpos_truth1Hz(i) =  sim_truth(i,4);
    
    Xvel_truth1Hz(i) =  sim_truth(i,5);
    Yvel_truth1Hz(i) =  sim_truth(i,6);
    Zvel_truth1Hz(i) =  sim_truth(i,7);
    
    Xacc_truth1Hz(i) = sim_truth(i,8);
    Yacc_truth1Hz(i) = sim_truth(i,9);
    Zacc_truth1Hz(i) = sim_truth(i,10);
        
    Xjerk_truth1Hz(i) = sim_truth(i,11);
    Yjerk_truth1Hz(i) = sim_truth(i,12);
    Zjerk_truth1Hz(i) = sim_truth(i,13);            
              
    Lat_truth1Hz(i) = sim_truth(i,14);
    Lon_truth1Hz(i) = sim_truth(i,15);
    Hgt_truth1Hz(i) = sim_truth(i,16);
    
     
    Heading_truth1Hz(i) = sim_truth(i,17);        %these are in local NED
    Elevation_truth1Hz(i) =  sim_truth(i,18);
    Bank_truth1Hz(i) = sim_truth(i,19);
    
    %dunk and lennon:
    %can't transform angles the same way u can do velocity (or position
    %error) going from ECEF to NED. 
    %but should be able to transform angular rates because its a unit
    %vector? 
    %so using transformation matrix should be valid
    
    X_ang_dot_truth1Hz(i) = sim_truth(i,20);   %these are angular rates about ECEF axes
    Y_ang_dot_truth1Hz(i) =  sim_truth(i,21);
    Z_ang_dot_truth1Hz(i) = sim_truth(i,22);
    
                                          
 %In motion_v1.csv

%Ang_vel X etc is euler rates. Not really p ,q, r body rates but will be close enough for what i want to do

%Heading, elevation and bank, is flight path heading elevation and bank. I
%can assume they are aircraft roll pitch and yaw but with zero angle of attack and
%sideslip
              
        
    %antenna data in the motion_v1.csv is there because u can include a
    %delta between the c of g and GPS antenna, in the simulation. 
    %so for this simulation dunc didnt add in a delta so assume the gps
    %antenna is at the c of g.       
    
    %assume that body axes of aircraft aligned with nav frame (i.e. no
    %angle of attack/sideslip, aligned with 
    Roll_truth1Hz(i) = Bank_truth1Hz(i);          
    Pitch_truth1Hz(i) = Elevation_truth1Hz(i);
    Yaw_truth1Hz(i) = Heading_truth1Hz(i);   
    
        
    
     [quat] = EulerToQuat([Roll_truth1Hz(i), Pitch_truth1Hz(i), Yaw_truth1Hz(i)]);
    
     Quaternions_truth1Hz(1,i) = quat(1);
     Quaternions_truth1Hz(2,i) = quat(2);
     Quaternions_truth1Hz(3,i) = quat(3);
     Quaternions_truth1Hz(4,i) = quat(4);
    
               
      
     %===========================
     %  %NED Velocities and accelerations   
     %===========================
      

TMatrix = T_ECEF2NED(Lon_truth1Hz(i),Lat_truth1Hz(i));

NEDangacc_truth1Hz(:,i) = TMatrix*[X_ang_dot_truth1Hz(i),Y_ang_dot_truth1Hz(i),Z_ang_dot_truth1Hz(i)]';
N_ang_dot_truth1Hz(i) = NEDangacc_truth1Hz(1,i);
E_ang_dot_truth1Hz(i) = NEDangacc_truth1Hz(2,i);
D_ang_dot_truth1Hz(i) = NEDangacc_truth1Hz(3,i);


NEDacc_truth1Hz(:,i) = TMatrix*[Xacc_truth1Hz(i),Yacc_truth1Hz(i),Zacc_truth1Hz(i)]';
A_n_truth1Hz(i) = NEDacc_truth1Hz(1,i);
A_e_truth1Hz(i) = NEDacc_truth1Hz(2,i);
A_d_truth1Hz(i) = NEDacc_truth1Hz(3,i);      
    
    
NEDvel_truth1Hz(:,i) = TMatrix*[Xvel_truth1Hz(i),Yvel_truth1Hz(i),Zvel_truth1Hz(i)]';  
V_n_truth1Hz(i) = NEDvel_truth1Hz(1,i);
V_e_truth1Hz(i) = NEDvel_truth1Hz(2,i);
V_d_truth1Hz(i) = NEDvel_truth1Hz(3,i);
 
 
    %Body velocities and accelerations
    
   TMatrix1Hz = T_Body2NED(Roll_truth1Hz(i),Pitch_truth1Hz(i), Yaw_truth1Hz(i));
        
   pqr_truth1Hz(:,i) = TMatrix1Hz'*[N_ang_dot_truth1Hz(i),E_ang_dot_truth1Hz(i),D_ang_dot_truth1Hz(i)]';
     
   p_truth1Hz(i) = pqr_truth1Hz(1,i);   % p
    q_truth1Hz(i) = pqr_truth1Hz(2,i);    %q
    r_truth1Hz(i) = pqr_truth1Hz(3,i);    %r
    
    
    Axyz_truth1Hz(:,i) = TMatrix1Hz'*[A_n_truth1Hz(i),A_e_truth1Hz(i),A_d_truth1Hz(i)]';
    
    A_x_truth1Hz(i) =   Axyz_truth1Hz(1,i);
    A_y_truth1Hz(i) =   Axyz_truth1Hz(2,i);
    A_z_truth1Hz(i) =   Axyz_truth1Hz(3,i);
     
         
    %velocities in body
    
    uvw_truth1Hz = TMatrix1Hz'*[V_n_truth1Hz(i),V_e_truth1Hz(i),V_d_truth1Hz(i)]';
    
    u_truth1Hz(i) = uvw_truth1Hz(1);
       v_truth1Hz(i) = uvw_truth1Hz(2);
          w_truth1Hz(i) = uvw_truth1Hz(3);
    
          
    
    
    %ABOVE THIS IS fromt eh SIMULATOR
    
    
    
    
    
    %Approximate INS outputs (adding noise).      
       
        
        
      %velocity truth

%     V_n_truth(i) = vel_truth(2,i);
%       V_e_truth(i) = vel_truth(3,i);
%         V_d_truth(i) = vel_truth(4,i);
%               
%         
%         Xvel_truth(i) =   vel_truthECEF(2,i);
%         Yvel_truth(i) =   vel_truthECEF(3,i);
%         Zvel_truth(i) =   vel_truthECEF(4,i);
%         

    %sensor measuremements
    %Read in INS sensor data
           
  
        

    %Attitude truth

    Roll_truth(i) = att_truth(2,i);   % bank angle(rad)  PHI
    Pitch_truth(i) = att_truth(3,i);    %pitch angle (rad) THETA
    Yaw_truth(i) = att_truth(4,i);    %heading (rad) PSI


    %Angular Accelerations

    p_truth(i) = AngRate(2,i);   % p
    q_truth(i) = AngRate(3,i);    %q
    r_truth(i) = AngRate(4,i);    %r


    p_dot_truth(i) = AngAcc(2,i);   % p dot
    q_dot_truth(i) = AngAcc(3,i);    %q dot
    r_dot_truth(i) = AngAcc(4,i);    %r dot


    %Air data
    %CORRECT

    Airspeed_truth(i) = AirData1(2,i);    %Airspeed (m/s)
    Beta_truth(i) = AirData1(3,i);     %Sideslip (rad)
    Alpha_truth(i) = AirData1(4,i);              %AOA(rad)

    % flight controls
    %CORRECT

    Flap_truth(i) = flight_controls(2,i);          % Flap (rad)
    Elevator_truth(i) = flight_controls(3,i);       %Elevator (rad)
    Aileron_truth(i) = flight_controls(4,i);        % Aileron (rad)
    Rudder_truth(i) = flight_controls(5,i);         % Rudder (rad)
    Throttle_truth(i) = flight_controls(6,i);          %Throttle (frac (0 to 1??))
    Mixture_truth(i) = flight_controls(7,i);       %Mixture (ratio)
    Ignition_truth(i) = flight_controls(8,i);      %Ignition (bool)


    %coefficients aerosim calculates
    %aerodynamic coeffs.

    CD_truth(i) = AeroCoeff(2,i);    %rad^-1
    CY_truth(i) = AeroCoeff(3,i);     %ditto
    CL_truth(i) = AeroCoeff(4,i);
    Cl_truth(i) = AeroCoeff(5,i);
    Cm_truth(i) = AeroCoeff(6,i);
    Cn_truth(i) = AeroCoeff(7,i);


    %prop coeffs.

    J_truth(i) = PropCoeff(2,i);   % J
    CT_truth(i) = PropCoeff(3,i);    %CT
    CP_truth(i) = PropCoeff(4,i);    %CP

    %Eng. coeffs.

    MAP_truthEng(i) = EngCoeff(2,i);   % MAP  Current manifold pressure for current throttle setting and altitude, in kPa
    m_dotAir_truthEng(i) = EngCoeff(3,i);    %m_dotAir current instantaneous mass air flow
    m_dotFuel_truthEng(i) = EngCoeff(4,i);    %m_dotFuel current instantaneous mass fuel flow
    BSFC_truthEng(i) = EngCoeff(5,i);    %BSFC  Brake specific fuel consumption
    P_truthEng(i) = EngCoeff(6,i);    %P  Current engine power


    OMEGA_truthEng(i) = EngParameters(2,i);

    Faero_truth(1,i) = Forces(2,i);   % all forces In body axes
    Faero_truth(2,i) = Forces(3,i);   % all forces In body axes
    Faero_truth(3,i) = Forces(4,i);   % all forces In body axes

    Fprop_truth(1,i) = Forces(5,i);    %
    Fprop_truth(2,i) = Forces(6,i);    %
    Fprop_truth(3,i) = Forces(7,i);    %


    Mass_truth(i) = Forces(8,i);    %mass used in total force calculation

    Maero_truth(1,i) = Moments(2,i);   % all moments are in body axes
    Maero_truth(2,i) = Moments(3,i);   % all moments are in body axes
    Maero_truth(3,i) = Moments(4,i);   % all moments are in body axes


    Mprop_truth(1,i) = Moments(5,i);    %
    Mprop_truth(2,i) = Moments(6,i);    %
    Mprop_truth(3,i) = Moments(7,i);    %

    Mcg_truth(1,i) = Moments(8,i);    %Moment about the c of g
    Mcg_truth(2,i) = Moments(9,i);    %Moment about the c of g
    Mcg_truth(3,i) = Moments(10,i);    %Moment about the c of g


    Jx_truth(i) = InertiaCG(2,i);
    Jy_truth(i)= InertiaCG(3,i);
    Jz_truth(i)= InertiaCG(4,i);
    Jxz_truth(i) = InertiaCG(5,i);
    CGxpos_truth(i) = InertiaCG(6,i);
    CGypos_truth(i) = InertiaCG(7,i);
    CGzpos_truth(i) = InertiaCG(8,i);

    AtmosTruth(1,i)  = AtmosGrav(2,i);   %this is P
    AtmosTruth(2,i)  = AtmosGrav(3,i);                          %this is T
    AtmosTruth(3,i)  = AtmosGrav(4,i);                 %this is rho
    AtmosTruth(4,i)  = AtmosGrav(5,i);              %this is a (speed of sound)
    GravityTruth(i) = AtmosGrav(6,i);

    Mach_truth(i) = Mach(2,i);

    %Quaternions

    Quaternions_truth(1,i) = Quaternions(2,i); %q0
    Quaternions_truth(2,i) = Quaternions(3,i); %q1
    Quaternions_truth(3,i) = Quaternions(4,i); %q2
    Quaternions_truth(4,i) = Quaternions(5,i); %q3


    %Load wind speeds
    u_wind_truth(i) = WindB(2,i);
    v_wind_truth(i) = WindB(3,i);
    w_wind_truth(i) = WindB(4,i);

    p_wind_truth(i) = WindRates(2,i);
    q_wind_truth(i) = WindRates(3,i);
    r_wind_truth(i) = WindRates(4,i);
    
    A_x_INS1Hz(i) = A_x_truth1Hz(i);
A_y_INS1Hz(i) =  A_y_truth1Hz(i);
A_z_INS1Hz(i) =  A_z_truth1Hz(i);

p_INS1Hz(i) = p_truth1Hz(i);
q_INS1Hz(i) = q_truth1Hz(i); 
r_INS1Hz(i) = r_truth1Hz(i);
    
    
end




%differentiate velocities to get acceleration. To check that it's correct. 

XaccCalc = diff(Xvel_truth1Hz);
YaccCalc = diff(Yvel_truth1Hz);
ZaccCalc = diff(Zvel_truth1Hz);


%Use state 0 of the random number generator

randn('state',0);

%interpolate acceleration data to
%get it from 1 Hz to 50 Hz

%Rate
%Frequency
Rate = 50; %Hertz INS output rate



f = 1/Rate;

x1 = 0:length(A_x_truth1Hz)-1;
y1 = A_x_truth1Hz;
xx1 = 0:f:length(A_x_truth1Hz)-1;  %1/50 = 0.02
A_x_truth50Hz = spline(x1,y1,xx1);

x2 = 0:length(A_y_truth1Hz)-1;
y2 = A_y_truth1Hz;
xx2 = 0:f:length(A_y_truth1Hz)-1;  %1/50 = 0.02
A_y_truth50Hz = spline(x1,y1,xx2);

x3 = 0:length(A_z_truth1Hz)-1;
y3 = A_z_truth1Hz;
xx3 = 0:f:length(A_z_truth1Hz)-1;  %1/50 = 0.02
A_z_truth50Hz = spline(x3,y3,xx3);

x4 = 0:length(p_truth1Hz)-1;
y4 = p_truth1Hz;
xx4 = 0:f:length(p_truth1Hz)-1;  %1/50 = 0.02
p_truth50Hz = spline(x4,y4,xx4);

x5 = 0:length(q_truth1Hz)-1;
y5 = q_truth1Hz;
xx5 = 0:f:length(q_truth1Hz)-1;  %1/50 = 0.02
q_truth50Hz = spline(x5,y5,xx5);


x6 = 0:length(r_truth1Hz)-1;
y6 = r_truth1Hz;
xx6 = 0:f:length(r_truth1Hz)-1;  %1/50 = 0.02
r_truth50Hz = spline(x6,y6,xx6);




%Simulate INS error with Gauss Markov process


Sigma_Xacc_noise = 0.0001; %m/s^2  %Power spectral density of X accelerometer
Sigma_Yacc_noise = 0.0001; %m/s^2  %Power spectral density of X accelerometer
Sigma_Zacc_noise = 0.0001; %m/s^2  %Power spectral density of X accelerometer

Sigma_p_gyro_noise = 0.0004848; %rad/s  %Power spectral density of X accelerometer  %100 deg/hr gyros
Sigma_q_gyro_noise = 0.0004848; %rad/s  %Power spectral density of X accelerometer
Sigma_r_gyro_noise = 0.0004848; %rad/s  %Power spectral density of X accelerometer


Beta_Xacc_noise = 300;
Beta_Yacc_noise = 300;
Beta_Zacc_noise = 300;

Beta_p_gyro_noise = 300;
Beta_q_gyro_noise = 300;
Beta_r_gyro_noise = 300;


Xacc_noise(1) = Sigma_Xacc_noise*randn(1); 
Yacc_noise(1) = Sigma_Yacc_noise*randn(1); 
Zacc_noise(1) = Sigma_Zacc_noise*randn(1); 

p_gyro_noise(i) = Sigma_p_gyro_noise*randn(1); 
q_gyro_noise(i) = Sigma_q_gyro_noise*randn(1); 
r_gyro_noise(i) = Sigma_r_gyro_noise*randn(1); 

for i = 2:length(A_x_truth50Hz)  
    
    [Xacc_noise(i)] = GaussMarkov_Process(Xacc_noise(i-1),Beta_Xacc_noise, Sigma_Xacc_noise);
      [Yacc_noise(i)] = GaussMarkov_Process(Yacc_noise(i-1),Beta_Yacc_noise, Sigma_Yacc_noise);
        [Zacc_noise(i)] = GaussMarkov_Process(Zacc_noise(i-1),Beta_Zacc_noise, Sigma_Zacc_noise);
        
          [p_gyro_noise(i)] = GaussMarkov_Process(p_gyro_noise(i-1),Beta_p_gyro_noise,Sigma_p_gyro_noise);
           [q_gyro_noise(i)] = GaussMarkov_Process(q_gyro_noise(i-1),Beta_q_gyro_noise,Sigma_q_gyro_noise);
           [r_gyro_noise(i)] = GaussMarkov_Process(r_gyro_noise(i-1),Beta_r_gyro_noise,Sigma_r_gyro_noise);
        
              
end




%simulated INS output
for i = 1:length(A_x_truth50Hz)  
    
    
A_x_INS50Hz(i) = A_x_truth50Hz(i); %+ Xacc_noise(i);
A_y_INS50Hz(i) =  A_y_truth50Hz(i); %+ Yacc_noise(i);
A_z_INS50Hz(i) =  A_z_truth50Hz(i);% + Zacc_noise(i);

p_INS50Hz(i) = p_truth50Hz(i);% + p_gyro_noise(i);
q_INS50Hz(i) =  q_truth50Hz(i); %+ q_gyro_noise(i);
r_INS50Hz(i) =  r_truth50Hz(i); %+ r_gyro_noise(i);

end









%     ax_b_INS_50Hz(i) = sensors(2,i)+ l_ax_INS_Error;  %m/s
%     ay_b_INS_50Hz(i) = sensors(3,i)+ l_ay_INS_Error;
%     az_b_INS_50Hz(i) = sensors(4,i)+ l_az_INS_Error;
% 
%     p_INS_50Hz(i) = sensors(5,i)+ l_p_INS_Error;  %rad/s
%     q_INS_50Hz(i) = sensors(6,i)+ l_q_INS_Error;
%     r_INS_50Hz(i) = sensors(7,i)+ l_r_INS_Error;


















