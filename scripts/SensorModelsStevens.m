
clear all
%function SensorModels

%Load data from aerosim 'truth'. 

%With GPS need to assume that the GPS antenna is at the centre of gravity of the airplane. Incorporate GPS antenna location from cofg in body axes later. 
    %  load the truth data
    
 %load ecef positions I calculated from lat lon and hgt
 
%  
%  
% if ~exist('pos_truth_llh')
%     load 'data\BrisbaneFlightWithWind12.10.06\pos_truth_llh.mat';
% end
% 
% %load ecef positions I calculated from lat lon and hgt
% if ~exist('pos_truth_ecef')
%     load 'data\BrisbaneFlightWithWind12.10.06\pos_truth_ecef.mat';
% end
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\BrisbaneFlightWithWind12.10.06\ECEF.mat';
% end
% 
% 
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\BrisbaneFlightWithWind12.10.06\att_truth.mat';
% end
% 
% %load air data
% if ~exist('AirData1')
%     load 'data\BrisbaneFlightWithWind12.10.06\AirData1.mat';
% end
% 
% %load flight controls
% if ~exist('flight_controls')
%     load 'data\BrisbaneFlightWithWind12.10.06\flight_controls.mat';
% end
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\BrisbaneFlightWithWind12.10.06\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\BrisbaneFlightWithWind12.10.06\AngRate.mat';
% end
% 
% %load engine coefficients
% if ~exist('EngCoeff')
%     load 'data\BrisbaneFlightWithWind12.10.06\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\BrisbaneFlightWithWind12.10.06\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\BrisbaneFlightWithWind12.10.06\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\BrisbaneFlightWithWind12.10.06\AeroCoeff.mat';
% end
% 
% %load forces
% if ~exist('Forces')
%     load 'data\BrisbaneFlightWithWind12.10.06\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\BrisbaneFlightWithWind12.10.06\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\BrisbaneFlightWithWind12.10.06\InertiaCG.mat';
% end
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\BrisbaneFlightWithWind12.10.06\AtmosGrav.mat';
% end
% 
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\BrisbaneFlightWithWind12.10.06\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\BrisbaneFlightWithWind12.10.06\Quaternions.mat';
% end
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     load('data\BrisbaneFlightWithWind12.10.06\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\BrisbaneFlightWithWind12.10.06\sensors_clean.mat');
%     sensors = sensors_clean;
%     clear 'sensors_clean';
% end
% 
% 
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\BrisbaneFlightWithWind12.10.06\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\BrisbaneFlightWithWind12.10.06\WindRates.mat';
% end
%  
    
    








% 
% 
% 
% 
% 
%  
% if ~exist('pos_truth_llh')
%     load 'data\Test3\pos_truth_llh.mat';
% end
% 
% %load ecef positions I calculated from lat lon and hgt
% if ~exist('pos_truth_ecef')
%     load 'data\Test3\pos_truth_ecef.mat';
% end
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Test3\ECEF.mat';
% end
% 
% 
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Test3\att_truth.mat';
% end
% 
% %load air data
% if ~exist('AirData1')
%     load 'data\Test3\AirData1.mat';
% end
% 
% %load flight controls
% if ~exist('flight_controls')
%     load 'data\Test3\flight_controls.mat';
% end
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Test3\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Test3\AngRate.mat';
% end
% 
% %load engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Test3\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Test3\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Test3\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Test3\AeroCoeff.mat';
% end
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Test3\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Test3\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Test3\InertiaCG.mat';
% end
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\Test3\AtmosGrav.mat';
% end
% 
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Test3\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Test3\Quaternions.mat';
% end
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     load('data\Test3\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Test3\sensors_clean.mat');
%     sensors = sensors_clean;
%     clear 'sensors_clean';
% end
% 
% 
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Test3\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Test3\WindRates.mat';
% end
%  
% 
% 
% 

%sTRAIGHT IN dESCEND
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Even Descent 2.5.07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Even Descent 2.5.07\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Even Descent 2.5.07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Even Descent 2.5.07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Even Descent 2.5.07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Even Descent 2.5.07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Even Descent 2.5.07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Even Descent 2.5.07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Even Descent 2.5.07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Even Descent 2.5.07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Even Descent 2.5.07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Even Descent 2.5.07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Even Descent 2.5.07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Even Descent 2.5.07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Even Descent 2.5.07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Even Descent 2.5.07\InertiaCG.mat';
% end
% 
% 
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     load('data\Even Descent 2.5.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Even Descent 2.5.07\sensors_clean.mat');
%     sensors = sensors_clean;
%     clear 'sensors_clean';
% end
% 
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\Even Descent 2.5.07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Even Descent 2.5.07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Even Descent 2.5.07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Even Descent 2.5.07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Even Descent 2.5.07\WindRates.mat';
% end

    
%     
%     
%  
% %RNA
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\TestFlight4.7.07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\TestFlight4.7.07\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\TestFlight4.7.07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\TestFlight4.7.07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\TestFlight4.7.07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\TestFlight4.7.07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\TestFlight4.7.07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\TestFlight4.7.07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\TestFlight4.7.07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\TestFlight4.7.07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\TestFlight4.7.07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\TestFlight4.7.07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\TestFlight4.7.07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\TestFlight4.7.07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\TestFlight4.7.07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\TestFlight4.7.07\InertiaCG.mat';
% end
% 
% 
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     load('data\TestFlight4.7.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\TestFlight4.7.07\sensors_clean.mat');
%     sensors = sensors_clean;
%     clear 'sensors_clean';
% end
% 
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\TestFlight4.7.07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\TestFlight4.7.07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\TestFlight4.7.07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\TestFlight4.7.07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\TestFlight4.7.07\WindRates.mat';
% end
% 
% 
%     
    




%RNAV

%POSITION

%load ecef positions aerosim calculates
if ~exist('ECEF')
    load 'data\Flight22Nov07\ECEF.mat';
end


if ~exist('pos_truth_llh')
    load 'data\Flight22Nov07\pos_truth_llh.mat';
end


%VELOCITY
if ~exist('vel_truth')
    load 'data\Flight22Nov07\vel_truth.mat';
end

if ~exist('vel_truthECEF')
    load 'data\Flight22Nov07\vel_truthECEF.mat';
end




%ATTITUDE
%load attitude truth
if ~exist('att_truth')
    load 'data\Flight22Nov07\att_truth.mat';
end


%load angular accelerations ( pdot qdot rdot)
if ~exist('AngAcc')
    load 'data\Flight22Nov07\AngAcc.mat';
end

%load angular rates ( p q r)
if ~exist('AngRate')
    load 'data\Flight22Nov07\AngRate.mat';
end


%AIR DATA
if ~exist('AirData1')
    load 'data\Flight22Nov07\AirData1.mat';
end



%FLIGHT CONTROLS
if ~exist('flight_controls')
    load 'data\Flight22Nov07\flight_controls.mat';
end


%ENGINE coefficients
%engine coefficients
if ~exist('EngCoeff')
    load 'data\Flight22Nov07\EngCoeff.mat';
end

%load extra engine coefficients like shaft rotation speed (related
%to RPM)
if ~exist('EngParameters')
    load 'data\Flight22Nov07\EngParameters.mat';
end

%load prop coefficients
if ~exist('PropCoeff')
    load 'data\Flight22Nov07\PropCoeff.mat';
end

%load aero coefficients
if ~exist('AeroCoeff')
    load 'data\Flight22Nov07\AeroCoeff.mat';
end




%LOAD AERODYNAMIC TRUTH

%load forces
if ~exist('Forces')
    load 'data\Flight22Nov07\Forces.mat';
end

%load moments
if ~exist('Moments')
    load 'data\Flight22Nov07\Moments.mat';
end


%load Inertia and CofG Truth
if ~exist('InertiaCG')
    load 'data\Flight22Nov07\InertiaCG.mat';
end




%Load Inertial sensor data

% load IMU measurements data
%this is from the 'sensors' output of the navion 6 DOF model

load('data\Flight22Nov07\sensors_noisy.mat');
load('data\Flight22Nov07\sensors_clean.mat');

UseNoisy = 0;

if UseNoisy == 1
    
    sensors = sensors_noisy;
    
    %clear 'sensors_noisy';
else
    
    sensors = sensors_clean;
    %clear 'sensors_clean';
end
%get copy of noisy sensor
sensorsNoisy = sensors_noisy;
sensorsClean = sensors_clean;



%ATMOSPHERE TRUTH

%load Atmosphere and Gravity Truth
if ~exist('AtmosGrav')
    load 'data\Flight22Nov07\AtmosGrav.mat';
end

%load Mach number Truth
if ~exist('Mach')
    load 'data\Flight22Nov07\Mach.mat';
end

%load Quaternions Truth
if ~exist('Quaternions')
    load 'data\Flight22Nov07\Quaternions.mat';
end

%Load Wind speed data

%load Wind speeds in body axes (m/s)
if ~exist('WindB')
    load 'data\Flight22Nov07\WindB.mat';
end

%load Wind on angular rates (p q r) in rad/s
if ~exist('WindRates')
    load 'data\Flight22Nov07\WindRates.mat';
end


    
            
    

    
 %simulate rate in hz , can't be greater than 50, must be a whole number
 Rate = 50; 
    
 a = 50/Rate;    
 
 
 %startepoch = 10*Rate;  %start at 10 seconds
 
 startepoch = 1;  %start epoch is 10 for 1 Hz Rate and 19 for 2 Hz
 %data is at 50 Hz so take every 50th for 1 second  
endepoch = 2730;
% endepoch = 386;
%if you want at 1 hz a = 50/1 = 50, therefore you use i*a - (a-1), ie i*50-49
%if you want at 5 hz you do a = 50/5 = 10 therefore you use i*10-9
%for 2 hz you use a = 50/2 = 25, therefore you use i*25-24
%etc
step_size = 1/Rate;  

 


    
for i = startepoch:endepoch    
           
  %sensor measuremements
    %Read in INS sensor data
        
%     ax_b_INS(i) = sensors(2,i*a-(a-1));  %m/s
%     ay_b_INS(i) = sensors(3,i*a-(a-1));
%     az_b_INS(i) = sensors(4,i*a-(a-1));
%     
%     p_INS(i) = sensors(5,i*a-(a-1));  %rad/s
%     q_INS(i) = sensors(6,i*a-(a-1));
%     r_INS(i) = sensors(7,i*a-(a-1));

        
    %read in the data and perform necessary conversions
    %this should be the GPS inputs but use the truth for now.
    %ECEF positions
    %Note that first row is simulation time data thats why the index starts at
    %2 in the following

    Lat_truth(i) = pos_truth_llh(2,i*a-(a-1));
    Lon_truth(i) = pos_truth_llh(3,i*a-(a-1));
    Hgt_truth(i) = pos_truth_llh(4,i*a-(a-1));


%     XposCalc(i) = pos_truth_ecef(2,i*a-(a-1));
%     YposCalc(i) = pos_truth_ecef(3,i*a-(a-1));
%     ZposCalc(i) = pos_truth_ecef(4,i*a-(a-1));


    %this is aerosim output from its airplane block. Note this is the same as
    %XposCalc.
    Xpos_truth(i) = ECEF(2,i*a-(a-1));
    Ypos_truth(i) = ECEF(3,i*a-(a-1));
    Zpos_truth(i) = ECEF(4,i*a-(a-1));


    %Attitude truth

    Roll_truth(i) = att_truth(2,i*a-(a-1));   % bank angle(rad)  PHI
    Pitch_truth(i) = att_truth(3,i*a-(a-1));    %pitch angle (rad) THETA
    Yaw_truth(i) = att_truth(4,i*a-(a-1));    %heading (rad) PSI


    %Angular Accelerations

    p_truth(i) = AngRate(2,i*a-(a-1));   % p
    q_truth(i) = AngRate(3,i*a-(a-1));    %q
    r_truth(i) = AngRate(4,i*a-(a-1));    %r


    p_dot_truth(i) = AngAcc(2,i*a-(a-1));   % p dot
    q_dot_truth(i) = AngAcc(3,i*a-(a-1));    %q dot
    r_dot_truth(i) = AngAcc(4,i*a-(a-1));    %r dot

    

    %Air data
    %CORRECT

    Airspeed_truth(i) = AirData1(2,i*a-(a-1));    %Airspeed (m/s)
    Beta_truth(i) = AirData1(3,i*a-(a-1));     %Sideslip (rad)
    Alpha_truth(i) = AirData1(4,i*a-(a-1));              %AOA(rad)

    % flight controls
    %CORRECT

    Flap_truth(i) = flight_controls(2,i*a-(a-1));          % Flap (rad)
    Elevator_truth(i) = flight_controls(3,i*a-(a-1));       %Elevator (rad)
    Aileron_truth(i) = flight_controls(4,i*a-(a-1));        % Aileron (rad)
    Rudder_truth(i) = flight_controls(5,i*a-(a-1));         % Rudder (rad)
    Throttle_truth(i) = flight_controls(6,i*a-(a-1));          %Throttle (frac (0 to 1??))
    Mixture_truth(i) = flight_controls(7,i*a-(a-1));       %Mixture (ratio)
    Ignition_truth(i) = flight_controls(8,i*a-(a-1));      %Ignition (bool)


    
    %coefficients aerosim calculates
    %aerodynamic coeffs.

    CD_truth(i) = AeroCoeff(2,i*a-(a-1));    %rad^-1
    CY_truth(i) = AeroCoeff(3,i*a-(a-1));     %ditto
    CL_truth(i) = AeroCoeff(4,i*a-(a-1));
    Cl_truth(i) = AeroCoeff(5,i*a-(a-1));
    Cm_truth(i) = AeroCoeff(6,i*a-(a-1));
    Cn_truth(i) = AeroCoeff(7,i*a-(a-1));


    %prop coeffs.

    J_truth(i) = PropCoeff(2,i*a-(a-1));   % J
    CT_truth(i) = PropCoeff(3,i*a-(a-1));    %CT
    CP_truth(i) = PropCoeff(4,i*a-(a-1));    %CP

    %Eng. coeffs.

    MAP_truthEng(i) = EngCoeff(2,i*a-(a-1));   % MAP  Current manifold pressure for current throttle setting and altitude, in kPa
    m_dotAir_truthEng(i) = EngCoeff(3,i*a-(a-1));    %m_dotAir current instantaneous mass air flow
    m_dotFuel_truthEng(i) = EngCoeff(4,i*a-(a-1));    %m_dotFuel current instantaneous mass fuel flow
    BSFC_truthEng(i) = EngCoeff(5,i*a-(a-1));    %BSFC  Brake specific fuel consumption
    P_truthEng(i) = EngCoeff(6,i*a-(a-1));    %P  Current engine power


    OMEGA_truthEng(i) = EngParameters(2,i*a-(a-1));

    Faero_truth(1,i) = Forces(2,i*a-(a-1));   % all forces In body axes
    Faero_truth(2,i) = Forces(3,i*a-(a-1));   % all forces In body axes
    Faero_truth(3,i) = Forces(4,i*a-(a-1));   % all forces In body axes

    Fprop_truth(1,i) = Forces(5,i*a-(a-1));    %
    Fprop_truth(2,i) = Forces(6,i*a-(a-1));    %
    Fprop_truth(3,i) = Forces(7,i*a-(a-1));    %


    Mass_truth(i) = Forces(8,i*a-(a-1));    %mass used in total force calculation

    Maero_truth(1,i) = Moments(2,i*a-(a-1));   % all moments are in body axes
    Maero_truth(2,i) = Moments(3,i*a-(a-1));   % all moments are in body axes
    Maero_truth(3,i) = Moments(4,i*a-(a-1));   % all moments are in body axes


    Mprop_truth(1,i) = Moments(5,i*a-(a-1));    %
    Mprop_truth(2,i) = Moments(6,i*a-(a-1));    %
    Mprop_truth(3,i) = Moments(7,i*a-(a-1));    %

    Mcg_truth(1,i) = Moments(8,i*a-(a-1));    %Moment about the c of g
    Mcg_truth(2,i) = Moments(9,i*a-(a-1));    %Moment about the c of g
    Mcg_truth(3,i) = Moments(10,i*a-(a-1));    %Moment about the c of g


    Jx_truth(i) = InertiaCG(2,i*a-(a-1));
    Jy_truth(i)= InertiaCG(3,i*a-(a-1));
    Jz_truth(i)= InertiaCG(4,i*a-(a-1));
    Jxz_truth(i) = InertiaCG(5,i*a-(a-1));
    CGxpos_truth(i) = InertiaCG(6,i*a-(a-1));
    CGypos_truth(i) = InertiaCG(7,i*a-(a-1));
    CGzpos_truth(i) = InertiaCG(8,i*a-(a-1));


    
    AtmosTruth(1,i)  = AtmosGrav(2,i*a-(a-1));   %this is P
    AtmosTruth(2,i)  = AtmosGrav(3,i*a-(a-1));                          %this is T
    AtmosTruth(3,i)  = AtmosGrav(4,i*a-(a-1));                 %this is rho
    AtmosTruth(4,i)  = AtmosGrav(5,i*a-(a-1));              %this is a (speed of sound)
    GravityTruth(i) =  AtmosGrav(6,i*a-(a-1));

    Mach_truth(i) = Mach(2,i*a-(a-1));

    %Quaternions

    Quaternions_truth(1,i) = Quaternions(2,i*a-(a-1)); %q0
    Quaternions_truth(2,i) = Quaternions(3,i*a-(a-1)); %q1
    Quaternions_truth(3,i) = Quaternions(4,i*a-(a-1)); %q2
    Quaternions_truth(4,i) = Quaternions(5,i*a-(a-1)); %q3    
       
        
%     %Load wind speeds
%     u_wind_truth(i) = WindB(2,i*a-(a-1));
%      v_wind_truth(i) = WindB(3,i*a-(a-1));
%       w_wind_truth(i) = WindB(4,i*a-(a-1));
%             
%         p_wind_truth(i) = WindRates(2,i*a-(a-1));
%         q_wind_truth(i) = WindRates(3,i*a-(a-1));
%         r_wind_truth(i) = WindRates(4,i*a-(a-1));     
%     
  

end;




%=====================================================================
%Gravity Model
%=====================================================================
% Calculate Earth Parameters (WGS-84)

%calculate gravity based on true position




%[Gravity] = Earth_Gravity(x_LLH(1:3)',0.0818191908426^2,9.7803267714,0.00193185138639);

%just advance it by 1 save me getting rid of all the -1's in the code below
%startepoch = startepoch+1;

Gravity = GravityTruth(startepoch);


%=====================================================================
%Atmospheric Parameters
%=====================================================================

%Standard Atmosphere in SI units

%Altitude % ASsume GA aircraft operating betweeen 0 and 15,000 feet. 

%Temperature




%[P_atm T_atm p_atm a_atm] = Atmosphere_STDAtmosphere(x_LLH(3),287.0531,1.40,9.80665);


P_atm = AtmosTruth(1,startepoch); 

T_atm = AtmosTruth(2,startepoch) ;

p_atm = AtmosTruth(3,startepoch) ;

a_atm = AtmosTruth(4,startepoch) ;




%absolute ambient temperature in Kelvins
T_amb = T_atm; %K
R_gas = 287.0531,; %m^2/(K.s^2)

%Pressure

%Density
rho = p_atm ; %kg/m^3 %density




Atmos = [T_amb R_gas rho];   %just use these parameters for now. 



%==============================================
%Convert from Fixed frame (ECEF) to body frame for the following:
%=================================================================


% vel_x = (Xpos(startepoch-1) - Xpos(startepoch-2))/step_size;  %m/s  %differentiate position to get velocities, need to divide by step_size to make it into m/s, otherwise it is m/stepsizeth second. 
% vel_y = (Ypos(startepoch-1) - Ypos(startepoch-2))/step_size; 
% vel_z = (Zpos(startepoch-1) - Zpos(startepoch-2))/step_size; 
% 

%convert ECEF velocities to NED. 


% [Lat,Lon,Hgt] = ECEF2LLH([Xpos(startepoch-2), Ypos(startepoch-2), Zpos(startepoch-2)]);   %I assume you use epoch 8 and not at epoch 9 as coordinate origin for the NED. 
% TMatrixECEF2NED = T_ECEF2NED(Lon,Lat);
% vel_NED = TMatrixECEF2NED*[vel_x, vel_y, vel_z]';
% %convert NED velocities to body 
% TMatrixBody2NED = T_Body2NED( Bank(startepoch-2),Pitch(startepoch-2), Heading(startepoch-2)); %phi, theta, psi ie roll pitch yaw
% vel_uvw = TMatrixBody2NED'*vel_NED; %remember use the transform of TMatrixBody2NED since we want NED 2 Body.



% EulerRates = [ (Bank(startepoch-1) - Bank(startepoch-2))/step_size ; (Pitch(startepoch-1)- Pitch(startepoch-2))/step_size; (Heading(startepoch-1)-Heading(startepoch-2))/step_size] %psi_dot, theta_dot, phi_dot

%convert from euler rates to angular velocities in the body frame (p q r). 
% phi =  Bank(startepoch);
% theta =  Pitch(startepoch);
 
%calculated pqr's
%pqrCalc = T_EulerRate2pqrbody(phi,theta)*EulerRates; 
 

%use truth for pqr
pqr = [p_truth(startepoch), q_truth(startepoch), r_truth(startepoch)];


%true airspeed 
V_T0 = Airspeed_truth(startepoch);  %I think Airspeed is the true airspeed - better check this
beta0 = Beta_truth(startepoch);   %sideslip angle
alpha0 = Alpha_truth(startepoch);
p0 =  pqr(1);      %roll axis angular rate in body frame
q0 = pqr(2);   %pitch angular rate in the body frame 
r0 =  pqr(3);      %yaw axis angular rate in body frame
phi0 = Roll_truth(startepoch);    %bank angle euler angle
theta0 = Pitch_truth(startepoch);  %pitch angle euler angle
psi0 = Yaw_truth(startepoch);    %yaw euler angle


%add 
%starting position in NED frame

% [Lat,Lon,Hgt] = ECEF2LLH([Xpos(startepoch), Ypos(startepoch), Zpos(startepoch)]);   %I assume you use epoch 8 and not at epoch 9 as coordinate origin for the NED. 
% TMatrixECEF2NED = T_ECEF2NED(Lon,Lat);
% 
% NEDTruth = TMatrixECEF2NED*[Xpos(startepoch), Ypos(startepoch), Zpos(startepoch)]';

%is this right??
% pn0 = NEDTruth(1);
% pe0 = NEDTruth(2);
% h0 = NEDTruth(3); 
pn0 = 0;  %just set them to 0 because they are not being used in the model anyway
pe0 = 0;
h0 = 0;
%calculate initial quaternions from euler angles

% Cphi2 = cos(phi0/2);
% Ctheta2 = cos(theta0/2);
% Cpsi2 = cos(psi0/2);
% 
% Sphi2 = sin(phi0/2);
% Stheta2 = sin(theta0/2);
% Spsi2 = sin(psi0/2);
% 
% q00 = +(Cphi2*Ctheta2*Cpsi2 + Sphi2*Stheta2*Spsi2);
% q10 = +(Sphi2*Ctheta2*Cpsi2 - Cphi2*Stheta2*Spsi2);
% q20 = +(Cphi2*Stheta2*Cpsi2 + Sphi2*Ctheta2*Spsi2);
% q30 = +(Cphi2*Ctheta2*Spsi2 - Sphi2*Stheta2*Cpsi2);


quat = EulerToQuat([phi0, theta0, psi0]);

q00 = quat(1);
q10 = quat(2);
q20 = quat(3);
q30 = quat(4);


Lat0 = Lat_truth(startepoch);
Lon0 = Lon_truth(startepoch);
Hgt0 = Hgt_truth(startepoch);


% u0 = vel_uvw(1);
% v0 = vel_uvw(2);
% w0 = vel_uvw(3);


%x0(startepoch-1,:) = [V_T0, beta0 ,alpha0, phi0, theta0, psi0, 0,0,0,pn0, pe0, h0, q00, q10, q20, q30];

%x0(startepoch-1,:) = [V_T0, beta0 ,alpha0, 0, 0, 0, p0,q0,r0,pn0, pe0, h0, q00, q10, q20, q30];
%x0(startepoch-1,:) = [V_T0, 0 ,0, 0, 0, 0, 0,0,0,pn0, pe0, h0, q00, q10, q20, q30];



%x00 = x0(startepoch,:);


%a_sos = a_atm;


M0 = V_T0/a_atm; 


%starting point for the model
XposModel(startepoch) = Xpos_truth(startepoch);
YposModel(startepoch) = Ypos_truth(startepoch);
ZposModel(startepoch) = Zpos_truth(startepoch);

ts = step_size;
%ts= 0.02;


AeroCoeffsTruth = [CD_truth(startepoch), CY_truth(startepoch), CL_truth(startepoch), Cl_truth(startepoch), Cm_truth(startepoch), Cn_truth(startepoch)];

EngCoeffsTruth = [MAP_truthEng(startepoch), m_dotAir_truthEng(startepoch),m_dotFuel_truthEng(startepoch),BSFC_truthEng(startepoch),P_truthEng(startepoch),  OMEGA_truthEng(startepoch)] ;

PropCoeffsTruth = [J_truth(startepoch), CT_truth(startepoch) ,CP_truth(startepoch)];


x0(startepoch,:) = [V_T0, beta0 ,alpha0, phi0, theta0, psi0, p0,q0,r0,pn0, pe0, h0, q00, q10, q20, q30, Lat0, Lon0, Hgt0];

%start point



 a = 6378137.0;   % semi-major axis (metres)
    f = 1/298.2572; % flattening
    e2 = f * (2-f); % eccentricity squared
    e = sqrt(e2);   % first eccentricity

    Rmeridian = a * (1 - e2) / (sqrt(1 - e2 * sin(Lat0)^2))^3;

    %find the normal radius of curvature
    Rnormal = a / sqrt(1 - e2 * sin(Lat0)^2);

    Rn = Rnormal;
    Re = Rmeridian;


for N = 1:10

    
x0(startepoch,:) = [V_T0+N*2, beta0+N*0.5*pi/180 ,alpha0+N*0.5*pi/180, phi0+N*0.5*pi/180, theta0+N*0.5*pi/180, psi0+N*0.5*pi/1802, p0+N*0.03*pi/180,q0+N*0.03*pi/180,r0+N*0.03*pi/180,pn0, pe0, h0, q00, q10, q20, q30, Lat0+N*2/Rnormal, Lon0+N*2/Rnormal, Hgt0+N*2];    
    
 



for i = startepoch:100         


 
  delta_x_control(i,:) = [Elevator_truth(i), Aileron_truth(i), Rudder_truth(i), Throttle_truth(i)]; %this throttle needs to be scaled, not suree what it means at the moment, but probably needs to represent a change in thrust, since it is multiplied by the Xdelta_t term in the model    %control inputs will be the current control inputs at i
   
  %delta_x_control(i,:) = [0, 0,0,0]; %this throttle needs to be scaled, not suree what it means at the moment, but probably needs to represent a change in thrust, since it is multiplied by the Xdelta_t term in the model    %control inputs will be the current control inputs at i
   
   AeroCoeffsTruth = [CD_truth(i), CY_truth(i), CL_truth(i), Cl_truth(i), Cm_truth(i), Cn_truth(i)];
   EngCoeffsTruth = [MAP_truthEng(i), m_dotAir_truthEng(i),m_dotFuel_truthEng(i),BSFC_truthEng(i),P_truthEng(i), OMEGA_truthEng(i)] ;
   PropCoeffsTruth = [J_truth(i), CT_truth(i) ,CP_truth(i)];
   
   ForcesTruth = [Faero_truth(1,i),Faero_truth(2,i),Faero_truth(3,i),Fprop_truth(1,i),Fprop_truth(2,i),Fprop_truth(3,i),];
   
   MomentsTruth = [  Mcg_truth(1,i),Mcg_truth(2,i),Mcg_truth(3,i),Mprop_truth(1,i),Mprop_truth(2,i),Mprop_truth(3,i)];
   %MomentsTruth = [  Maero_truth(1,i),Maero_truth(2,i),Maero_truth(3,i),Mprop_truth(1,i),Mprop_truth(2,i),Mprop_truth(3,i)];
   
   
   MassTruth = Mass_truth(i);
    
   InertiaTruth = [Jx_truth(i), Jy_truth(i), Jz_truth(i), Jxz_truth(i)];
   CGposTruth = [CGxpos_truth(i), CGypos_truth(i), CGzpos_truth(i)];
   MachTruth =  Mach_truth(i);
    
  
% 
% TMatrixBodytoNED = T_Body2NED(Roll_truth(i),Pitch_truth(i), Yaw_truth(i));
% VWind_Truth = TMatrixBodytoNED*[u_wind_truth(i), v_wind_truth(i), w_wind_truth(i)]'; 
% 
% VWind_North_Truth(i) = VWind_Truth(1);
% VWind_East_Truth(i) = VWind_Truth(2);
% VWind_Down_Truth(i) =  VWind_Truth(3);
% 
% 
% Wind_NED = [VWind_North_Truth(i) VWind_East_Truth(i) VWind_Down_Truth(i)];
% 
% 
% Wind_Body = [u_wind_truth(i), v_wind_truth(i), w_wind_truth(i)]; 


Wind_Body = [0 0 0];


    %set alpha and beta to zero
%     x0(i,2) = 0;
%     x0(i,3) = 0;

%calculate state velocity dx


  %calculate MACH number for next time around
      M0 = x0(i,1)/AtmosTruth(4,i); 
      
      LatModel(i) = x0(i,17);
      LonModel(i) = x0(i,18);
      HgtModel(i) = x0(i,19);  
      
      
       [ECEFModel] = LLH2ECEF(LatModel(i),LonModel(i),HgtModel(i));

        %current position
        XposModel(i) = ECEFModel(1);
        YposModel(i) = ECEFModel(2);
        ZposModel(i) = ECEFModel(3);
      

  
     %integrate 
%     [k1dx,AeroCoeffs,Forces, Moments] = AerodynamicModelStevens(M0,GravityTruth(i),AtmosTruth(:,i),x0(i,:),delta_x_control(i,:), AeroCoeffsTruth, EngCoeffsTruth,PropCoeffsTruth,ForcesTruth, MomentsTruth, MassTruth,InertiaTruth,CGposTruth,MachTruth);
%     [k2dx,AeroCoeffs,Forces, Moments] = AerodynamicModelStevens(M0,GravityTruth(i),AtmosTruth(:,i),x0(i,:)+(ts/2)*(k1dx),delta_x_control(i,:), AeroCoeffsTruth, EngCoeffsTruth,PropCoeffsTruth,ForcesTruth, MomentsTruth, MassTruth,InertiaTruth,CGposTruth,MachTruth);
%     [k3dx,AeroCoeffs,Forces, Moments] = AerodynamicModelStevens(M0,GravityTruth(i),AtmosTruth(:,i),x0(i,:)+(ts/2)*(k2dx),delta_x_control(i,:), AeroCoeffsTruth, EngCoeffsTruth,PropCoeffsTruth,ForcesTruth, MomentsTruth, MassTruth,InertiaTruth,CGposTruth,MachTruth);
%     [k4dx,AeroCoeffs,Forces, Moments] = AerodynamicModelStevens(M0,GravityTruth(i),AtmosTruth(:,i),x0(i,:)+ts*(k3dx),delta_x_control(i,:), AeroCoeffsTruth, EngCoeffsTruth,PropCoeffsTruth,ForcesTruth, MomentsTruth, MassTruth,InertiaTruth,CGposTruth,MachTruth);
%     x0(i+1,:)  =  x0(i,:) + (ts/6)*(k1dx+2*k2dx+2*k3dx+k4dx) ;   
    
    
         AeroCoeffCorrections = zeros(1,6);
        Parameter_Noise = zeros(1,6);
        Wind_NED = [0,0,0];
        %GravityTruth(i) = 9.81;
    [dx, AeroCoeffs,Forces, Moments] = AerodynamicModelStevens(M0,GravityTruth(i),AtmosTruth(:,i),x0(i,:),delta_x_control(i,:), AeroCoeffsTruth, EngCoeffsTruth,PropCoeffsTruth,ForcesTruth, MomentsTruth, MassTruth,InertiaTruth,CGposTruth,MachTruth, Parameter_Noise, Wind_NED,AeroCoeffCorrections);

     
     
      x0(i+1,:) = x0(i,:) + ts*dx(1:19);
           
      
 
        %use quaternion truth
      %x0(i+1,13:16) = Quaternions_truth(:,i+1);
      

   
%get copy of state vector


     
 
     
     
     
     
     %set aoa and beta to zero
%      dx(2) = 0;
%      dx(3) = 0;
     
     
%       dx(7) = p_dottruth(i);
%      dx(8) = q_dottruth(i);
%      dx(9) = r_dottruth(i);
        
     
   
    
    
    
%     dx(1) = Airspeed_dot(i);
%      dx(2) = Sideslip_dot(i);
%       dx(3) = AOA_dot(i);
%     
%     
    
     %make dx the truth
     
      
     
%     dxtruth = [Airspeed_dot(i), Sideslip_dot(i), AOA_dot(i), Bank_dot(i), Pitch_dot(i) , Heading_dot(i),p_dot(i), q_dot(i), r_dot(i),dx(10), dx(11),dx(12),dx(13),dx(14),dx(15),dx(16)];
    
    
         
     %x0(i,:) =  x00+ 0.02*dx;
    
     
      
      %x0(i,:) =  x00+ 0.02*dx;      
      
      
     %save dx
%      
%      dx_save(:,i) = dx;
     
     
    % try passing the rates to my model and see if that fixes things?
     
         
         
     
     %set some of these parameters to the truth and see the effect on the
     %output
%      x0(i,1) = Airspeed(i);
%      x0(i,2) = Sideslip(i);
%      x0(i,3) = AOA(i);
%      
%      
%      x0(i,4) = Bank(i);
%      x0(i,5) = Pitch(i);
%      x0(i,6) = Heading(i);
%      
%      
%        x0(i,7) =  p_truth(i);
%         x0(i,8) = q_truth(i);
%         x0(i,9) = r_truth(i); 
%      
     
     
     
      
      
%    %  x0(i,1) = 30;
%     %wrap pitch , roll and yaw between min and max values    
%     x0(i,4:5) = wrap1(x0(i,4:5),-pi,pi);  %phi and theta, roll and pitch
%     x0(i,6) = wrap1(x0(i,6),0,2*pi);          %psi, yaw angle
%        
%     
%     
    

    
       
    
     
     %convert quaternions back to euler angles 
     
 %need to re-normalise the quaternions , after integrating them.    

q0 = x0(i+1,13);
q1 = x0(i+1,14);
q2 = x0(i+1,15);
q3 = x0(i+1,16);

length_q = sqrt(q0^2 + q1^2 + q2^2 + q3^2);

q0 = q0/length_q;
q1 = q1/length_q;
q2 = q2/length_q;
q3 = q3/length_q;

%copy them back into state vector
x0(i+1,13) = q0;
x0(i+1,14)= q1;
x0(i+1,15)= q2;
x0(i+1,16)= q3;


%should be equal to 1
blah = q0^2 + q1^2 + q2^2 + q3^2;


%     c11 = q0^2 + q1^2 - q2^2 - q3^2;
%     c12 = 2*(q1*q2 + q0*q3); 
%     c13 = 2*(q1*q3 - q0*q2);
%     c21 = 2*(q1*q2 - q0*q3);
%     c22 = q0^2 - q1^2 + q2^2 - q3^2;
%     c23 = 2*(q2*q3 + q0*q1);
%     c31 = 2*(q1*q3 + q0*q2);
%     c32 = 2*(q2*q3 - q0*q1);
%     c33 = q0^2 - q1^2 - q2^2 +q3^2;
% 
%  C_bn = [c11, c12, c13; c21, c22, c23; c31, c32, c33];


%check with Titterton, this Cbn matrix is correct

   
   
   %should be equal to 1
   %test
blah = q0^2 + q1^2 + q2^2 + q3^2;

quat = [q0, q1, q2, q3];
C_bn = GARDSim_DCMfromQuat(quat);

 %this is what I used to have
%   phi_quat = atan2(c32,c33);           %phi
%   
%   %stuck a - in front of these and it gives me the similar results as
%   %aerosim. note that in the model i must use Hgt_dot = h_dot;  not Hgt_dot
%   %= -h_dot; as well, for this to work.   
%   theta_quat = -asin(-c31);    %theta   
%   psi_quat = -atan2(c21,c11);   %psi   
%      
     
  
  
  %computer euler angles from quaternions using the method given in aerosim
  %p 148  - Euler angles from quaternions block. 
  

%   phi_quat = atan2( 2*(q0*q1 + q2*q3), q0^2 + q3^2 - q1^2 -q2^2);
%   theta_quat = asin(2*(q0*q2 - q1*q3));
%   psi_quat = atan2( 2*(q0*q3 + q1*q2), q0^2 + q1^2 - q2^2 - q3^2);
%   
  
  euler = QuatToEuler(quat);
  
  phi_quat = euler(1);
   theta_quat = euler(2);
    psi_quat = euler(3);
  
    x0(i+1,4) = phi_quat;
    x0(i+1,5) = theta_quat;
    x0(i+1,6) = psi_quat;   
    
    
    


%           
%     %wrap pitch , roll and yaw between min and max values    
%     x0(i,4:5) = wrap1(x0(i,4:5),-pi,pi);  %phi and theta, roll and pitch
%     x0(i,6) = wrap1(x0(i,6),0,2*pi);          %psi, yaw angle 
    
    
    

    
    
    
    
    %For Display
         
   
        
    %save aeromodel coefficients , forces and moments
    AeroCoeffsModel(i,:) = AeroCoeffs;  %CD, CY, CL, Cl, Cm, Cn
    
    ForcesModel(i,:) = Forces;
    MomentsModel(i,:) = Moments; 
   
   

   
%==============================================================
%update aircraft's position in ECEF based on the aerodynamic model output
%=================================================================

%convert u,v,w, velocities in body system (xb, yb, zb) into NED velocities 

% TMatrixBody2Wind = T_Body2Wind(alpha, beta);
% 
% %velocity in body axes
% V_b = TMatrixBody2Wind'*[V_T, 0, 0]';
% 
% u = V_b(1);
% v = V_b(2);
% w = V_b(3);
% 
% 
% BankOut = x0(i,4);
% PitchOut = x0(i,5);
% HeadingOut = x0(i,6); 
% 
% %PitchOut = Pitch(i);
% 
% 
% 
% 
% TMatrixBody2NED = T_Body2NED( BankOut,PitchOut,HeadingOut); %phi theta psi
% 
% 
% vel_NED = TMatrixBody2NED*[u v w]';
% 
% 
% %convert position in NED into ECEF
% 
% 
% [Lat,Lon,Hgt] = ECEF2LLH([Xpos(i), Ypos(i), Zpos(i)]);   %I assume you use epoch 8 and not at epoch 9 as coordinate origin for the NED. 
% 
% TMatrixECEF2NED = T_ECEF2NED(Lon,Lat);
% 
% ECEFModel = TMatrixECEF2NED'*[x0(i,10), x0(i,11), x0(i,12)]';




% 
% 
% 
% [LatTruth(i),LonTruth(i),HgtTruth(i)] = ECEF2LLH([Xpos(i), Ypos(i), Zpos(i)]);   %calculate longitude and latitude for truth
% 
% 
% x_LLH(1) = LatTruth(i);
% x_LLH(2) = LonTruth(i);
% x_LLH(3) = HgtTruth(i);

% %recalculate gravity and atmosphere for next time around..
% 
% [Gravity] = Earth_Gravity(x_LLH(1:3)',0.0818191908426^2,9.7803267714,0.00193185138639);
% 
% %=====================================================================
% %Atmospheric Parameters
% %=====================================================================
% 
% %Standard Atmosphere in SI units
% 
% %Altitude % ASsume GA aircraft operating betweeen 0 and 15,000 feet. 
% 
% %Temperature
% 
% 
% %[airDens,airPres,temp,soundSpeed] = Atmos(x_NED(3))
% 
% [P_atm T_atm p_atm a_atm] = Atmosphere_STDAtmosphere(x_LLH(3),287.0531,1.40,9.80665);
% 
% 
% %absolute ambient temperature in Kelvins
% T_amb = T_atm; %K
% R_gas = 287.0531,; %m^2/(K.s^2)
% 
% %Pressure
% 
% %Density
% rho = p_atm ; %kg/m^3 %density
% 
% Atmos = [T_amb R_gas rho];   %just use these parameters for now. 
% 
% 
% 
% a_atm = AtmosTruth(4,i);

% [LatModel(i),LonModel(i),HgtModel(i)] = ECEF2LLH([XposModel(i), YposModel(i), ZposModel(i)]); 
  
end  %end i

x0ALL(:,:,N) =  x0(:,:) ;

end %end N



%modify these if want to plot for a different range
startepochplot = startepoch;
endepochplot = 100;



for N = 1:10
    for i = startepoch:100
        
      LatModelPlot(i,N) = x0ALL(i,17,N);
      LonModelPlot(i,N) = x0ALL(i,18,N);
      HgtModelPlot(i,N) = x0ALL(i,19,N);  
      
      
       [ECEFModelPlot] = LLH2ECEF(LatModelPlot(i,N),LonModelPlot(i,N),HgtModelPlot(i,N));

        %current position
        XposModelPlot(i,N) = ECEFModelPlot(1);
        YposModelPlot(i,N) = ECEFModelPlot(2);
        ZposModelPlot(i,N) = ECEFModelPlot(3);        
    
        
    end
                        
        
XErrorECEFPlot(:,N) = Xpos_truth(startepochplot:endepochplot)' - XposModelPlot(startepochplot:endepochplot,N);
YErrorECEFPlot(:,N) = Ypos_truth(startepochplot:endepochplot)'  - YposModelPlot(startepochplot:endepochplot,N);
ZErrorECEFPlot(:,N) = Zpos_truth(startepochplot:endepochplot)'  - ZposModelPlot(startepochplot:endepochplot,N);

HeightErrorPlot(:,N) = Hgt_truth(startepochplot:endepochplot)' - HgtModelPlot(startepochplot:endepochplot,N); 
     
    
    
end



%plot3(XposModel(100:5000), YposModel(100:5000), ZposModel(100:5000))
%compare output ECEF from model with the true state.









%Plot ECEF position 


     subplot(2,1,1); plot(Xpos_truth(startepochplot:endepochplot)); title 'X truth (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
     subplot(2,1,2); plot(XposModelPlot(startepochplot:endepochplot,:)); title 'X model (m)'; xlabel('Time'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;


     subplot(2,1,1); plot(Ypos_truth(startepochplot:endepochplot)); title 'Y truth (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
     subplot(2,1,2); plot(YposModelPlot(startepochplot:endepochplot,:)); title 'Y model (m)'; xlabel('Time'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

     pause;
     clf;
     
     
     subplot(2,1,1); plot(Zpos_truth(startepochplot:endepochplot)); title 'Z truth (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
     subplot(2,1,2); plot(ZposModelPlot(startepochplot:endepochplot,:)); title 'Z model (m)'; xlabel('Time'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;
    
     
     
%plot ECEF position errors

plot(XErrorECEFPlot); title 'X Error (m)';xlabel('Time');
pause;

plot(YErrorECEFPlot); title 'Y Error (m)';xlabel('Time');
pause;

plot(ZErrorECEFPlot); title 'Z Error (m)';xlabel('Time');

pause;

%plot height

     subplot(2,1,1); plot(Hgt_truth(startepochplot:endepochplot)); title 'Height truth (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
     subplot(2,1,2); plot(HgtModelPlot(startepochplot:endepochplot,:)); title 'Height model (m)'; xlabel('Time'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);


pause;
clf



%plot height Error

plot(HeightErrorPlot); title 'Height Error (m)';xlabel('Time');

pause;

clf;




%compare outputs of model with truth. 

subplot(2,1,1); plot(Airspeed_truth(startepochplot:endepochplot)); title 'True Airspeed (V_T) Truth m/s'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,1)); title 'True Airspeed (V_T) Model m/s'; xlabel('Time');
pause;
clf;
AirspeedError = Airspeed_truth(startepochplot:endepochplot)' - x0(startepochplot:endepochplot,1);
plot(AirspeedError); title 'True Airspeed Error (m/s)';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(Beta_truth(startepochplot:endepochplot)*180/pi); title 'Beta Truth(theta) (deg)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,2)*180/pi); title 'Beta Model(theta) (deg)'; xlabel('Time');
pause;
clf;
SideslipError = Beta_truth(startepochplot:endepochplot)'*180/pi - x0(startepochplot:endepochplot,2)*180/pi;
plot(SideslipError); title 'Beta Error (deg)';xlabel('Time');
pause;
clf



subplot(2,1,1); plot(Alpha_truth(startepochplot:endepochplot)*180/pi); title 'Alpha Truth(theta) (deg)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,3)*180/pi); title 'Alpha Model(theta) (deg)'; xlabel('Time');
pause;
clf;
AOAError = Alpha_truth(startepochplot:endepochplot)'*180/pi - x0(startepochplot:endepochplot,3)*180/pi;
plot(AOAError); title 'Alpha Error (deg)';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(Roll_truth(startepochplot:endepochplot)*180/pi); title 'Roll Angle Truth(theta) (deg)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,4)*180/pi); title 'Roll Angle Model(theta) (deg)'; xlabel('Time');
pause;
clf;
BankError = Roll_truth(startepochplot:endepochplot)'*180/pi - x0(startepochplot:endepochplot,4)*180/pi;
plot(BankError); title 'BankError Error (deg)';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(Pitch_truth(startepochplot:endepochplot)*180/pi); title 'Pitch Angle Truth(theta) (deg)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,5)*180/pi); title 'Pitch Angle Model(theta) (deg)'; xlabel('Time');
pause;
clf;
PitchError = Pitch_truth(startepochplot:endepochplot)'*180/pi - x0(startepochplot:endepochplot,5)*180/pi;
plot(PitchError); title 'Pitch Error (deg)';xlabel('Time');
pause;
clf



for i = startepochplot: endepochplot

    Yaw_truth1Hz_wrap(i) = wrap1(Yaw_truth(i),-pi,pi);

end

YawErrorINS = Yaw_truth1Hz_wrap(startepochplot:endepochplot)' - (x0(startepochplot:endepochplot,6));

subplot(2,1,1); plot(Yaw_truth1Hz_wrap(startepochplot:endepochplot)*180/pi); title 'Yaw Angle Truth(theta) (deg)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,6)*180/pi); title 'Yaw Angle Model(theta) (deg)'; xlabel('Time');
pause;
clf;
HeadingError = YawErrorINS*180/pi;

%wrap heading error between +/-180 deg so the error can be seen clearly


plot(HeadingError); title 'Heading Error (deg)';xlabel('Time');
pause;
clf



subplot(2,1,1); plot(p_truth(startepochplot:endepochplot)*180/pi); title 'p Truth (deg/s)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,7)*180/pi); title 'p Model(deg/s)'; xlabel('Time');
pause;
clf;
pError = p_truth(startepochplot:endepochplot)'*180/pi - x0(startepochplot:endepochplot,7)*180/pi;
plot(pError); title 'p Error (deg/s)';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(q_truth(startepochplot:endepochplot)*180/pi); title 'q Truth (deg/s)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,8)*180/pi); title 'q Model (deg/s)'; xlabel('Time');
pause;
clf;
qError = q_truth(startepochplot:endepochplot)'*180/pi - x0(startepochplot:endepochplot,8)*180/pi;
plot(qError); title 'q Error (deg/s)';xlabel('Time');
pause;
clf

subplot(2,1,1); plot(r_truth(startepochplot:endepochplot)*180/pi); title 'r Truth (deg/s)'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,9)*180/pi); title 'r Model (deg/s)'; xlabel('Time');
pause;
clf;
rError = r_truth(startepochplot:endepochplot)'*180/pi - x0(startepochplot:endepochplot,9)*180/pi;
plot(rError); title 'r Error (deg/s)';xlabel('Time');
pause;
clf


%quaternions


subplot(2,1,1); plot(Quaternions_truth(1,startepochplot:endepochplot)); title 'q0 Truth'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,13)); title 'q0 Model'; xlabel('Time');
pause;
clf;
rError = Quaternions_truth(1,startepochplot:endepochplot)' - x0(startepochplot:endepochplot,13);
plot(rError); title 'q0 Error';xlabel('Time');
pause;
clf

subplot(2,1,1); plot(Quaternions_truth(2,startepochplot:endepochplot)); title 'q1 Truth'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,14)); title 'q0 Model'; xlabel('Time');
pause;
clf;
rError = Quaternions_truth(2,startepochplot:endepochplot)' - x0(startepochplot:endepochplot,14);
plot(rError); title 'q1 Error';xlabel('Time');
pause;
clf

subplot(2,1,1); plot(Quaternions_truth(3,startepochplot:endepochplot)); title 'q2 Truth'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,15)); title 'r Model (deg/s)'; xlabel('Time');
pause;
clf;
rError = Quaternions_truth(3,startepochplot:endepochplot)' - x0(startepochplot:endepochplot,15);
plot(rError); title 'q2 Error';xlabel('Time');
pause;
clf

subplot(2,1,1); plot(Quaternions_truth(4,startepochplot:endepochplot)); title 'q3 Truth'; xlabel('Time');
subplot(2,1,2); plot(x0(startepochplot:endepochplot,16)); title 'r Model (deg/s)'; xlabel('Time');
pause;
clf;
rError = Quaternions_truth(4,startepochplot:endepochplot)' - x0(startepochplot:endepochplot,16);
plot(rError); title 'q3 Error';xlabel('Time');
pause;
clf




checktruth = Quaternions_truth(1,:).^2 + Quaternions_truth(2,:).^2 + Quaternions_truth(3,:).^2 + Quaternions_truth(4,:).^2;
checkmodel = x0(startepochplot:endepochplot,13).^2 + x0(startepochplot:endepochplot,14).^2 + x0(startepochplot:endepochplot,15).^2 + x0(startepochplot:endepochplot,16).^2;


%Error = checktruth - checkmodel'; 

%the quaternions summed squared is not exactly equal to 1. 














%======================================================================
%Coefficients
%==========================

%note that this is the coefficients and forces and moments which WERE used
%to calculate the state velocities at time i. 

subplot(2,1,1); plot(CD_truth(startepochplot:endepochplot)); title 'CD Truth'; xlabel('Time');
subplot(2,1,2); plot(AeroCoeffsModel(startepochplot:endepochplot,1)); title 'CD Model'; xlabel('Time');
pause;
clf;
CDError = CD_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,1);
plot(CDError); title 'CD Error';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(CY_truth(startepochplot:endepochplot)); title 'CY Truth'; xlabel('Time');
subplot(2,1,2); plot(AeroCoeffsModel(startepochplot:endepochplot,2)); title 'CY Model'; xlabel('Time');
pause;
clf;
CYError = CY_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,2);
plot(CYError); title 'CY Error (m/s)';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(CL_truth(startepochplot:endepochplot)); title 'CL Truth'; xlabel('Time');
subplot(2,1,2); plot(AeroCoeffsModel(startepochplot:endepochplot,3)); title 'CL Model'; xlabel('Time');
pause;
clf;
CLError = CL_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,3);
plot(CLError); title 'CL Error';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(Cl_truth(startepochplot:endepochplot)); title 'Cl Truth'; xlabel('Time');
subplot(2,1,2); plot(AeroCoeffsModel(startepochplot:endepochplot,4)); title 'Cl Model'; xlabel('Time');
pause;
clf;
ClError = Cl_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,4);
plot(ClError); title 'Cl Error';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(Cm_truth(startepochplot:endepochplot)); title 'Cm Truth'; xlabel('Time');
subplot(2,1,2); plot(AeroCoeffsModel(startepochplot:endepochplot,5)); title 'Cm Model'; xlabel('Time');
pause;
clf;
CmError = Cm_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,5);
plot(CmError); title 'Cm Error';xlabel('Time');
pause;
clf


subplot(2,1,1); plot(Cn_truth(startepochplot:endepochplot)); title 'Cn Truth m/s'; xlabel('Time');
subplot(2,1,2); plot(AeroCoeffsModel(startepochplot:endepochplot,6)); title 'Cn Model'; xlabel('Time');
pause;
clf;
CnError = Cn_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,6);
plot(CnError); title 'Cn Error';xlabel('Time');
pause;
clf



%compare aerosim forces and moments and mine

%FORCES

subplot(2,1,1); plot(ForcesModel(startepochplot:endepochplot,1)); title 'Fxb Aero Aerosim '; xlabel('Time');
subplot(2,1,2); plot(ForcesModel(startepochplot:endepochplot,5)); title 'Fxb Aero Model'; xlabel('Time');
pause;
clf;
% CnError = Cn_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,6);
% plot(CnError); title 'Cn Error';xlabel('Time');
% pause;
% clf

subplot(2,1,1); plot(ForcesModel(startepochplot:endepochplot,2)); title 'Fy Aerosim '; xlabel('Time');
subplot(2,1,2); plot(ForcesModel(startepochplot:endepochplot,6)); title 'Fy Model'; xlabel('Time');
pause;
clf;

subplot(2,1,1); plot(ForcesModel(startepochplot:endepochplot,3)); title 'Fz Aerosim '; xlabel('Time');
subplot(2,1,2); plot(ForcesModel(startepochplot:endepochplot,7)); title 'Fz Model'; xlabel('Time');
pause;
clf;


subplot(2,1,1); plot(ForcesModel(startepochplot:endepochplot,4)); title 'Fprop Aerosim'; xlabel('Time');
subplot(2,1,2); plot(ForcesModel(startepochplot:endepochplot,8)); title 'Fprop Model'; xlabel('Time');
pause;
clf;

%MOMENTs

subplot(2,1,1); plot(MomentsModel(startepochplot:endepochplot,1)); title 'L Aerosim '; xlabel('Time');
subplot(2,1,2); plot(MomentsModel(startepochplot:endepochplot,5)); title 'L Model'; xlabel('Time');
pause;
clf;

subplot(2,1,1); plot(MomentsModel(startepochplot:endepochplot,2)); title 'M Aerosim '; xlabel('Time');
subplot(2,1,2); plot(MomentsModel(startepochplot:endepochplot,6)); title 'M Model'; xlabel('Time');
pause;
clf;


subplot(2,1,1); plot(MomentsModel(startepochplot:endepochplot,3)); title 'N Aerosim '; xlabel('Time');
subplot(2,1,2); plot(MomentsModel(startepochplot:endepochplot,7)); title 'N Model'; xlabel('Time');
pause;
clf;


subplot(2,1,1); plot(MomentsModel(startepochplot:endepochplot,4)); title 'Mprop Aerosim '; xlabel('Time');
subplot(2,1,2); plot(MomentsModel(startepochplot:endepochplot,8)); title 'Mprop Model'; xlabel('Time');
pause;
clf;











% 
% CnError = Cn_truth(startepochplot:endepochplot)' - AeroCoeffsModel(startepochplot:endepochplot,6);
% plot(CnError); title 'Cn Error';xlabel('Time');
% pause;
% clf























