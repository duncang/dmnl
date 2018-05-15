%this script loads the data for use in simulation

% 
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Flight22.2.07noJoystick\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Flight22.2.07noJoystick\pos_truth_llh.mat';
% end
% % 
% % %load ecef positions I calculated from lat lon and hgt
% % if ~exist('pos_truth_ecef')
% %     load 'data\Flight20.2.07\pos_truth_ecef.mat';
% % end
% 
% 
% 
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Flight22.2.07noJoystick\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Flight22.2.07noJoystick\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Flight22.2.07noJoystick\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Flight22.2.07noJoystick\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Flight22.2.07noJoystick\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Flight22.2.07noJoystick\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Flight22.2.07noJoystick\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Flight22.2.07noJoystick\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Flight22.2.07noJoystick\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Flight22.2.07noJoystick\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Flight22.2.07noJoystick\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Flight22.2.07noJoystick\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Flight22.2.07noJoystick\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Flight22.2.07noJoystick\InertiaCG.mat';
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
% UseNoisy = 1;
% 
% if UseNoisy == 1
%     load('data\Flight22.2.07noJoystick\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Flight22.2.07noJoystick\sensors_clean.mat');
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
%     load 'data\Flight22.2.07noJoystick\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Flight22.2.07noJoystick\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Flight22.2.07noJoystick\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Flight22.2.07noJoystick\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Flight22.2.07noJoystick\WindRates.mat';
% end

% 



%==================================================================
%==================================================================
%this one has no wind, use the one below for wind.
% 
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Flight22.2.07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Flight22.2.07\pos_truth_llh.mat';
% end
% % 
% % %load ecef positions I calculated from lat lon and hgt
% % if ~exist('pos_truth_ecef')
% %     load 'data\Flight20.2.07\pos_truth_ecef.mat';
% % end
% 
% 
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Flight22.2.07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Flight22.2.07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Flight22.2.07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Flight22.2.07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Flight22.2.07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Flight22.2.07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Flight22.2.07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Flight22.2.07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Flight22.2.07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Flight22.2.07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Flight22.2.07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Flight22.2.07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Flight22.2.07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Flight22.2.07\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% UseNoisy = 1;
% 
% if UseNoisy == 1
%     load('data\Flight22.2.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Flight22.2.07\sensors_clean.mat');
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
%     load 'data\Flight22.2.07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Flight22.2.07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Flight22.2.07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Flight22.2.07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Flight22.2.07\WindRates.mat';
% end




% 
% %==================================================================
% %==================================================================
% %this one has wind error and poorer INS
% 
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Flight27.3.07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Flight27.3.07\pos_truth_llh.mat';
% end
% % 
% % %load ecef positions I calculated from lat lon and hgt
% % if ~exist('pos_truth_ecef')
% %     load 'data\Flight20.2.07\pos_truth_ecef.mat';
% % end
% 
% 
% 
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Flight27.3.07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Flight27.3.07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Flight27.3.07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Flight27.3.07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Flight27.3.07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Flight27.3.07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Flight27.3.07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Flight27.3.07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Flight27.3.07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Flight27.3.07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Flight27.3.07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Flight27.3.07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Flight27.3.07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Flight27.3.07\InertiaCG.mat';
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
% UseNoisy = 1;
% 
% if UseNoisy == 1
%     load('data\Flight27.3.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Flight27.3.07\sensors_clean.mat');
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
%     load 'data\Flight27.3.07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Flight27.3.07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Flight27.3.07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Flight27.3.07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Flight27.3.07\WindRates.mat';
%end



% 
% 
%no wind, 100 deg/hr gyros, straight and level and two 180 degree turns, then straight and level sort of

%==================================================================
%==================================================================


% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Flight29.3.07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Flight29.3.07\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Flight29.3.07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Flight29.3.07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Flight29.3.07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Flight29.3.07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Flight29.3.07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Flight29.3.07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Flight29.3.07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Flight29.3.07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Flight29.3.07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Flight29.3.07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Flight29.3.07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Flight29.3.07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Flight29.3.07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Flight29.3.07\InertiaCG.mat';
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
%     load('data\Flight29.3.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Flight29.3.07\sensors_clean.mat');
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
%     load 'data\Flight29.3.07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Flight29.3.07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Flight29.3.07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Flight29.3.07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Flight29.3.07\WindRates.mat';
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %=============================================
% %sTRAIGHT IN dESCEND
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Holding Pattern test S&L 2.5.07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Holding Pattern test S&L 2.5.07\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Holding Pattern test S&L 2.5.07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Holding Pattern test S&L 2.5.07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Holding Pattern test S&L 2.5.07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Holding Pattern test S&L 2.5.07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Holding Pattern test S&L 2.5.07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Holding Pattern test S&L 2.5.07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Holding Pattern test S&L 2.5.07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Holding Pattern test S&L 2.5.07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Holding Pattern test S&L 2.5.07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Holding Pattern test S&L 2.5.07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Holding Pattern test S&L 2.5.07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Holding Pattern test S&L 2.5.07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Holding Pattern test S&L 2.5.07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Holding Pattern test S&L 2.5.07\InertiaCG.mat';
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
% UseNoisy = 1;
% 
% if UseNoisy == 1
%     load('data\Holding Pattern test S&L 2.5.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Holding Pattern test S&L 2.5.07\sensors_clean.mat');
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
%     load 'data\Holding Pattern test S&L 2.5.07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Holding Pattern test S&L 2.5.07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Holding Pattern test S&L 2.5.07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Holding Pattern test S&L 2.5.07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Holding Pattern test S&L 2.5.07\WindRates.mat';
% end
% 


%CONTINUOUS TURNING AROUND AND AROUND IN THE ONE SPOT


%=============================================


% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Holding Pattern test Turn 2.5.07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Holding Pattern test Turn 2.5.07\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Holding Pattern test Turn 2.5.07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Holding Pattern test Turn 2.5.07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Holding Pattern test Turn 2.5.07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Holding Pattern test Turn 2.5.07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Holding Pattern test Turn 2.5.07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Holding Pattern test Turn 2.5.07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Holding Pattern test Turn 2.5.07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Holding Pattern test Turn 2.5.07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Holding Pattern test Turn 2.5.07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Holding Pattern test Turn 2.5.07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Holding Pattern test Turn 2.5.07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Holding Pattern test Turn 2.5.07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Holding Pattern test Turn 2.5.07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Holding Pattern test Turn 2.5.07\InertiaCG.mat';
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
% UseNoisy = 1;
% 
% if UseNoisy == 1
%     load('data\Holding Pattern test Turn 2.5.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\Holding Pattern test Turn 2.5.07\sensors_clean.mat');
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
%     load 'data\Holding Pattern test Turn 2.5.07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Holding Pattern test Turn 2.5.07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Holding Pattern test Turn 2.5.07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Holding Pattern test Turn 2.5.07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Holding Pattern test Turn 2.5.07\WindRates.mat';
% end




% 
% %sTRAIGHT IN dESCEND
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
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% UseNoisy = 1;
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
% %RNA
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\RNAVGNNSYBBN\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\RNAVGNNSYBBN\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\RNAVGNNSYBBN\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\RNAVGNNSYBBN\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\RNAVGNNSYBBN\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\RNAVGNNSYBBN\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\RNAVGNNSYBBN\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\RNAVGNNSYBBN\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\RNAVGNNSYBBN\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\RNAVGNNSYBBN\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\RNAVGNNSYBBN\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\RNAVGNNSYBBN\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\RNAVGNNSYBBN\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\RNAVGNNSYBBN\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\RNAVGNNSYBBN\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\RNAVGNNSYBBN\InertiaCG.mat';
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
%     load('data\RNAVGNNSYBBN\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\RNAVGNNSYBBN\sensors_clean.mat');
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
%     load 'data\RNAVGNNSYBBN\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\RNAVGNNSYBBN\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\RNAVGNNSYBBN\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\RNAVGNNSYBBN\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\RNAVGNNSYBBN\WindRates.mat';
% end
% 
% 
% 
% 
% 
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
% UseNoisy = 1;
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


% %RNAV
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
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\TestFlight4.7.07\sensors_noisy.mat');
% load('data\TestFlight4.7.07\sensors_clean.mat');
% 
% UseNoisy = 1;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
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
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\100HzData1\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\100HzData1\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\100HzData1\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\100HzData1\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\100HzData1\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\100HzData1\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\100HzData1\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\100HzData1\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\100HzData1\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\100HzData1\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\100HzData1\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\100HzData1\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\100HzData1\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\100HzData1\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\100HzData1\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\100HzData1\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\100HzData1\sensors_noisy.mat');
% load('data\100HzData1\sensors_clean.mat');
% 
% UseNoisy = 1;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\100HzData1\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\100HzData1\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\100HzData1\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\100HzData1\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\100HzData1\WindRates.mat';
% end
% 





% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\100HzData2ForINSMODELFilter\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\100HzData2ForINSMODELFilter\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\100HzData2ForINSMODELFilter\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\100HzData2ForINSMODELFilter\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\100HzData2ForINSMODELFilter\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\100HzData2ForINSMODELFilter\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\100HzData2ForINSMODELFilter\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\100HzData2ForINSMODELFilter\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\100HzData2ForINSMODELFilter\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\100HzData2ForINSMODELFilter\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\100HzData2ForINSMODELFilter\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\100HzData2ForINSMODELFilter\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\100HzData2ForINSMODELFilter\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\100HzData2ForINSMODELFilter\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\100HzData2ForINSMODELFilter\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\100HzData2ForINSMODELFilter\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\100HzData2ForINSMODELFilter\sensors_noisy.mat');
% load('data\100HzData2ForINSMODELFilter\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\100HzData2ForINSMODELFilter\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\100HzData2ForINSMODELFilter\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\100HzData2ForINSMODELFilter\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\100HzData2ForINSMODELFilter\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\100HzData2ForINSMODELFilter\WindRates.mat';
% end



% 
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Flight22Nov07\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Flight22Nov07\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Flight22Nov07\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Flight22Nov07\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Flight22Nov07\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Flight22Nov07\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Flight22Nov07\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Flight22Nov07\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Flight22Nov07\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Flight22Nov07\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Flight22Nov07\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Flight22Nov07\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Flight22Nov07\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Flight22Nov07\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Flight22Nov07\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Flight22Nov07\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\Flight22Nov07\sensors_noisy.mat');
% load('data\Flight22Nov07\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\Flight22Nov07\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Flight22Nov07\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Flight22Nov07\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Flight22Nov07\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Flight22Nov07\WindRates.mat';
% end




% 
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\DebugAeroModel7.1.08\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\DebugAeroModel7.1.08\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\DebugAeroModel7.1.08\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\DebugAeroModel7.1.08\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\DebugAeroModel7.1.08\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\DebugAeroModel7.1.08\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\DebugAeroModel7.1.08\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\DebugAeroModel7.1.08\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\DebugAeroModel7.1.08\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\DebugAeroModel7.1.08\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\DebugAeroModel7.1.08\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\DebugAeroModel7.1.08\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\DebugAeroModel7.1.08\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\DebugAeroModel7.1.08\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\DebugAeroModel7.1.08\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\DebugAeroModel7.1.08\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\DebugAeroModel7.1.08\sensors_noisy.mat');
% load('data\DebugAeroModel7.1.08\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\DebugAeroModel7.1.08\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\DebugAeroModel7.1.08\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\DebugAeroModel7.1.08\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\DebugAeroModel7.1.08\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\DebugAeroModel7.1.08\WindRates.mat';
% end
% 
% 
% 

% % 
% % 
% % %THIS IS THE ONE I WAS USING BEFORE the Feb 01 Flight one at the bottom
% % 
% 
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\DebugAeromodelpqrchanged11.1.08\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\DebugAeromodelpqrchanged11.1.08\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\DebugAeromodelpqrchanged11.1.08\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\DebugAeromodelpqrchanged11.1.08\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\DebugAeromodelpqrchanged11.1.08\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\DebugAeromodelpqrchanged11.1.08\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\DebugAeromodelpqrchanged11.1.08\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\DebugAeromodelpqrchanged11.1.08\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\DebugAeromodelpqrchanged11.1.08\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\DebugAeromodelpqrchanged11.1.08\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\DebugAeromodelpqrchanged11.1.08\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\DebugAeromodelpqrchanged11.1.08\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\DebugAeromodelpqrchanged11.1.08\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\DebugAeromodelpqrchanged11.1.08\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\DebugAeromodelpqrchanged11.1.08\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\DebugAeromodelpqrchanged11.1.08\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\DebugAeromodelpqrchanged11.1.08\sensors_noisy.mat');
% load('data\DebugAeromodelpqrchanged11.1.08\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\DebugAeromodelpqrchanged11.1.08\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\DebugAeromodelpqrchanged11.1.08\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\DebugAeromodelpqrchanged11.1.08\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\DebugAeromodelpqrchanged11.1.08\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\DebugAeromodelpqrchanged11.1.08\WindRates.mat';
% end
% 
% 
% 







% 
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\EulerInt\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\EulerInt\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\EulerInt\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\EulerInt\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\EulerInt\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\EulerInt\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\EulerInt\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\EulerInt\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\EulerInt\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\EulerInt\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\EulerInt\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\EulerInt\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\EulerInt\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\EulerInt\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\EulerInt\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\EulerInt\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\EulerInt\sensors_noisy.mat');
% load('data\EulerInt\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\EulerInt\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\EulerInt\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\EulerInt\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\EulerInt\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\EulerInt\WindRates.mat';
% end

% 
% 
% 
% 
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\EulerInt2\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\EulerInt2\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\EulerInt2\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\EulerInt2\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\EulerInt2\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\EulerInt2\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\EulerInt2\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\EulerInt2\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\EulerInt2\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\EulerInt2\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\EulerInt2\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\EulerInt2\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\EulerInt2\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\EulerInt2\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\EulerInt2\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\EulerInt2\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\EulerInt2\sensors_noisy.mat');
% load('data\EulerInt2\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\EulerInt2\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\EulerInt2\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\EulerInt2\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\EulerInt2\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\EulerInt2\WindRates.mat';
% end








% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\RK4\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\RK4\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\RK4\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\RK4\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\RK4\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\RK4\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\RK4\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\RK4\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\RK4\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\RK4\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\RK4\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\RK4\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\RK4\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\RK4\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\RK4\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\RK4\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\RK4\sensors_noisy.mat');
% load('data\RK4\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\RK4\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\RK4\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\RK4\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\RK4\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\RK4\WindRates.mat';
% end
% 
% 






% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\RK42\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\RK42\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\RK42\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\RK42\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\RK42\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\RK42\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\RK42\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\RK42\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\RK42\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\RK42\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\RK42\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\RK42\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\RK42\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\RK42\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\RK42\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\RK42\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\RK42\sensors_noisy.mat');
% load('data\RK42\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\RK42\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\RK42\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\RK42\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\RK42\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\RK42\WindRates.mat';
% end
% 
% 
% 
% 
% 
% 

% 
% %Feb0108Flight
% 
% %aim of this flight was to see if i could realy ignore control inputs
% %during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% %turns and s&L flight and some climbs. 
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Feb0108Flight\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Feb0108Flight\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Feb0108Flight\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Feb0108Flight\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Feb0108Flight\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Feb0108Flight\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Feb0108Flight\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Feb0108Flight\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Feb0108Flight\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Feb0108Flight\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Feb0108Flight\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Feb0108Flight\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Feb0108Flight\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Feb0108Flight\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Feb0108Flight\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Feb0108Flight\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\Feb0108Flight\sensors_noisy.mat');
% load('data\Feb0108Flight\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\Feb0108Flight\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Feb0108Flight\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Feb0108Flight\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Feb0108Flight\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Feb0108Flight\WindRates.mat';
% end
% 





%Feb0108Flight
% 
% %aim of this flight was to see if i could realy ignore control inputs
% %during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% %turns and s&L flight and some climbs. 
% %RNAV
% 
% %POSITION
% 
% %load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\Feb0108FlightNoWind\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\Feb0108FlightNoWind\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\Feb0108FlightNoWind\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\Feb0108FlightNoWind\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\Feb0108FlightNoWind\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\Feb0108FlightNoWind\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\Feb0108FlightNoWind\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\Feb0108FlightNoWind\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\Feb0108FlightNoWind\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\Feb0108FlightNoWind\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\Feb0108FlightNoWind\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\Feb0108FlightNoWind\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\Feb0108FlightNoWind\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\Feb0108FlightNoWind\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\Feb0108FlightNoWind\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\Feb0108FlightNoWind\InertiaCG.mat';
% end
% 
% 
% 
% 
% %Load Inertial sensor data
% 
% % load IMU measurements data
% %this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\Feb0108FlightNoWind\sensors_noisy.mat');
% load('data\Feb0108FlightNoWind\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%     %clear 'sensors_noisy';
% else
%     
%     sensors = sensors_clean;
%     %clear 'sensors_clean';
% end
% %get copy of noisy sensor
% sensorsNoisy = sensors_noisy;
% sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\Feb0108FlightNoWind\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\Feb0108FlightNoWind\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\Feb0108FlightNoWind\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\Feb0108FlightNoWind\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\Feb0108FlightNoWind\WindRates.mat';
% end
% 


% 
% 
% % Journal Flight
% % 
% % aim of this flight was to see if i could realy ignore control inputs
% % during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% % turns and s&L flight and some climbs. 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\JournalFlight\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\JournalFlight\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\JournalFlight\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\JournalFlight\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\JournalFlight\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\JournalFlight\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\JournalFlight\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\JournalFlight\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\JournalFlight\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\JournalFlight\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\JournalFlight\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\JournalFlight\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\JournalFlight\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\JournalFlight\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\JournalFlight\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\JournalFlight\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\JournalFlight\sensors_noisy.mat');
% load('data\JournalFlight\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\JournalFlight\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\JournalFlight\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\JournalFlight\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\JournalFlight\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\JournalFlight\WindRates.mat';
% end








% % JoystickTurningFlight
% % 
% % aim of this flight was to see if the glitches in rnavapproach were there
% % when i did turns using the joystick
% %i also did some loops and rolls and hard manouvres 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\JoystickTurningFlight\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\JoystickTurningFlight\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\JoystickTurningFlight\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\JoystickTurningFlight\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\JoystickTurningFlight\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\JoystickTurningFlight\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\JoystickTurningFlight\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\JoystickTurningFlight\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\JoystickTurningFlight\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\JoystickTurningFlight\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\JoystickTurningFlight\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\JoystickTurningFlight\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\JoystickTurningFlight\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\JoystickTurningFlight\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\JoystickTurningFlight\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\JoystickTurningFlight\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\JoystickTurningFlight\sensors_noisy.mat');
% load('data\JoystickTurningFlight\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\JoystickTurningFlight\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\JoystickTurningFlight\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\JoystickTurningFlight\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\JoystickTurningFlight\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\JoystickTurningFlight\WindRates.mat';
% end
% 
% 


% 
% % JoystickTurningFlightODE4
% % 
% % aim of this flight was to see if the glitches in rnavapproach were there
% % when i did turns using the joystick
% %i also did some loops and rolls and hard manouvres 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\JoystickTurningFlightODE4\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\JoystickTurningFlightODE4\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\JoystickTurningFlightODE4\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\JoystickTurningFlightODE4\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\JoystickTurningFlightODE4\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\JoystickTurningFlightODE4\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\JoystickTurningFlightODE4\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\JoystickTurningFlightODE4\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\JoystickTurningFlightODE4\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\JoystickTurningFlightODE4\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\JoystickTurningFlightODE4\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\JoystickTurningFlightODE4\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\JoystickTurningFlightODE4\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\JoystickTurningFlightODE4\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\JoystickTurningFlightODE4\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\JoystickTurningFlightODE4\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\JoystickTurningFlightODE4\sensors_noisy.mat');
% load('data\JoystickTurningFlightODE4\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\JoystickTurningFlightODE4\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\JoystickTurningFlightODE4\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\JoystickTurningFlightODE4\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\JoystickTurningFlightODE4\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\JoystickTurningFlightODE4\WindRates.mat';
% end
% 
% 
% 
% 
% 
% 






%load Dunc's RNAV Approach Data


% 
% load 'data\rnav_approach\att_truth';
% 
% 
% load 'data\rnav_approach\pos_truth_ecef';
% 
% load 'data\rnav_approach\pos_truth_llh';
% 
% load 'data\rnav_approach\sensors_clean';
% 
% load 'data\rnav_approach\gravity';
% 
% 
% load 'data\rnav_approach\vel_truth';
% 
% 
% 
% 
% AtmosGrav
% 
% 
% 
% GravityTruth(i) =    AtmosGrav(6,i);
% 
% 
% data\DebugAeromodelpqrchanged11.1.08
% 
% 








% 
% 
% % Journal Flight
% % 
% % aim of this flight was to see if i could realy ignore control inputs
% % during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% % turns and s&L flight and some climbs. 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\RudderInputsFlight\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\RudderInputsFlight\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\RudderInputsFlight\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\RudderInputsFlight\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\RudderInputsFlight\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\RudderInputsFlight\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\RudderInputsFlight\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\RudderInputsFlight\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\RudderInputsFlight\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\RudderInputsFlight\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\RudderInputsFlight\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\RudderInputsFlight\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\RudderInputsFlight\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\RudderInputsFlight\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\RudderInputsFlight\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\RudderInputsFlight\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\RudderInputsFlight\sensors_noisy.mat');
% load('data\RudderInputsFlight\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\RudderInputsFlight\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\RudderInputsFlight\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\RudderInputsFlight\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\RudderInputsFlight\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\RudderInputsFlight\WindRates.mat';
% end
% 



% 
% % Journal Flight
% % 
% % aim of this flight was to see if i could realy ignore control inputs
% % during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% % turns and s&L flight and some climbs. 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\StraightandLevelFlight10Mins\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\StraightandLevelFlight10Mins\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\StraightandLevelFlight10Mins\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\StraightandLevelFlight10Mins\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\StraightandLevelFlight10Mins\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\StraightandLevelFlight10Mins\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\StraightandLevelFlight10Mins\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\StraightandLevelFlight10Mins\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\StraightandLevelFlight10Mins\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\StraightandLevelFlight10Mins\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\StraightandLevelFlight10Mins\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\StraightandLevelFlight10Mins\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\StraightandLevelFlight10Mins\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\StraightandLevelFlight10Mins\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\StraightandLevelFlight10Mins\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\StraightandLevelFlight10Mins\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\StraightandLevelFlight10Mins\sensors_noisy.mat');
% load('data\StraightandLevelFlight10Mins\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\StraightandLevelFlight10Mins\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\StraightandLevelFlight10Mins\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\StraightandLevelFlight10Mins\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\StraightandLevelFlight10Mins\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\StraightandLevelFlight10Mins\WindRates.mat';
% end


% 
% 
% 
% 
% % Journal Flight
% % 
% % aim of this flight was to see if i could realy ignore control inputs
% % during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% % turns and s&L flight and some climbs. 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\TurnsandRudder10mins\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\TurnsandRudder10mins\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\TurnsandRudder10mins\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\TurnsandRudder10mins\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\TurnsandRudder10mins\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\TurnsandRudder10mins\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\TurnsandRudder10mins\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\TurnsandRudder10mins\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\TurnsandRudder10mins\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\TurnsandRudder10mins\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\TurnsandRudder10mins\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\TurnsandRudder10mins\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\TurnsandRudder10mins\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\TurnsandRudder10mins\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\TurnsandRudder10mins\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\TurnsandRudder10mins\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\TurnsandRudder10mins\sensors_noisy.mat');
% load('data\TurnsandRudder10mins\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1
%     
%     sensors = sensors_noisy;
%     
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\TurnsandRudder10mins\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\TurnsandRudder10mins\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\TurnsandRudder10mins\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\TurnsandRudder10mins\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\TurnsandRudder10mins\WindRates.mat';
% end
% 



% % Journal Flight
% % 
% % aim of this flight was to see if i could realy ignore control inputs
% % during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% % turns and s&L flight and some climbs. 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\TurnsandRudder10minsnopqrEarth\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\TurnsandRudder10minsnopqrEarth\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\TurnsandRudder10minsnopqrEarth\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\TurnsandRudder10minsnopqrEarth\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\TurnsandRudder10minsnopqrEarth\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\TurnsandRudder10minsnopqrEarth\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\TurnsandRudder10minsnopqrEarth\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\TurnsandRudder10minsnopqrEarth\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\TurnsandRudder10minsnopqrEarth\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\TurnsandRudder10minsnopqrEarth\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\TurnsandRudder10minsnopqrEarth\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\TurnsandRudder10minsnopqrEarth\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\TurnsandRudder10minsnopqrEarth\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\TurnsandRudder10minsnopqrEarth\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\TurnsandRudder10minsnopqrEarth\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\TurnsandRudder10minsnopqrEarth\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\TurnsandRudder10minsnopqrEarth\sensors_noisy.mat');
% load('data\TurnsandRudder10minsnopqrEarth\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\TurnsandRudder10minsnopqrEarth\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\TurnsandRudder10minsnopqrEarth\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\TurnsandRudder10minsnopqrEarth\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\TurnsandRudder10minsnopqrEarth\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\TurnsandRudder10minsnopqrEarth\WindRates.mat';
% end
% 
% 
% 
% 
% 



% 
% 
% %navionwithFlapTest
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\navionwithFlapTest\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\navionwithFlapTest\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\navionwithFlapTest\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\navionwithFlapTest\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\navionwithFlapTest\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\navionwithFlapTest\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\navionwithFlapTest\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\navionwithFlapTest\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\navionwithFlapTest\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\navionwithFlapTest\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\navionwithFlapTest\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\navionwithFlapTest\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\navionwithFlapTest\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\navionwithFlapTest\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\navionwithFlapTest\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\navionwithFlapTest\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\navionwithFlapTest\sensors_noisy.mat');
% load('data\navionwithFlapTest\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\navionwithFlapTest\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\navionwithFlapTest\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\navionwithFlapTest\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\navionwithFlapTest\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\navionwithFlapTest\WindRates.mat';
% end
% 


% 
% % Journal Flight
% % 
% % aim of this flight was to see if i could realy ignore control inputs
% % during a descent phase (about 500ft/min), its mostly descent for first 60 seconds with some
% % turns and s&L flight and some climbs. 
% % RNAV
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\StraightWithRightthenLeftTurn\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\StraightWithRightthenLeftTurn\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\StraightWithRightthenLeftTurn\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\StraightWithRightthenLeftTurn\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\StraightWithRightthenLeftTurn\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\StraightWithRightthenLeftTurn\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\StraightWithRightthenLeftTurn\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\StraightWithRightthenLeftTurn\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\StraightWithRightthenLeftTurn\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\StraightWithRightthenLeftTurn\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\StraightWithRightthenLeftTurn\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\StraightWithRightthenLeftTurn\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\StraightWithRightthenLeftTurn\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\StraightWithRightthenLeftTurn\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\StraightWithRightthenLeftTurn\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\StraightWithRightthenLeftTurn\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\StraightWithRightthenLeftTurn\sensors_noisy.mat');
% load('data\StraightWithRightthenLeftTurn\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\StraightWithRightthenLeftTurn\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\StraightWithRightthenLeftTurn\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\StraightWithRightthenLeftTurn\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\StraightWithRightthenLeftTurn\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\StraightWithRightthenLeftTurn\WindRates.mat';
% end
% 
% 
% 









% % Test Flap 
% % 
% % aim of this flight was to see if the ADM was accurate with flap. Flap is
% % on the whole time, and no wind, and no earth rotation error in aerosim.
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\TestFlap\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\TestFlap\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\TestFlap\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\TestFlap\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\TestFlap\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\TestFlap\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\TestFlap\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\TestFlap\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\TestFlap\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\TestFlap\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\TestFlap\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\TestFlap\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\TestFlap\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\TestFlap\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\TestFlap\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\TestFlap\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\TestFlap\sensors_noisy.mat');
% load('data\TestFlap\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\TestFlap\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\TestFlap\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\TestFlap\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\TestFlap\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\TestFlap\WindRates.mat';
% end
% 

% 
% % Test Flap 
% % 
% % aim of this flight was to see if the ADM was accurate with flap. Flap is
% % on the whole time, and no wind, and no earth rotation error in aerosim.
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\TestFlapWithEarthRotation\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\TestFlapWithEarthRotation\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\TestFlapWithEarthRotation\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\TestFlapWithEarthRotation\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\TestFlapWithEarthRotation\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\TestFlapWithEarthRotation\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\TestFlapWithEarthRotation\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\TestFlapWithEarthRotation\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\TestFlapWithEarthRotation\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\TestFlapWithEarthRotation\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\TestFlapWithEarthRotation\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\TestFlapWithEarthRotation\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\TestFlapWithEarthRotation\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\TestFlapWithEarthRotation\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\TestFlapWithEarthRotation\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\TestFlapWithEarthRotation\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\TestFlapWithEarthRotation\sensors_noisy.mat');
% load('data\TestFlapWithEarthRotation\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\TestFlapWithEarthRotation\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\TestFlapWithEarthRotation\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\TestFlapWithEarthRotation\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\TestFlapWithEarthRotation\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\TestFlapWithEarthRotation\WindRates.mat';
% end




% % Test Flap 
% % 
% % aim of this flight was to see if the ADM was accurate with flap. Flap is
% % on the whole time, and no wind, and no earth rotation error in aerosim.
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\TestFlapWithEarthRotationAgain\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\TestFlapWithEarthRotationAgain\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\TestFlapWithEarthRotationAgain\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\TestFlapWithEarthRotationAgain\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\TestFlapWithEarthRotationAgain\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\TestFlapWithEarthRotationAgain\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\TestFlapWithEarthRotationAgain\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\TestFlapWithEarthRotationAgain\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\TestFlapWithEarthRotationAgain\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\TestFlapWithEarthRotationAgain\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\TestFlapWithEarthRotationAgain\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\TestFlapWithEarthRotationAgain\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\TestFlapWithEarthRotationAgain\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\TestFlapWithEarthRotationAgain\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\TestFlapWithEarthRotationAgain\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\TestFlapWithEarthRotationAgain\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\TestFlapWithEarthRotationAgain\sensors_noisy.mat');
% load('data\TestFlapWithEarthRotationAgain\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\TestFlapWithEarthRotationAgain\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\TestFlapWithEarthRotationAgain\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\TestFlapWithEarthRotationAgain\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\TestFlapWithEarthRotationAgain\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\TestFlapWithEarthRotationAgain\WindRates.mat';
% end
% 
% 







% Thesis Flight
% 


% POSITION
% 
% load ecef positions aerosim calculates
if ~exist('ECEF')
    load 'data/ThesisFlight/ECEF.mat';
end


if ~exist('pos_truth_llh')
    load 'data/ThesisFlight/pos_truth_llh.mat';
end


%VELOCITY
if ~exist('vel_truth')
    load 'data/ThesisFlight/vel_truth.mat';
end

if ~exist('vel_truthECEF')
    load 'data/ThesisFlight/vel_truthECEF.mat';
end




%ATTITUDE
%load attitude truth
if ~exist('att_truth')
    load 'data/ThesisFlight/att_truth.mat';
end


%load angular accelerations ( pdot qdot rdot)
if ~exist('AngAcc')
    load 'data/ThesisFlight/AngAcc.mat';
end

%load angular rates ( p q r)
if ~exist('AngRate')
    load 'data/ThesisFlight/AngRate.mat';
end


%AIR DATA
if ~exist('AirData1')
    load 'data/ThesisFlight/AirData1.mat';
end



%FLIGHT CONTROLS
if ~exist('flight_controls')
    load 'data/ThesisFlight/flight_controls.mat';
end


%ENGINE coefficients
%engine coefficients
if ~exist('EngCoeff')
    load 'data/ThesisFlight/EngCoeff.mat';
end

%load extra engine coefficients like shaft rotation speed (related
%to RPM)
if ~exist('EngParameters')
    load 'data/ThesisFlight/EngParameters.mat';
end

%load prop coefficients
if ~exist('PropCoeff')
    load 'data/ThesisFlight/PropCoeff.mat';
end

%load aero coefficients
if ~exist('AeroCoeff')
    load 'data/ThesisFlight/AeroCoeff.mat';
end




%LOAD AERODYNAMIC TRUTH

%load forces
if ~exist('Forces')
    load 'data/ThesisFlight/Forces.mat';
end

%load moments
if ~exist('Moments')
    load 'data/ThesisFlight/Moments.mat';
end


%load Inertia and CofG Truth
if ~exist('InertiaCG')
    load 'data/ThesisFlight/InertiaCG.mat';
end


% 
% 
% Load Inertial sensor data
% 
% load IMU measurements data
% this is from the 'sensors' output of the navion 6 DOF model

load('data/ThesisFlight/sensors_noisy.mat');
load('data/ThesisFlight/sensors_clean.mat');

UseNoisy = 0;

if UseNoisy == 1    
    
    sensors = sensors_noisy;   
   
else
    
    sensors = sensors_clean;
    
end
%get copy of noisy sensor
 sensorsNoisy = sensors_noisy;
 sensorsClean = sensors_clean;



%ATMOSPHERE TRUTH

%load Atmosphere and Gravity Truth
if ~exist('AtmosGrav')
    load 'data/ThesisFlight/AtmosGrav.mat';
end

%load Mach number Truth
if ~exist('Mach')
    load 'data/ThesisFlight/Mach.mat';
end

%load Quaternions Truth
if ~exist('Quaternions')
    load 'data/ThesisFlight/Quaternions.mat';
end

%Load Wind speed data

%load Wind speeds in body axes (m/s)
if ~exist('WindB')
    load 'data/ThesisFlight/WindB.mat';
end

%load Wind on angular rates (p q r) in rad/s
if ~exist('WindRates')
    load 'data/ThesisFlight/WindRates.mat';
end





% 
% % DescentShear
% % A descent , with wind shear, ie change in wind magnitude and direction
% 
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\DescentShear\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\DescentShear\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\DescentShear\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\DescentShear\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\DescentShear\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\DescentShear\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\DescentShear\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\DescentShear\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\DescentShear\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\DescentShear\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\DescentShear\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\DescentShear\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\DescentShear\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\DescentShear\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\DescentShear\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\DescentShear\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\DescentShear\sensors_noisy.mat');
% load('data\DescentShear\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\DescentShear\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\DescentShear\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\DescentShear\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\DescentShear\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\DescentShear\WindRates.mat';
% end




% 
% 
% % Thesis Flight
% % 
% 
% % 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data\VerticalShear\ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data\VerticalShear\pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data\VerticalShear\vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data\VerticalShear\vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data\VerticalShear\att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data\VerticalShear\AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data\VerticalShear\AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data\VerticalShear\AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data\VerticalShear\flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data\VerticalShear\EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data\VerticalShear\EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data\VerticalShear\PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data\VerticalShear\AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data\VerticalShear\Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data\VerticalShear\Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data\VerticalShear\InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data\VerticalShear\sensors_noisy.mat');
% load('data\VerticalShear\sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data\VerticalShear\AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data\VerticalShear\Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data\VerticalShear\Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data\VerticalShear\WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data\VerticalShear\WindRates.mat';
% end
% 
% 
% 
% 






% 
% % Thesis Flight
% % 
% 
% 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data/ThesisFlightFaster/ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data/ThesisFlightFaster/pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data/ThesisFlightFaster/vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data/ThesisFlightFaster/vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data/ThesisFlightFaster/att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data/ThesisFlightFaster/AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data/ThesisFlightFaster/AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data/ThesisFlightFaster/AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data/ThesisFlightFaster/flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data/ThesisFlightFaster/EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data/ThesisFlightFaster/EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data/ThesisFlightFaster/PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data/ThesisFlightFaster/AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data/ThesisFlightFaster/Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data/ThesisFlightFaster/Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data/ThesisFlightFaster/InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data/ThesisFlightFaster/sensors_noisy.mat');
% load('data/ThesisFlightFaster/sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data/ThesisFlightFaster/AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data/ThesisFlightFaster/Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data/ThesisFlightFaster/Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data/ThesisFlightFaster/WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data/ThesisFlightFaster/WindRates.mat';
% end
% 

% 
% 
% % Thesis Flight
% % 
% 
% 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data/ThesisFlightFast/ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data/ThesisFlightFast/pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data/ThesisFlightFast/vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data/ThesisFlightFast/vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data/ThesisFlightFast/att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data/ThesisFlightFast/AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data/ThesisFlightFast/AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data/ThesisFlightFast/AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data/ThesisFlightFast/flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data/ThesisFlightFast/EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data/ThesisFlightFast/EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data/ThesisFlightFast/PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data/ThesisFlightFast/AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data/ThesisFlightFast/Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data/ThesisFlightFast/Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data/ThesisFlightFast/InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data/ThesisFlightFast/sensors_noisy.mat');
% load('data/ThesisFlightFast/sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data/ThesisFlightFast/AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data/ThesisFlightFast/Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data/ThesisFlightFast/Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data/ThesisFlightFast/WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data/ThesisFlightFast/WindRates.mat';
% end
% 



%CRCSI data


% 
% % Thesis Flight
% % 
% 
% 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data/Aerosonde/ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data/Aerosonde/pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data/Aerosonde/vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data/Aerosonde/vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data/Aerosonde/att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data/Aerosonde/AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data/Aerosonde/AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data/Aerosonde/AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data/Aerosonde/flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data/Aerosonde/EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data/Aerosonde/EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data/Aerosonde/PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data/Aerosonde/AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data/Aerosonde/Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data/Aerosonde/Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data/Aerosonde/InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data/Aerosonde/sensors_noisy.mat');
% load('data/Aerosonde/sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data/Aerosonde/AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data/Aerosonde/Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data/Aerosonde/Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data/Aerosonde/WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data/Aerosonde/WindRates.mat';
% end
% 



% %for testing right turns
% 
% % Thesis Flight
% % 
% 
% 
% % POSITION
% % 
% % load ecef positions aerosim calculates
% if ~exist('ECEF')
%     load 'data/Testturn/ECEF.mat';
% end
% 
% 
% if ~exist('pos_truth_llh')
%     load 'data/Testturn/pos_truth_llh.mat';
% end
% 
% 
% %VELOCITY
% if ~exist('vel_truth')
%     load 'data/Testturn/vel_truth.mat';
% end
% 
% if ~exist('vel_truthECEF')
%     load 'data/Testturn/vel_truthECEF.mat';
% end
% 
% 
% 
% 
% %ATTITUDE
% %load attitude truth
% if ~exist('att_truth')
%     load 'data/Testturn/att_truth.mat';
% end
% 
% 
% %load angular accelerations ( pdot qdot rdot)
% if ~exist('AngAcc')
%     load 'data/Testturn/AngAcc.mat';
% end
% 
% %load angular rates ( p q r)
% if ~exist('AngRate')
%     load 'data/Testturn/AngRate.mat';
% end
% 
% 
% %AIR DATA
% if ~exist('AirData1')
%     load 'data/Testturn/AirData1.mat';
% end
% 
% 
% 
% %FLIGHT CONTROLS
% if ~exist('flight_controls')
%     load 'data/Testturn/flight_controls.mat';
% end
% 
% 
% %ENGINE coefficients
% %engine coefficients
% if ~exist('EngCoeff')
%     load 'data/Testturn/EngCoeff.mat';
% end
% 
% %load extra engine coefficients like shaft rotation speed (related
% %to RPM)
% if ~exist('EngParameters')
%     load 'data/Testturn/EngParameters.mat';
% end
% 
% %load prop coefficients
% if ~exist('PropCoeff')
%     load 'data/Testturn/PropCoeff.mat';
% end
% 
% %load aero coefficients
% if ~exist('AeroCoeff')
%     load 'data/Testturn/AeroCoeff.mat';
% end
% 
% 
% 
% 
% %LOAD AERODYNAMIC TRUTH
% 
% %load forces
% if ~exist('Forces')
%     load 'data/Testturn/Forces.mat';
% end
% 
% %load moments
% if ~exist('Moments')
%     load 'data/Testturn/Moments.mat';
% end
% 
% 
% %load Inertia and CofG Truth
% if ~exist('InertiaCG')
%     load 'data/Testturn/InertiaCG.mat';
% end
% 
% 
% % 
% % 
% % Load Inertial sensor data
% % 
% % load IMU measurements data
% % this is from the 'sensors' output of the navion 6 DOF model
% 
% load('data/Testturn/sensors_noisy.mat');
% load('data/Testturn/sensors_clean.mat');
% 
% UseNoisy = 0;
% 
% if UseNoisy == 1    
%     
%     sensors = sensors_noisy;   
%    
% else
%     
%     sensors = sensors_clean;
%     
% end
% %get copy of noisy sensor
%  sensorsNoisy = sensors_noisy;
%  sensorsClean = sensors_clean;
% 
% 
% 
% %ATMOSPHERE TRUTH
% 
% %load Atmosphere and Gravity Truth
% if ~exist('AtmosGrav')
%     load 'data/Testturn/AtmosGrav.mat';
% end
% 
% %load Mach number Truth
% if ~exist('Mach')
%     load 'data/Testturn/Mach.mat';
% end
% 
% %load Quaternions Truth
% if ~exist('Quaternions')
%     load 'data/Testturn/Quaternions.mat';
% end
% 
% %Load Wind speed data
% 
% %load Wind speeds in body axes (m/s)
% if ~exist('WindB')
%     load 'data/Testturn/WindB.mat';
% end
% 
% %load Wind on angular rates (p q r) in rad/s
% if ~exist('WindRates')
%     load 'data/Testturn/WindRates.mat';
% end
% 
% 
% 
% 




%========================================================================
%END OF LOADING DATA
%========================================================================