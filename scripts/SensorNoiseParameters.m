%This is a script of most of the various noise parameters used in the
%simulation. 


global Sigma_PR_noise Beta_PR_noise 
global SVwithError TimeLow TimeHigh BiasError
global GRAS_OUTAGE_START_TIME GRAS_OUTAGE_STOP_TIME

global l_ax_INS_Error l_ay_INS_Error l_az_INS_Error l_p_INS_Error l_q_INS_Error l_r_INS_Error

global Sigma_Elevator_noise Beta_Elevator_noise Sigma_Aileron_noise Beta_Aileron_noise Sigma_Rudder_noise Beta_Rudder_noise Sigma_Throttle_noise Beta_Throttle_noise
 
global Sigma_CL_noise Beta_CL_noise Sigma_CY_noise Beta_CY_noise Sigma_CD_noise Beta_CD_noise
global Sigma_Cl_noise Beta_Cl_noise Sigma_Cm_noise Beta_Cm_noise Sigma_Cn_noise Beta_Cn_noise

global Sigma_Alpha_noise Beta_Alpha_noise Sigma_Beta_noise Beta_Beta_noise Sigma_Altimeter_noise Beta_Altimeter_noise

global Sigma_Temp_noise Beta_Temp_noise

global g_para1 g_para2 g_para3
 
    


%Use state 0 of the random number generator

%randn('state',0);


%*************
%GPS
%*************

%Receiver Noise

% Sigma_PR_noise = 7; %metres. %note this is all the other noises not including the clock bias , sat bias, earth rotation, atmospheric etc. 
% Beta_PR_noise = 5;

%putting a fault on the GPS

%bias parameters
% Type = 0   %type = 0 is constant bias, 1 = ramp, 3 = gauss markov process
% 
% Lower time limit
% Upper time limit
% 
% Satellite
SVwithError1 = 1; %SV number 
SVwithError2 = 3; 

TimeLow = 432000+10;
TimeHigh = 432000+100;
BiasError = 100;

%seconds into simulation when Gras outage occurs:
GRAS_OUTAGE_START_TIME = 0;

GRAS_OUTAGE_STOP_TIME = 0;

%*************
%INS
%*************

%Accel and gyro constant bias due to misalignment

% l_ax_INS_Error = 100e-6; % radians
% l_ay_INS_Error = 100e-6; % radians
% l_az_INS_Error = 100e-6; % radians
% 
% l_p_INS_Error = 100e-6; % radians
% l_q_INS_Error = 100e-6; % radians
% l_r_INS_Error = 100e-6; % radians


l_ax_INS_Error = 0; % radians
l_ay_INS_Error = 0; % radians
l_az_INS_Error = 0; % radians

l_p_INS_Error = 0; % radians
l_q_INS_Error = 0; % radians
l_r_INS_Error = 0; % radians


%white noise, corresponding to other things.. 



%error in knowledge of cofg location.

%*************
%Aero model
%*************
%Not sure what the right beta values are for these
%control inputs
% Sigma_Elevator_noise = 0*pi/180;  %radians
% Beta_Elevator_noise = 1/150; %1/ seconds
% 
% Sigma_Aileron_noise = 0*pi/180;   %radians
% Beta_Aileron_noise = 1/150; %seconds
% 
% Sigma_Rudder_noise = 0*pi/180;    %radians
% Beta_Rudder_noise = 1/150; %seconds
% 
% Sigma_Throttle_noise = 0;  %throttle setting in %
% Beta_Throttle_noise = 1/150; %seconds




Sigma_Elevator_noise = 0.02*pi/180;  %radians
Beta_Elevator_noise = 1/150e-8; %seconds

Sigma_Aileron_noise = 0.02*pi/180;   %radians
Beta_Aileron_noise = 1/150e-8; %seconds

Sigma_Rudder_noise = 0.02*pi/180;    %radians
Beta_Rudder_noise = 1/150e-8; %seconds

Sigma_Throttle_noise = 0.02;  %throttle setting in %
Beta_Throttle_noise = 1/150e-8; %seconds





%Propulsion

%mass

%Inertia

%Parameter Uncertainties all values in % of true value

Sigma_CL_noise = 0.068;   %this is to give roughly +/- 10% of true value , 0.068 is standard deviation of course for gaussian distribution
Beta_CL_noise = 1/120;

Sigma_CY_noise = 0.068;   %this is to give roughly +/- 10% of true value
Beta_CY_noise = 1/120;

Sigma_CD_noise = 0.068;   %this is to give roughly +/- 10% of true value
Beta_CD_noise = 1/120;

Sigma_Cl_noise = 0.068;   %this is to give roughly +/- 10% of true value
Beta_Cl_noise = 1/120;

Sigma_Cm_noise = 0.068;   %this is to give roughly +/- 10% of true value
Beta_Cm_noise = 1/120;

Sigma_Cn_noise = 0.068;   %this is to give roughly +/- 10% of true value
Beta_Cn_noise = 1/120;

% 
% Sigma_CL_noise = 0;   %this is to give roughly +/- 10% of true value
% Beta_CL_noise = 1/120;
% 
% Sigma_CY_noise = 0;   %this is to give roughly +/- 10% of true value
% Beta_CY_noise = 1/120;
% 
% Sigma_CD_noise = 0;   %this is to give roughly +/- 10% of true value
% Beta_CD_noise = 1/120;
% 
% Sigma_Cl_noise = 0;   %this is to give roughly +/- 10% of true value
% Beta_Cl_noise = 1/120;
% 
% Sigma_Cm_noise = 0;   %this is to give roughly +/- 10% of true value
% Beta_Cm_noise = 1/120;
% 
% Sigma_Cn_noise = 0;   %this is to give roughly +/- 10% of true value
% Beta_Cn_noise = 1/120;





%*************
%Air Data Sensors
%*************

% Angle of attack 
Sigma_Alpha_noise = 0.68*0.25*pi/180;  %this is 1 sigma .. %radians
Beta_Alpha_noise = 1/150e-8; %seconds; %assume its white

%-----------------------------------------------------
%Angle of Sideslip Sensor
%-----------------------------------------------------
  
Sigma_Beta_noise = 0.68*0.25*pi/180;% 3*pi/180;  %radians
Beta_Beta_noise = 1/150e-8; %seconds; %assume its white


% AOS
% 
% Airspeed indicator

  
Sigma_Airspeed_noise = 0.68*2.5;% 2;  %m/s
Beta_Airspeed_noise = 1/150e-8; %seconds; %assume its white




 

% 
% Magnetic Compass
% 
% Baro Altimeter

Sigma_Altimeter_noise = 34; %metres , roughly 100 feet. 
Beta_Altimeter_noise = 1/100; %seconds
    
%Temperature

%Random Noise randn is gaussian distributed
    
Sigma_Temp_noise = 2; %deg C
Beta_Temp_noise = 1/100; %seconds    
   

%*************
%Environmental 
%*************

% Gravity  %use gravity model , parameters

%[Gravity] = Earth_Gravity(x_LLH(1:3)',0.0818191908426^2,9.7803267714,0.00193185138639);

g_para1 = 0.0818191908426^2;
g_para2 = 9.7803267714;
g_para3 = 0.00193185138639;


%[Gravity] = Earth_Gravity(x_LLH(1:3)',g_para1,g_para2,g_para3);

% Pressure
% 
% 
% An

% 

% 
% %plot the noise and autocorrelation functions
% 
% 
% M=50;
% Rx=Rx_est(Temp_noise,M);  
% plot(Temp_noise)
% title('Gaussian-Markov Random Process')
% pause
% plot([-M:M],Rx)
% title('Autocorrelation function')
% 
% 
% %with matlab xcorr function need to take the envelope because its all noisy
% plot(xcorr(Temp_noise,Temp_noise)/1000); %need to divide by the length if using matlabs autocorrelation function






